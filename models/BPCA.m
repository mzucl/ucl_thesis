classdef BPCA < handle
    properties 
        view            % ViewHandler
        K               % Number of latent dimensions/principal components

        % Model parameters with prior ds
        Z               % [K x N] GaussianDistributionContainer      [size: N; for each latent variable zn]
        
        mu              % [D x 1] GaussianDistribution               [D x 1; all observations have the same 'mu' parameter]
        
        W               % [D x K] GaussianDistributionContainer      [size: D; for each row in W matrix]
                        % Prior over W is defined per columns (each column
                        % has its own precision parameter, but update
                        % equations are defined by rows, so we are
                        % representing W as a size D container in a row
                        % format.

        alpha           % [K x 1] GammaDistributionContainer         [size: K]
        tau             % [scalar]GammaDistribution                  [scalar]

        % Optimization parameters
        maxIter
        tol

        %% Dependent properties
        %   They are not declared as dependent because they never change
        %   upon the initialization
        N   % Number of observations/latent variables
        D   % Dimensionality
    end


    properties (Constant)
        DEBUG = Constants.DEBUG;
    end

    
    methods
        function obj = BPCA(X, maxIter, tol)
            % Optional parameters: maxIter, tol
            if nargin < 1
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
            elseif nargin > 3
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too many arguments passed.']);
            end
            
            obj.view = ViewHandler(X, false);
            obj.N = obj.view.N;
            obj.D = obj.view.D;

            % BPCA can infer the right number of components, thus K is not
            % passed in as a parameter
            obj.K = obj.D - 1;

            % Set default values
            obj.maxIter = Constants.DEFAULT_MAX_ITER;
            obj.tol = Constants.DEFAULT_TOL;

            if nargin > 1
                obj.maxIter = maxIter;
                if nargin > 2
                    obj.tol = tol;
                end
            end

            %% Model setup and initialization
            % Z
            % ------------------------------------------------------ %
            % Initialize the model - set random values for the 'mu'
            % This means we will run the update equation for W first and
            % that we should set some values for all the moments that are
            % in those update equations.
            initZMu = randn(obj.K, 1);
            zPrior = GaussianDistribution(initZMu, eye(obj.K));
            obj.Z = GaussianDistributionContainer(obj.N, zPrior, true);            % cols = true

            % mu
            % empiricalMean = mean(obj.view.X, 2);
            muPrior = GaussianDistribution(0, 1/Constants.DEFAULT_GAUSS_PRECISION * eye(obj.D));
            obj.mu = GaussianDistribution(muPrior);

            % alpha
            alphaPrior = GammaDistribution(Constants.DEFAULT_GAMMA_A, Constants.DEFAULT_GAMMA_B);
            obj.alpha = GammaDistributionContainer(repmat(alphaPrior, obj.K, 1));
            

            % W 
            % Sample from obj.alpha for the prior
            %       Should we do this? The values for alpha are so small!!!
            % wPrior = GaussianDistribution(0, diag(1./obj.alpha.Val));

            % Use PPCA result as an initial point
            [W_PPCA, sigmaSq] = PPCA(obj.view.X', obj.D - 1);
            
            wPriors = repmat(GaussianDistribution(), obj.D, 1); % Preallocate
            wPrior = GaussianDistribution(randn(obj.K, 1), 1e-3 * eye(obj.K));
            for d = 1:obj.D
                wPriors(d) = GaussianDistribution(W_PPCA(d, :), 1e-3 * eye(obj.K));
            end

            obj.W = GaussianDistributionContainer(obj.D, wPrior, false);    % cols = false
            

            % tau
            tauPrior = GammaDistribution(Constants.DEFAULT_GAMMA_A, Constants.DEFAULT_GAMMA_B);
            obj.tau = GammaDistribution(tauPrior);

            % Model initialization - second part
            % The first update equation is for W, so we need to initialize
            % everything that is used in those two equations and those
            % initilizations are given below.
            %   obj.tau.expInit
            %   obj.alpha.expCInit
            %   obj.mu.expInit
            % ----------------------------------------------------------------
            obj.tau.setExpInit(1 / sigmaSq);
            obj.alpha.setExpCInit(repmat(1e-3, obj.K, 1));
            obj.mu.setExpInit(randn(obj.D, 1));

            obj.qConstantUpdates();
        end


        
        %% Update methods
        function obj = qConstantUpdates(obj)
            % alpha.a
            obj.alpha.updateAllDistributionsA(obj.alpha.ds(1).prior.a + obj.D/2);
            % tau.a
            obj.tau.updateA(obj.tau.prior.a + (obj.N * obj.D)/2);
        end


        function obj = qZUpdate(obj, it)
            tauExp = Utility.ternary(it == 1, obj.tau.getExpInit(), obj.tau.E);
            muExp = Utility.ternary(it == 1, obj.mu.getExpInit(), obj.mu.E);

            covNew = Utility.matrixInverse(eye(obj.K) + tauExp * obj.W.E_CtC);
            obj.Z.updateAllDistributionsCovariance(covNew);
    
            obj.Z.updateAllDistributionsMu(tauExp * covNew * obj.W.E_Ct * ...
                (obj.view.X - muExp));
        end

        
        function obj = qWUpdate(obj, it)
            alphaExp = Utility.ternary(it == 1, obj.alpha.getExpCInit(), obj.alpha.EC);
            tauExp = Utility.ternary(it == 1, obj.tau.getExpInit(), obj.tau.E);
            muExp = Utility.ternary(it == 1, obj.mu.getExpInit(), obj.mu.E);

            covNew = Utility.choleskyInverse(diag(alphaExp) + tauExp * obj.Z.E_CCt);
    
            obj.W.updateAllDistributionsCovariance(covNew);

            obj.W.updateAllDistributionsMu(tauExp * covNew * obj.Z.EC * ...
                (obj.view.X' - muExp'));
        end


        function obj = qAlphaUpdate(obj)
            obj.alpha.updateAllDistributionsB(obj.alpha.ds(1).prior.b + ...
                1/2 * obj.W.getExpectationOfColumnsNormSq());
        end


        function obj = qMuUpdate(obj)
            tauExp = obj.tau.E;

            newCov = (1/(obj.mu.PPrec + obj.N * tauExp)) * eye(obj.D);
            newMu = tauExp * newCov * (obj.view.X - obj.W.EC * obj.Z.EC) * ones(obj.N, 1);

            obj.mu.updateParameters(newMu, newCov);

            % obj.mu.clearCache();
        end


        function obj = qTauUpdate(obj)
            expWtW_tr = obj.W.E_CtC';
            expZZt = obj.Z.E_CCt;
            expWZ = obj.W.EC * obj.Z.EC;
            
            newBVal = obj.tau.prior.b + 1/2 * obj.view.Tr_XtX + obj.N/2 * obj.mu.E_XtX + ...
                1/2 * dot(expWtW_tr(:), expZZt(:)) - dot(obj.view.X(:), expWZ(:)) + ...
                obj.mu.E' * (expWZ - obj.view.X) * ones(obj.N, 1);
        
            obj.tau.updateB(newBVal);
        end
        
        

        %% fit() and ELBO
        function [elboVals, it, resArr] = fit(obj, elboIterStep)
            if nargin < 2
                elboIterStep = 1;
            end
            elboVals = -Inf(1, obj.maxIter);
            
            for it = 1:obj.maxIter
                obj.qZUpdate(it);
                obj.qWUpdate(it);
                obj.qMuUpdate();
                obj.qAlphaUpdate();
                obj.qTauUpdate();

                if mod(it, elboIterStep) ~= 0
                    continue;
                end

                % Compute elbo
                if obj.DEBUG
                    resArr = cell(1, obj.maxIter);
                    [currElbo, res] = obj.computeELBO();
                    resArr{it} = res;

                    if it ~= 1
                        disp(['======= ELBO increased by: ', num2str(currElbo - elboVals(it - 1))]);
                    end
                else
                    currElbo = obj.computeELBO();
                end

                % ELBO has to increase from iteration to iteration
                if it ~= 1 && currElbo < elboVals(it - 1)
                    fprintf(2, 'ELBO decreased in iteration %d\n!!!', it);
                end 

                elboVals(it) = currElbo;

                % Check for convergence
                if it ~= 1 && abs(currElbo - elboVals(it - 1)) / abs(currElbo) < obj.tol
                    disp(['Convergence at iteration: ', num2str(it)]);
                    elboVals = elboVals(1:it); % cut the -Inf values at the end
                    if obj.DEBUG
                        resArr = resArr(1:it);
                    end
                    break;
                end
            end
        end
        

        function [elbo, res] = computeELBO(obj)
            elbo = obj.getExpectationLnPX() + obj.Z.E_LnPC + obj.getExpectationLnPW() + ... % p(.)
                obj.alpha.E_LnPC + obj.mu.E_LnP + obj.tau.E_LnP + ... % p(.)
                obj.Z.HC + obj.W.HC + obj.alpha.HC + obj.mu.H + obj.tau.H; % q(.)
            if obj.DEBUG
                res = {};
                res.pX = obj.getExpectationLnPX();
                res.pZ = obj.Z.E_LnPC;
                res.pW = obj.getExpectationLnPW();
                res.pAlpha = obj.alpha.E_LnPC;
                res.pMu = obj.mu.E_LnP;
                res.pTau = obj.tau.E_LnP;
                
                res.qZ = obj.Z.HC;
                res.qW = obj.W.HC;
                res.qAlpha = obj.alpha.HC;
                res.qMu = obj.mu.H;
                res.qTau = obj.tau.H;

                res.elbo = elbo;
            end
        end


        function value = getExpectationLnPX(obj)
            % Setup
            expWtW_tr = obj.W.E_CtC';
            expZZt = obj.Z.E_CCt;
            expW_tr = obj.W.EC';
            expZXt = obj.Z.EC * obj.view.X';

            value = obj.N * obj.D/2 * (obj.tau.E_LnX - log(2 * pi)) - obj.tau.E/2 * ( ...
                obj.view.Tr_XtX - 2 * dot(expW_tr(:), expZXt(:)) ...
                - 2 * obj.mu.E_Xt * (obj.view.X * ones(obj.N, 1)) ...
                + 2 * (expW_tr * obj.mu.E)' * (obj.Z.EC * ones(obj.N, 1)) ...
                + dot(expWtW_tr(:), expZZt(:)) + obj.N * obj.mu.E_XtX);
        end


        function value = getExpectationLnPW(obj)
            value = -1/2 * obj.W.getExpectationOfColumnsNormSq()' * obj.alpha.EC + ...
                + obj.D/2 * (obj.alpha.E_LnC - obj.K * log(2*pi));
        end
    end
end