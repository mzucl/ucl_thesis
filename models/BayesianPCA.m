classdef BayesianPCA < handle
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
    end

    properties (Dependent)
        N   % Number of observations/latent variables
        D   % Dimensionality
    end
    
    methods
        function obj = BayesianPCA(X, maxIter, tol)
            % Optional parameters: maxIter, tol
            if nargin < 1
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
            elseif nargin > 3
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too many arguments passed.']);
            end
            
            % Set obj.X right away, so it can be used below to set obj.K
            obj.view = ViewHandler(X, false);

            % BPCA can infer the right number of components
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
        end


        
        %% Update methods
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

        function obj = qAlphaUpdate(obj, it)
            % Alpha is updated to the same value through the iterations, so
            % it is enough to update it once
            if it == 1
                newAVal = obj.alpha.ds(1).prior.a + obj.D/2;
                obj.alpha.updateAllDistributionsA(newAVal);
            end
            
            newBVals = obj.alpha.ds(1).prior.b + 1/2 * obj.W.getExpectationOfColumnsNormSq();

            obj.alpha.updateAllDistributionsB(newBVals);
        end

        function obj = qMuUpdate(obj)
            tauExp = obj.tau.E;

            newCov = (1/(obj.mu.PPrec + obj.N * tauExp)) * eye(obj.D);
            newMu = tauExp * newCov * (obj.view.X - obj.W.EC * obj.Z.EC) * ones(obj.N, 1);

            obj.mu.updateParameters(newMu, newCov);
        end

        % obj.tau is GammaDistribution
        function obj = qTauUpdate(obj, it)
            % Alpha is updated to the same value through the iterations, so
            % it is enough to update it once
            if it == 1
                newAVal = obj.tau.prior.a + (obj.N * obj.D)/2;
                obj.tau.updateA(newAVal);
            end

            expWtW_Tr = obj.W.E_CtC';
            expZZt = obj.Z.E_CCt;
            expWZ = obj.W.EC * obj.Z.EC;
            
            newBVal = obj.tau.prior.b + 1/2 * obj.view.Tr_XtX + obj.N/2 * obj.mu.E_XtX + ...
                1/2 * dot(expWtW_Tr(:), expZZt(:)) - dot(obj.view.X(:), expWZ(:)) + ...
                obj.mu.E' * (expWZ - obj.view.X) * ones(obj.N, 1);
        
            obj.tau.updateB(newBVal);
        end
        
        

        %% fit() and ELBO
        function [elboVals, it, resArr] = fit(obj)
            elboVals = -Inf(1, obj.maxIter);
            resArr = cell(1, obj.maxIter);
        
            for it = 1:obj.maxIter
                obj.qZUpdate(it);
                obj.qWUpdate(it);
                obj.qMuUpdate();
                obj.qAlphaUpdate(it);
                obj.qTauUpdate(it);

                [currElbo, res] = obj.computeELBO();

                resArr{it} = res;

                % CHECK: ELBO has to increase from iteration to iteration
                if it ~= 1 && currElbo < elboVals(it - 1)
                    fprintf(2, 'ELBO decreased in iteration %d\n', it);
                end 

                elboVals(it) = currElbo;

                if it ~= 1
                    disp(['======= ELBO increased by: ', num2str(currElbo - elboVals(it - 1))]);
                end

                % Check for convergence
                if it ~= 1 && abs(currElbo - elboVals(it - 1)) / abs(currElbo) < obj.tol
                    disp(['Convergence at iteration: ', num2str(it)]);
                    elboVals = elboVals(1:it); % cut the -Inf values at the end
                    resArr = resArr(1:it);
                    break;
                end
            end
        end
        
        function [elbo, res] = computeELBO(obj)
            % DEBUG
            res = {};
            res.pX = obj.getExpectationLnPX();
            res.pZ = obj.Z.E_LnPC;
            res.pW = obj.getExpectationLnW();
            res.pAlpha = obj.alpha.E_LnPC;
            res.pMu = obj.mu.E_LnP;
            res.pTau = obj.tau.E_LnP;
            
            res.qZ = obj.Z.HC;
            res.qW = obj.W.HC;
            res.qAlpha = obj.alpha.HC;
            res.qMu = obj.mu.H;
            res.qTau = obj.tau.H;
            % DEBUG

            elbo = obj.getExpectationLnPX() + obj.Z.E_LnPC + obj.getExpectationLnW() + ... % p(.)
                obj.alpha.E_LnPC + obj.mu.E_LnP + obj.tau.E_LnP + ... % p(.)
                obj.Z.HC + obj.W.HC + obj.alpha.HC + obj.mu.H + obj.tau.H; % q(.)

            % DEBUG
            res.elbo = elbo;
            % DEBUG
        end

        function value = getExpectationLnPX(obj)
            value = obj.N * obj.D/2 * (obj.tau.E_Ln - log(2 * pi)) - obj.tau.E/2 * ( ...
                obj.view.Tr_XtX - 2 * trace(obj.W.EC * obj.Z.EC * obj.view.X') ...
                - 2 * obj.mu.E_Xt * obj.view.X * ones(obj.N, 1) ...
                + 2 * obj.mu.E_Xt * obj.W.EC * obj.Z.EC * ones(obj.N, 1) ...
                + trace(obj.W.E_CtC * obj.Z.E_CCt) + obj.N * obj.mu.E_XtX);

        end

        function value = getExpectationLnW(obj)
            value = -1/2 * obj.W.getExpectationOfColumnsNormSq()' * obj.alpha.EC + ...
                + obj.D/2 * (obj.alpha.E_LnC - obj.K * log(2*pi));
        end

        %% Getters
        function value = get.N(obj)
            value = obj.view.N;
        end
        
        function value = get.D(obj)
            value = obj.view.D;
        end
    end
end
