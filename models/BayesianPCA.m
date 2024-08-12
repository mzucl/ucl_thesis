classdef BayesianPCA < handle
    properties 
        X               % ViewHandler
        K               % Number of latent dimensions/principal components

        % Model parameters with prior distributions
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
        function obj = BayesianPCA(X, K, maxIter, tol)
            if nargin < 1
                error(['Error in class ' class(obj) ': Too few arguments passed.']);
            end
            
            % Set obj.X right away, so it can be used below to set obj.K
            obj.X = ViewHandler(X);
            % totalVar = sum(diag(cov(X')));

            % Set parameters to default values that will be updated if value is
            % provided
            obj.maxIter = Constants.DEFAULT_MAX_ITER;
            obj.tol = Constants.DEFAULT_TOL;

            % Optional parameters: maxIter, tol
            switch nargin
                case 1
                    obj.K = obj.D - 1;
                case 2
                    obj.K = K;
                case 3
                    obj.K = K;
                    obj.maxIter = maxIter;
                case 4
                    obj.K = K;
                    obj.maxIter = maxIter;
                    obj.tol = tol;
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
            obj.Z = GaussianDistributionContainer(obj.N, zPrior, true);

            % mu
            muPrior = GaussianDistribution(0, 1/Constants.DEFAULT_GAUSS_PRECISION * eye(obj.D));
            obj.mu = GaussianDistribution(muPrior);

            % alpha
            alphaPrior = GammaDistribution(Constants.DEFAULT_GAMMA_A, Constants.DEFAULT_GAMMA_B);
            obj.alpha = GammaDistributionContainer(repmat(alphaPrior, K, 1));
            

            % W; sample from obj.alpha for the prior
            %       Should we do this? The values for alpha are so small!!!
            % wPrior = GaussianDistribution(0, diag(1./obj.alpha.Value));
            wPrior = GaussianDistribution(0, eye(K));
            obj.W = GaussianDistributionContainer(obj.D, wPrior, false);
            
            
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
            obj.tau.setExpInit(1e-3);
                        
            obj.alpha.setExpCInit(repmat(1e-1, obj.K, 1));
            % The line below is from GFA code
            % obj.alpha.setExpCInit(repmat(obj.K * obj.D / (totalVar - 1 / obj.tau.getExpInit()), obj.K, 1));
            obj.mu.setExpInit(randn(obj.D, 1));
        end


        
        %% Update methods
        % obj.Z is GaussianDistributionContainer(cols = true)
        function obj = qZUpdate(obj)
            % Update mu
            for n = 1:obj.N
                % All latent variables have the same covariance
                muNew = obj.tau.Expectation * obj.Z.distributions(n).cov * obj.W.ExpectationCt * ...
                    (obj.X.getObservation(n) - obj.mu.Expectation);
                obj.Z.updateDistributionMu(n, muNew);
            end

            % Update covariance
            covNew = Utility.matrixInverse(eye(obj.K) + obj.tau.Expectation * ... 
                obj.W.ExpectationCtC);
            obj.Z.updateAllDistributionsCovariance(covNew);
        end
        
        % obj.W is GaussianDistributionContainer(cols = false)
        % [NOTE] Initial distribution is defined per columns of W, but
            % update equations are defined per rows
        function obj = qWUpdate(obj, it)
            disp(['Min W value: ', num2str(min(obj.W.ExpectationC, [], 'all'))]);
            disp(['Max W value: ', num2str(max(obj.W.ExpectationC, [], 'all'))]);
            
            % In the first iteration we perform the update based on the
            % initialized moments of tau, mu and alpha, and in every
            % subsequent iteration we use the 'normal' update equations

            % TODO (medium): DRY this code: define the vars for the values
            % that are different based on the value of 'it' before the
            % update equations. Use Utility.ternary(it == 1, ...)
            %   obj.alpha.expCInit
            %   obj.tau.expInit
            %   obj.mu.expInit
            if it > 1
                % Covariance update
                covNew = Utility.matrixInverse((diag(obj.alpha.ExpectationC)) + ...
                    obj.tau.Expectation * obj.Z.ExpectationCCt);
    
                obj.W.updateAllDistributionsCovariance(covNew);
    
                % Mean update
                for d = 1:obj.D
                    muNew = zeros(obj.K, 1);
                    for n = 1:obj.N
                        muNew = muNew + obj.Z.Expectation{n} * (obj.X.getObservationEntry(n, d) - obj.mu.Expectation(d));
                    end
                    muNew = obj.tau.Expectation * obj.W.distributions(d).cov * muNew;
                    obj.W.updateDistributionMu(d, muNew);
                end   
            % First iteration - use 'expInit' instead of read expectations
            else
                % Covariance update
                covNew = Utility.matrixInverse((diag(obj.alpha.getExpCInit())) + ...
                    obj.tau.getExpInit() * obj.Z.ExpectationCCt);
    
                obj.W.updateAllDistributionsCovariance(covNew);
    
                % Mean update
                for d = 1:obj.D
                    muNew = zeros(obj.K, 1);
                    for n = 1:obj.N
                        muExpInit = obj.mu.getExpInit();
                        muNew = muNew + obj.Z.Expectation{n} * (obj.X.getObservationEntry(n, d) - muExpInit(d));
                    end
                    muNew = obj.tau.getExpInit() * obj.W.distributions(d).cov * muNew;
                    obj.W.updateDistributionMu(d, muNew);
                end   
            end
        end

        % obj.alpha is GammaDistributionContainer
        function obj = qAlphaUpdate(obj)
            newAVal = obj.alpha.distributions(1).prior.a + obj.D/2; % All 'a' values are the same
            newBVals = obj.alpha.distributions(1).prior.b + 1/2 * obj.W.getExpectationOfColumnsNormSq();

            obj.alpha.updateAllDistributionsParams(newAVal, newBVals);
        end

        % obj.mu is GaussianDistribution
        function obj = qMuUpdate(obj)
            newMu = zeros(obj.D, 1);
           
            for n = 1:obj.N
                newMu = newMu + obj.X.getObservation(n) - obj.W.ExpectationC * obj.Z.Expectation{n};
            end
            newMu = obj.tau.Expectation * obj.mu.cov * newMu;
            newCov = 1/(obj.mu.PriorPrecision + obj.N * obj.tau.Expectation) * eye(obj.D);

            obj.mu.updateParameters(newMu, newCov);
        end

        % obj.tau is GammaDistribution
        function obj = qTauUpdate(obj)
            deltaA = obj.N * obj.D/2;
            deltaB = 0;

            for n = 1:obj.N
                deltaB = deltaB + obj.X.getObservationNormSq(n) ...
                    + obj.mu.ExpectationXtX ...
                    + trace(obj.W.ExpectationCtC * obj.Z.ExpectationXXt{n}) ...
                    + 2 * obj.mu.ExpectationXt * obj.W.ExpectationC * obj.Z.Expectation{n} ...
                    - 2 * obj.X.getObservation(n, true) * obj.W.ExpectationC * obj.Z.Expectation{n} ...
                    - 2 * obj.X.getObservation(n, true) * obj.mu.Expectation;
            end

            deltaB = 1/2 * deltaB;

            obj.tau.updateParameters(obj.tau.prior.a + deltaA, obj.tau.prior.b + deltaB); 
        end
        
        

        %% fit() and ELBO
        function [elboVals, it, resArr] = fit(obj)
            elboVals = -Inf(1, obj.maxIter);
            resArr = cell(1, obj.maxIter);
        
            for it = 1:obj.maxIter
                obj.qWUpdate(it);
                obj.qAlphaUpdate();
                obj.qZUpdate();
                obj.qTauUpdate();
                obj.qMuUpdate();
                

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
            res.pZ = obj.Z.ExpectationLnPC;
            res.pW = obj.getExpectationLnW();
            res.pAlpha = obj.alpha.ExpectationLnPC;
            res.pMu = obj.mu.ExpectationLnP;
            res.pTau = obj.tau.ExpectationLnP;
            
            res.qZ = obj.Z.HC;
            res.qW = obj.W.HC;
            res.qAlpha = obj.alpha.HC;
            res.qMu = obj.mu.H;
            res.qTau = obj.tau.H;
            % DEBUG

            elbo = obj.getExpectationLnPX() + obj.Z.ExpectationLnPC + obj.getExpectationLnW() + ... % p(.)
                obj.alpha.ExpectationLnPC + obj.mu.ExpectationLnP + obj.tau.ExpectationLnP + ... % p(.)
                obj.Z.HC + obj.W.HC + obj.alpha.HC + obj.mu.H + obj.tau.H; % q(.)

            % DEBUG
            res.elbo = elbo;
            % DEBUG
        end

        function value = getExpectationLnPX(obj)
            value = obj.N * obj.D/2 * (obj.tau.ExpectationLn - log(2 * pi)) - obj.tau.Expectation/2 * ( ...
                obj.X.TrXtX - 2 * trace(obj.W.ExpectationC * obj.Z.ExpectationC * obj.X.data') ...
                - 2 * obj.mu.ExpectationXt * obj.X.data * ones(obj.N, 1) ...
                + 2 * obj.mu.ExpectationXt * obj.W.ExpectationC * obj.Z.ExpectationC * ones(obj.N, 1) ...
                + trace(obj.W.ExpectationCtC * obj.Z.ExpectationCCt) + obj.N * obj.mu.ExpectationXtX);

        end

        function value = getExpectationLnW(obj)
            value = 0;
            colsNormSq = obj.W.getExpectationOfColumnsNormSq();
            for k = 1:obj.K % TODO (high): This can be implemented as a dot product
                value = value + obj.alpha.Expectation{k} * colsNormSq(k);
            end
            value = -1/2 * value + obj.D/2 * (obj.alpha.ExpectationLnC - obj.K * log(2*pi));
        end

        %% Getters
        function value = get.N(obj)
            value = obj.X.N;
        end
        
        function value = get.D(obj)
            value = obj.X.D;
        end
    end
end
