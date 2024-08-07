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
   
            % Z
            % ------------------------------------------------------ %
            % Initialize the model - set random values for the 'mu'
            zPrior = GaussianDistribution(randn(K, 1), eye(obj.K));
            obj.Z = GaussianDistributionContainer(obj.N, zPrior, true);

            % mu
            muPrior = GaussianDistribution(0, 1/Constants.DEFAULT_GAUSS_PRECISION * eye(obj.D));
            obj.mu = GaussianDistribution(muPrior);

            % alpha
            alphaPrior = GammaDistribution(Constants.DEFAULT_GAMMA_A, Constants.DEFAULT_GAMMA_B);
            obj.alpha = GammaDistributionContainer(repmat(alphaPrior, K, 1));
            
            % W; sample from obj.alpha for the prior
            % disp(obj.alpha.Value)
            wPrior = GaussianDistribution(0, diag(1./obj.alpha.Value));
            obj.W = GaussianDistributionContainer(obj.D, wPrior, false);
            
            % tau
            tauPrior = GammaDistribution(Constants.DEFAULT_GAMMA_A, Constants.DEFAULT_GAMMA_B);
            obj.tau = GammaDistribution(tauPrior);
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
        function obj = qWUpdate(obj)
            % [NOTE] Initial distribution is defined per columns of W, but
            % update equations are defined per rows!

            % % DEBUG
            % if ~all(obj.W.ExpectationC(:) == 0)
            %     disp('W is not zero!');
            % end
            % % DEBUG

            % Covariance update
            covNew = obj.Z.ExpectationCCt;

            covNew = Utility.matrixInverse((diag(obj.alpha.ExpectationC)) + obj.tau.Expectation * covNew);

            obj.W.updateAllDistributionsCovariance(covNew);

            % Mean update
            for k = 1:obj.K
                muNew = zeros(obj.K, 1);
                for n = 1:obj.N
                    muNew = muNew + obj.Z.Expectation{n} * (obj.X.getObservationEntry(n, k) - obj.mu.Expectation(k));
                end
                muNew = obj.tau.Expectation * obj.W.distributions(k).cov * muNew;
                obj.W.updateDistributionMu(k, muNew);
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
            newCov = 1./(obj.mu.PriorPrecision + obj.N * obj.tau.Expectation) * eye(obj.D);

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
        function obj = fit(obj)
            % elboVals = -Inf(1, obj.maxIter);
        
            for it = 1:obj.maxIter
                obj.qWUpdate();
                obj.qZUpdate();
                obj.qAlphaUpdate();
                obj.qMuUpdate();
                obj.qTauUpdate();
                
                
                
                % ------------------------------------------------------------
                % !!! ELBO is not implemented completely -> ignore the code
                % below
                % ------------------------------------------------------------
                % currElbo = obj.computeELBO();

                % CHECK: ELBO has to increase from iteration to iteration
                % if it ~= 1 && currElbo < obj.elboVals(it - 1)
                %     error(['ELBO decreased in iteration ' num2str(it)]);
                % end
                
                % Check for convergence
                % if it ~= 1 && abs(currElbo - obj.elboVals(it - 1)) < obj.tol
                %     obj.elboVals = obj.elboVals(1:obj.convIter); % cut the -Inf values at the end
                %     break;
                % end
            end
        end
        function elbo = computeELBO(obj)
            % Compute the Evidence Lower Bound
            elbo = 0;
            % PART 1: p(.)
            elbo = elbo + obj.tau.ExpectationLnP + obj.alpha.ExpectationLnPC;

            % PART 2: q(.)
            elbo = elbo + obj.Z.HC + obj.W.HC + obj.alpha.HC + obj.mu.H + obj.tau.H;
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
