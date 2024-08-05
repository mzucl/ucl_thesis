classdef BayesianPCA < handle
    properties 
        X               % ViewHandler
        K               % Number of latent dimensions/principal components

        % Model parameters with prior distributions
        Z               % GaussianDistributionContainer      [size: N; for each latent variable zn]
        mu              % Gaussian                           [D x 1; all observations have the same 'mu' parameter]
        W               % GaussianDistributionContainer      [size: D; for each row in W matrix]
                        % Prior over W is defined per columns (each column
                        % has its own precision parameter, but update
                        % equations are defined by rows, so we are
                        % representing W as a size D container in a row
                        % format.

        alpha           % GammaDistributionContainer         [size: K]
        tau             % Gamma                              [scalar]

        % Optimization parameters: 
        %   maximum number of iterations; 
        %   convergence tolerance;
        %   number of iterations until convergance;
        maxIter
        tol

        % Hyperparameters
        betaParam       % scalar
    end

    properties (Dependent)
        N
        D
    end
    
    methods
        function obj = BayesianPCA(X, K, maxIter, tol)
            % TODO: K can be optional as well, set to the D - 1, so we
            % should use obj.K instead of K below; Also, implement this!

            % Optional parameters: maxIter, tol
            switch nargin
                case {0, 1}
                    error(['Error in class ' class(obj) ': Too few arguments passed.']);
                case 2
                    % Set other parameters to default values
                    obj.maxIter = Constants.DEFAULT_MAX_ITER;
                    obj.tol = Constants.DEFAULT_TOL;
                case 3
                    obj.maxIter = maxIter;
                    obj.tol = Constants.DEFAULT_TOL;
                case 4
                    obj.maxIter = maxIter;
                    obj.tol = tol;
            end

            obj.X = ViewHandler(X);
            obj.K = K;
            
            % Init hyperparameters
            obj.betaParam = Constants.DEFAULT_GAUSS_PRECISION;

                
            % Z
            % dim, cols, numOfDistributions
            obj.Z = GaussianDistributionContainer(K, true, obj.N);

            % mu
            % Constructor parameters: mu, cov, priorPrec, dim
            obj.mu = GaussianDistribution(0, 1/Constants.DEFAULT_GAUSS_PRECISION, Constants.DEFAULT_GAUSS_PRECISION, obj.D);

            % alpha
            alphaPrior = GammaDistribution(Constants.DEFAULT_GAMMA_A, Constants.DEFAULT_GAMMA_B);
            obj.alpha = GammaDistributionContainer(Constants.DEFAULT_GAMMA_A, Constants.DEFAULT_GAMMA_B, alphaPrior, K);
            
            % W: 
            % #distributions: D 
            % dim: K
            % Each distribution describes a row in a matrix W (cols = false);
            obj.W = GaussianDistributionContainer(obj.K, false, obj.D, diag(obj.alpha.Value));
            
            % tau
            tauPrior = GammaDistribution(Constants.DEFAULT_GAMMA_A, Constants.DEFAULT_GAMMA_B);
            obj.tau = GammaDistribution(tauPrior);
        end
        
        function obj = fit(obj)
            elboVals = -Inf(1, obj.maxIter);
        
            for it = 1:obj.maxIter
                % obj.qWUpdate();
                % obj.qZUpdate();
                obj.qTauUpdate();
                obj.qAlphaUpdate();
                % obj.qMuUpdate();
                

                
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
            disp(obj.computeELBO());
            
            % % Display the principal components
            % disp('Principal Components:');
            % disp(obj.W);
        end
        
        function obj = qZUpdate(obj)
            % Update variational parameters for q(z)
            for n = 1:obj.N
                % All latent variables have the same covariance
                muNew = obj.tau.Expectation * obj.Z.distributions(n).cov * obj.W.ExpectationCt * (obj.X.getObservation(n) - obj.mu.Expectation);
                obj.Z.updateDistributionMu(n, muNew);
            end

            % Covariance update
            covNew = Utility.matrixInverse(eye(obj.K) + obj.tau.Expectation * ... 
                obj.W.ExpectationCtC);
            obj.Z.updateAllDistributionsCovariance(covNew);
            
        end
        
        function obj = qWUpdate(obj)
            % Update variational parameters for q(W)

            %% [NOTE] Initial distribution is defined per columns of W, but
            % update equations are defined per rows!

            % Covariance update
            covNew = zeros(obj.K);
            for n = 1:obj.N
                covNew = covNew + obj.Z.ExpectationXXt{n};
            end

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

        function obj = qAlphaUpdate(obj)
            % Update variational parameters for q(alpha)
            newAVals = arrayfun(@(d) d.prior.a, obj.alpha.distributions) + obj.D/2;
            newBVals = arrayfun(@(d) d.prior.b, obj.alpha.distributions) + 1/2 * obj.W.getExpectationOfColumnNormSq();
            
            obj.alpha.updateAllDistributionsParams(newAVals, newBVals, false);
        end

        function obj = qMuUpdate(obj)
            % Update variational parameters for q(mu)
            newMu = zeros(obj.D, 1);
           
            for n = 1:obj.N
                newMu = newMu + obj.X.getObservation(n) - obj.W.ExpectationC * obj.Z.Expectation{n};
            end
            newMu = obj.tau.Expectation * obj.mu.cov * newMu;
            newCov = 1./(obj.betaParam + obj.N * obj.tau.Expectation) * eye(obj.D);

            obj.mu.updateParameters(newMu, newCov);
            
        end

        function obj = qTauUpdate(obj)
            % Update variational parameters for q(tau)
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

            obj.tau.updateParameters(obj.tau.prior.a + deltaA, obj.tau.prior.b + deltaB, false); 
        end
        
        
        %% ELBO
        

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
