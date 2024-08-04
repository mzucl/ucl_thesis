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
        convIter    
        elboVals

        % Hyperparameters
        alphaParams     % struct('a', 'b')
        betaParam       % scalar
    end

    properties (Dependent)
        N
        D
    end
    
    methods
        function obj = BayesianPCA(X, K)
            obj.X = ViewHandler(X);
            obj.K = K;
            
            % Init hyperparameters
            obj.alphaParams = struct('a', Constants.DEFAULT_GAMMA_A, 'b', Constants.DEFAULT_GAMMA_B);
            obj.betaParam = Constants.DEFAULT_GAUSS_PRECISION;

            % Init optimization parameters
            obj.maxIter = 15;
            obj.tol = 1e-6;
            obj.convIter = 0;
            obj.elboVals = -Inf(1, obj.maxIter);
                
            % Z
            % dim, cols, numOfDistributions
            obj.Z = GaussianDistributionContainer(K, true, obj.N);

            % mu
            obj.mu = GaussianDistribution(0, 1./obj.betaParam, obj.D);

            % alpha
            obj.alpha = GammaDistributionContainer(obj.alphaParams.a, obj.alphaParams.b, K);
            
            % W: 
            % #distributions: D 
            % dim: K
            % Each distribution describes a row in a matrix W (cols = false);
            obj.W = GaussianDistributionContainer(obj.K, false, obj.D, diag(obj.alpha.Value));
            
            % tau
            tauPrior = GammaDistribution(Constants.DEFAULT_GAMMA_A, Constants.DEFAULT_GAMMA_B);
            obj.tau = GammaDistribution(Constants.DEFAULT_GAMMA_A, Constants.DEFAULT_GAMMA_B, tauPrior);
        end
        
        function obj = fit(obj)
            for it = 1:obj.maxIter
                obj.convIter = obj.convIter + 1;

                obj.qWUpdate();
                obj.qZUpdate();
                obj.qTauUpdate();
                obj.qAlphaUpdate();
                obj.qMuUpdate();
                

                
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
            obj.alpha.updateAllDistributionsParams(obj.D / 2, 1/2 * obj.W.getExpectationOfColumnNormSq());
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

            obj.tau.updateParameters(deltaA, deltaB); 
        end
        
        function ELBO = computeELBO(obj)
            % Compute the Evidence Lower Bound
            % This is a simplified placeholder for actual calculations
            ELBO = -sum(sum((obj.X - obj.mu * obj.W').^2));
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
