classdef BayesianPCA
    properties
        % Observations; number of samples; data dimensionality;
        X
        N
        D
        
        % Number of latent dimensions/principal components
        K

        % Model parameters with prior distributions
        z       % GaussianDistributionContainer      [size: N; for each latent variable zn]
        mu      % Gaussian                           [D x 1]
        W       % GaussianDistributionContainer      [size: K; for each column in W matrix]
        alpha   % GammaDistributionContainer         [size: K]
        tau     % Gamma                              [scalar]

        % Optimization parameters: 
        %   maximum number of iterations; 
        %   convergence tolerance;
        %   number of iterations until convergance;
        maxIter
        tol
        iter
        elboVals

        % Hyperparameters
        alphaParams     % struct('a', 'b')
        tauParams       % struct('a', 'b')
        betaParam       % scalar
    end
    
    methods
        function obj = BayesianPCA(X, K)
            obj.X = X;
            [obj.N, obj.D] = size(X);
            obj.K = K;
            
            % Init hyperparameters
            obj.alphaParams = struct('a', Constants.DEFAULT_GAMMA_A, 'b', Constants.DEFAULT_GAMMA_B);
            obj.tauParams = struct('a', Constants.DEFAULT_GAMMA_A, 'b', Constants.DEFAULT_GAMMA_B);
            obj.betaParam = Constants.DEFAULT_GAUSS_PRECISION;
                
            % Z
            obj.z = GaussianDistributionContainer(K, obj.N);

            % mu
            obj.mu = GaussianDistribution(0, 1./obj.betaParam, obj.D);

            % alpha
            obj.alpha = GammaDistributionContainer(obj.alphaParams.a, obj.alphaParams.b, K);
            
            disp(obj.alpha.Value);
            % W
            obj.W = GaussianDistributionContainer(obj.D, K, obj.alpha.Value);
            
            % tau
            obj.tau = GammaDistribution(obj.tauParams.a, obj.tauParams.b); % tau is a scalar

            % Init optimization parameters
            obj.maxIter = 100;
            obj.tol = 1e-6;
            obj.iter = 0;
            obj.elboVals = -Inf(1, obj.maxIter);
        end
        
        function obj = fit(obj)
            for it = 1:obj.maxIter
                obj.iter = obj.iter + 1;


                % E-step: Update q(Z)
                obj = obj.update_qZ();
                
                % M-step: Update q(W), q(tau)
                obj = obj.update_qWT();
                
                currElbo = obj.computeELBO();

                % CHECK: ELBO has to increase from iteration to iteration
                if it ~= 1 && currElbo < obj.elboVals(it - 1)
                    error(['ELBO decreased in iteration ' num2str(it)]);
                end
                
                % Check for convergence
                if it ~= 1 && abs(currElbo - obj.elboVals(it - 1)) < obj.tol
                    obj.elboVals = obj.elboVals(1:obj.iter); % cut the -Inf values at the end
                    break;
                end
            end
            
            % Display the principal components
            disp('Principal Components:');
            disp(obj.W);
        end
        
        function obj = qZUpdate(obj)
            % Update variational parameters for q(z)
            for i = 1:obj.N
                muNew = obj.tau.Expectation * obj.z
            end

            % Covariance update
            covNew = Utility.matrixInverse(eye(dim.K) + obj.tau.Expectation * ... 
                obj.W.ExpectationCtC);
            obj.z.updateAllDistributionsCovariance(covNew);
            
        end
        
        function obj = qWUpdate(obj)
            % Update variational parameters for q(W)
            cov = Utility.matrixInverse((diag(obj.alpha.Expectation)) + obj.tau.Expectation); %* TODO; FInish 
            
        end

        function obj = qAlphaUpdate(obj)
            % Update variational parameters for q(alpha)
            deltaB = 1; % TODO
            obj.alpha.updateAllDistributionsParams(obj.D/2, deltaB, 1);
        end

        function obj = qMuUpdate(obj)
            % Update variational parameters for q(mu)
            cov = 1./(obj.betaParam + obj.N * obj.tau.Expectation) * eye(obj.D);
            
        end

        function obj = qTauUpdate(obj)
            % Update variational parameters for q(tau)
            obj.tau.a = obj.tau.a + obj.N * obj.D/2;
            % obj.tau.b = ...
           
            
        end
        
        function ELBO = computeELBO(obj)
            % Compute the Evidence Lower Bound
            % This is a simplified placeholder for actual calculations
            ELBO = -sum(sum((obj.X - obj.mu * obj.W').^2));
        end
    end
end
