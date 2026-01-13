classdef BPCA < handle
    properties
        % View and latent space
        view            % ViewHandler                     -- wrapper for observations
        K               % Number of latent dimensions/principal components
        Z               % [K x N] GaussianContainer       -- size: N (one per latent variable z_n)
    
        % Mean parameter
        mu              % [D x 1] Gaussian                -- shared mean across all observations
    
        % Factor loadings
        W               % [D x K] GaussianContainer       -- size: D (one per row of W)
                        % Prior is column-wise (per latent dimension),
                        % but updates are done row-wise, hence row-based container
    
        % Precision parameters
        alpha           % [K x 1] GammaContainer          -- ARD prior precision for each latent dimension
        tau             % scalar  Gamma                   -- Noise precision (shared)
    
        % Optimization settings
        maxIter         % Maximum number of iterations
        tol             % Convergence tolerance
    end



    % The next two sections of variables support `Dependent properties` that 
    % are initialized in the constructor and should remain immutable
    % afterward.
    properties (Dependent, SetAccess = private)
        N               % Number of observations/latent variables
        D               % Dimensionality of the dataset
    end

    % Backing variables
    properties (Access = private) 
        N_
        D_
    end



    methods
        function obj = BPCA(X, W_init, maxIter, tol)
            CustomError.validateNumberOfParameters(nargin, 1, 4);

            % rng(42);
            obj.view = ViewHandler(X, false);
            obj.N_ = obj.view.N;
            obj.D_ = obj.view.D;

            % BPCA can infer the right number of components, thus K is not
            % passed in as a parameter
            obj.K = obj.D - 1;


            % Set default values
            if nargin < 4, tol = Utility.getConfigValue('Optimization', 'DEFAULT_TOL'); end
            if nargin < 3, maxIter = Utility.getConfigValue('Optimization', 'DEFAULT_MAX_ITER'); end
            if nargin < 2, W_init = randn(obj.D, obj.K); end

            % TODO: Why transpose above: mean of W is stored in columns even thought cols = false for obj.W
            
            obj.maxIter = maxIter;
            obj.tol = tol;
           

            %% Model setup and initialization
            %                         type, size_, cols, dim,     mu, cov, priorPrec
            obj.Z = GaussianContainer("DS", obj.N, true, obj.K, zeros(obj.K, 1));

            %                  dim, mu,    cov,  priorPrec
            obj.mu = Gaussian(obj.D, 0, eye(obj.D), 10^3);

            %                           type, size_, a, b, prior
            obj.alpha = GammaContainer("SD", obj.K);
            
            %                         type, size_, cols,   dim,   mu, cov, priorPrec
            obj.W = GaussianContainer("DS", obj.D, false, obj.K, W_init');
            
            %
            tauPrior = Gamma( ...
                Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_A'), ...
                Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_B'));

            obj.tau = Gamma(tauPrior);

            % Init model
            obj.tau.setExpInit(10);
            obj.alpha.setExpInit(repmat(1e2, obj.K, 1));
            obj.mu.setExpInit(randn(obj.D, 1));

            % Performed only once
            obj.qConstantUpdates();
        end



        %% Getters
        function value = get.N(obj)
            value = obj.N_;
        end

        function value = get.D(obj)
            value = obj.D_;
        end
        


        %% Update methods
        function obj = qConstantUpdates(obj)
            % alpha.a
            obj.alpha.updateAllDistributionsA(obj.alpha.prior.a + obj.D/2);
            % tau.a
            obj.tau.updateA(obj.tau.prior.a + (obj.N * obj.D)/2);
        end

        % [NOTE] This method computes (1/2) * E[sum(|x_n - (W * z_n + mu)|^2)], 
        % where the sum is over all observations `x_n`. This corresponds to the 
        % expected squared reconstruction error term (scaled!). It appears in both the 
        % `qTauUpdate` and `getExpectationLnPX` methods, so it is factored out 
        % into a separate method for clarity and reuse.
        function val = expectedReconstructionLoss(obj)
            expWtW_tr = obj.W.E_XtX';
            expZZt    = obj.Z.E_XXt;
            expWZ     = obj.W.E * obj.Z.E;

            val = 1/2 * obj.view.Tr_XtX + obj.N/2 * obj.mu.E_XtX + ...
                1/2 * sum(expWtW_tr(:) .* expZZt(:)) - sum(obj.view.X(:) .* expWZ(:)) + ...
                obj.mu.E' * (expWZ - obj.view.X) * ones(obj.N, 1);
        end

        function obj = qZUpdate(obj, it)
            tauExp = LogicUtils.ternary(it == 1, obj.tau.getExpInit(), obj.tau.E);
            muExp  = LogicUtils.ternary(it == 1, obj.mu.getExpInit(), obj.mu.E);

            covNew = Utility.matrixInverse(eye(obj.K) + tauExp * obj.W.E_XtX);
            muNew  = tauExp * covNew * obj.W.E_Xt * (obj.view.X - muExp);

            obj.Z.updateDistributionsParameters(muNew, covNew);
        end

        function obj = qWUpdate(obj, it)
            alphaExp = LogicUtils.ternary(it == 1, obj.alpha.getExpInit(), obj.alpha.E);
            tauExp   = LogicUtils.ternary(it == 1, obj.tau.getExpInit(), obj.tau.E);
            muExp    = LogicUtils.ternary(it == 1, obj.mu.getExpInit(), obj.mu.E);

            covNew = Utility.matrixInverse(diag(alphaExp) + tauExp * obj.Z.E_XXt);
            muNew = tauExp * covNew * obj.Z.E * (obj.view.X' - muExp');
            
            obj.W.updateDistributionsParameters(muNew, covNew);
        end

        function obj = qAlphaUpdate(obj)
            obj.alpha.updateAllDistributionsB(obj.alpha.prior.b + 1/2 * obj.W.E_SNC);
        end

        function obj = qMuUpdate(obj, it)
            tauExp = LogicUtils.ternary(it == 1, obj.tau.getExpInit(), obj.tau.E);

            covNew = (1/(obj.mu.priorPrec + obj.N * tauExp)) * eye(obj.D);
            muNew  = tauExp * covNew * (obj.view.X - obj.W.E * obj.Z.E) * ones(obj.N, 1);
            
            obj.mu.updateParameters(muNew, covNew);
        end

        function obj = qTauUpdate(obj)
            obj.tau.updateB(obj.tau.prior.b + obj.expectedReconstructionLoss());
        end



        %% fit() and ELBO
        function [elboVals, it] = fit(obj, elboRecalcInterval)
            % Fits the model using variational inference and optionally tracks ELBO values.
            %
            % Parameters:
            %   elboRecalcInterval - Specifies how often the ELBO should be computed.
            %                        For example, if elboRecalcInterval = 5, ELBO will be
            %                        computed at iterations {1, 5, 10, ...}.
            %
            % Returns:
            %   elboVals - A vector of computed ELBO values at the specified intervals.
            %   it       - The total number of iterations executed during fitting.
            CustomError.validateNumberOfParameters(nargin, 1, 2);
            if nargin < 2
                elboRecalcInterval = RunConfig.getInstance().elboRecalcInterval;
            end

            elboVals = -Inf(1, obj.maxIter);
            
            % [NOTE] When `elboRecalcInterval` ≠ 1, ELBO is not computed at every iteration.
            % In that case, `elboVals` should be indexed using `iter / elboRecalcInterval`,
            % with an additional `+1` to account for the initial ELBO computed at iteration 1.
            % A separate counter `elboIdx` is used here instead, for clarity.
            elboIdx = 1;
            for it = 1:obj.maxIter
                obj.qZUpdate(it);
                obj.qWUpdate(it);
                obj.qMuUpdate(it);
                obj.qAlphaUpdate();
                obj.qTauUpdate();

                if it ~= 1 && mod(it, elboRecalcInterval) ~= 0
                    continue;
                end

                currElbo = obj.computeELBO();
                elboVals(elboIdx) = currElbo;
                
                prevElbo = LogicUtils.ternaryOpt(elboIdx == 1, @()nan, @()elboVals(elboIdx - 1));
                
                % [NOTE] The ELBO must increase with each iteration. This is a critical error,
                % so it is logged regardless of the `RunConfig.getInstance().verbose` setting.
                if ~isnan(prevElbo)
                    if currElbo < prevElbo
                        fprintf(2, '[ERROR] ELBO decreased in iteration %d by %f!\n', it, abs(currElbo - prevElbo));
                    else
                        if RunConfig.getInstance().verbose
                            fprintf('------ ELBO increased by: %.4f\n', currElbo - prevElbo);
                        end
                    end
                end

                % Check for convergence
                if ~isnan(prevElbo) && abs(currElbo - prevElbo) / abs(currElbo) < obj.tol
                    fprintf('### Convergence at iteration: %d\n', it);
                    elboVals = elboVals(1:elboIdx); % cut the -Inf values at the end;
                    break;
                end

                elboIdx = elboIdx + 1;

                if it == obj.maxIter
                    fprintf(2, '[ERROR] Model did not converge in %d iterations!\n', obj.maxIter);
                end
            end
        end
        
        function elbo = computeELBO(obj)
            elbo = obj.getExpectationLnPX() + obj.Z.E_LnP + obj.getExpectationLnPW() + ... % p(.)
                obj.alpha.E_LnP + obj.mu.E_LnP + obj.tau.E_LnP + ... % p(.)
                obj.Z.H + obj.W.H + obj.alpha.H + obj.mu.H + obj.tau.H; % q(.)
        end

        function value = getExpectationLnPX(obj)
            value = obj.N * obj.D/2 * (obj.tau.E_LnX - log(2 * pi)) - ...
                obj.tau.E * obj.expectedReconstructionLoss();
        end

        function value = getExpectationLnPW(obj)
            value = -1/2 * dot(obj.W.E_SNC, obj.alpha.E) + ...
                + obj.D/2 * (obj.alpha.E_LnX - obj.K * log(2*pi));
        end
    end
end