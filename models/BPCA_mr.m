classdef BPCA_mr < handle
    properties 
        view            % ViewHandler
        K               % Number of latent dimensions/principal components

        % Model parameters
        Z               % [K x N] GaussianContainer      [size: N; for each latent variable zn]
        
        mu              % [D x 1] Gaussian               [D x 1; all observations have the same 'mu' parameter]
        
        W               % [D x K] GaussianContainer      [size: D; for each row in W matrix]
                        % Prior over W is defined per columns (each column
                        % has its own precision parameter, but update
                        % equations are defined by rows, so we are
                        % representing W as a size D container in a row
                        % format.

        alpha           % [K x 1]   GammaContainer         [size: K]
        tau             % [scalar]  Gamma                  [scalar]

        % Optimization parameters
        maxIter
        tol

        %% Constant dependent properties
        %   They are not declared as dependent because they never change upon the initialization
        N   % Number of observations/latent variables
        D   % Dimensionality
    end

    methods
        function obj = BPCA_mr(X, W_init, maxIter, tol)
            CustomError.validateNumberOfParameters(nargin, 1, 4);

            obj.view = ViewHandler(X, false);
            obj.N = obj.view.N;
            obj.D = obj.view.D;

            % BPCA can infer the right number of components, thus K is not
            % passed in as a parameter
            obj.K = obj.D - 1;


            % Set default values
            if nargin < 4, tol = Utility.getConfigValue('Optimization', 'DEFAULT_TOL'); end
            if nargin < 3, maxIter = Utility.getConfigValue('Optimization', 'DEFAULT_MAX_ITER'); end
            if nargin < 2, W_init = randn(obj.K, obj.D)'; end

            % TODO: Why transpose above: mean of W is stored in columns even thought cols = false for obj.W
            
            obj.maxIter = maxIter;
            obj.tol = tol;
           

            %% Model setup and initialization
            %                         type, size_, cols, dim,     mu, cov, priorPrec
            % obj.Z = GaussianContainer("DS", obj.N, true, obj.K, zeros(obj.K, 1));

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


        
        %% Update methods
        function obj = qConstantUpdates(obj)
            % alpha.a
            obj.alpha.updateAllDistributionsA(obj.alpha.prior.a + obj.D/2);
            % tau.a
            obj.tau.updateA(obj.tau.prior.a + (obj.N * obj.D)/2);
        end


        % function obj = qZUpdate(obj, it)
        %     tauExp = Utility.ternary(it == 1, obj.tau.getExpInit(), obj.tau.E);
        %     muExp = Utility.ternary(it == 1, obj.mu.getExpInit(), obj.mu.E);
        % 
        %     covNew = Utility.matrixInverse(eye(obj.K) + tauExp * obj.W.E_XtX);
        %     muNew = tauExp * covNew * obj.W.E_Xt * (obj.view.X - muExp);
        % 
        %     obj.Z.updateDistributionsParameters(muNew, covNew);
        % end

        
        function obj = qWUpdate(obj, it)
            alphaExp = Utility.ternary(it == 1, obj.alpha.getExpInit(), obj.alpha.E);
            tauExp = Utility.ternary(it == 1, obj.tau.getExpInit(), obj.tau.E);
            muExp = Utility.ternary(it == 1, obj.mu.getExpInit(), obj.mu.E);

            covNew = Utility.matrixInverse(diag(alphaExp) + tauExp * obj.Z.E_XXt);
            muNew = tauExp * covNew * obj.Z.E * (obj.view.X' - muExp');
            
            obj.W.updateDistributionsParameters(muNew, covNew);
        end


        function obj = qAlphaUpdate(obj)
            obj.alpha.updateAllDistributionsB(obj.alpha.prior.b + 1/2 * obj.W.E_SNC);
        end


        function obj = qMuUpdate(obj)
            covNew = (1/(obj.mu.priorPrec + obj.N * obj.tau.E)) * eye(obj.D);
            muNew = obj.tau.E * covNew * (obj.view.X - obj.W.E * obj.Z.E) * ones(obj.N, 1);
            
            obj.mu.updateParameters(muNew, covNew);
        end


        function obj = qTauUpdate(obj)
            expWtW_tr = obj.W.E_XtX';
            expZZt = obj.Z.E_XXt;
            expWZ = obj.W.E * obj.Z.E;
            
            bNew = obj.tau.prior.b + 1/2 * obj.view.Tr_XtX + obj.N/2 * obj.mu.E_XtX + ...
                1/2 * sum(expWtW_tr(:) .* expZZt(:)) - sum(obj.view.X(:) .* expWZ(:)) + ...
                obj.mu.E' * (expWZ - obj.view.X) * ones(obj.N, 1);
        
            obj.tau.updateB(bNew);
        end
        
        

        %% fit() and ELBO
        % elboIterStep - specifies the interval at which the ELBO should 
        % be computed; e.g. if elboIterStep = 2 elbo will be computed every
        % second iteration.
        function [elboVals, it] = fit(obj, elboIterStep)
            if nargin < 2
                elboIterStep = 1;
            end

            elboVals = -Inf(1, obj.maxIter);
            % When elboIterStep ~= 1, indexing into elbo array is not done
            % using 'iter'; iter / elboIterStep + 1, but having independent
            % counter is cleaner; '+ 1' because we compute elbo in the
            % first iteration.
            elboIdx = 1;
            
            for it = 1:obj.maxIter
                obj.qZUpdate(it);
                obj.qWUpdate(it);
                obj.qMuUpdate();
                obj.qAlphaUpdate();
                obj.qTauUpdate();

                if it ~= 1 && mod(it, elboIterStep) ~= 0
                    continue;
                end

                currElbo = obj.computeELBO();
                elboVals(elboIdx) = currElbo;
                
                if RunConfig.getInstance().enableLogging
                    if elboIdx ~= 1
                        disp(['======= ELBO increased by: ', num2str(currElbo - elboVals(elboIdx - 1))]);
                    end
                end

                % ELBO has to increase from iteration to iteration
                if elboIdx ~= 1 && currElbo < elboVals(elboIdx - 1)
                    fprintf(2, 'ELBO decreased in iteration %d by %f\n!!!', it, abs(currElbo - elboVals(elboIdx - 1)));
                end 

                % Check for convergence
                if elboIdx ~= 1 && abs(currElbo - elboVals(elboIdx - 1)) / abs(currElbo) < obj.tol
                    disp(['Convergence at iteration: ', num2str(it)]);
                    elboVals = elboVals(1:elboIdx); % cut the -Inf values at the end
                    break;
                end
                elboIdx = elboIdx + 1;
                if it == obj.maxIter
                    fprintf(2, 'Model did not converge in %d\n!!!', obj.maxIter);
                end
            end
        end
        

        function elbo = computeELBO(obj)
            elbo = obj.getExpectationLnPX() + obj.Z.E_LnP + obj.getExpectationLnPW() + ... % p(.)
                obj.alpha.E_LnP + obj.mu.E_LnP + obj.tau.E_LnP + ... % p(.)
                obj.Z.H + obj.W.H + obj.alpha.H + obj.mu.H + obj.tau.H; % q(.)
        end


        function value = getExpectationLnPX(obj)
            % Setup
            expWtW_tr = obj.W.E_XtX';
            expZZt = obj.Z.E_XXt;
            expW_tr = obj.W.E';
            expZXt = obj.Z.E * obj.view.X';

            value = obj.N * obj.D/2 * (obj.tau.E_LnX - log(2 * pi)) - obj.tau.E/2 * ( ...
                obj.view.Tr_XtX - 2 * sum(expW_tr(:) .* expZXt(:)) ...
                - 2 * obj.mu.E_Xt * (obj.view.X * ones(obj.N, 1)) ...
                + 2 * (expW_tr * obj.mu.E)' * (obj.Z.E * ones(obj.N, 1)) ...
                + sum(expWtW_tr(:) .* expZZt(:)) + obj.N * obj.mu.E_XtX);
        end


        function value = getExpectationLnPW(obj)
            value = -1/2 * dot(obj.W.E_SNC, obj.alpha.E) + ...
                + obj.D/2 * (obj.alpha.E_LnX - obj.K * log(2*pi));
        end








        % REDUCER (sum all emitted stats)
        function modelReducer(~, statsList, outKV)
            sumStats = statsList{1};
            for i = 2:length(statsList)
                sumStats = sumStats + statsList{i};
            end
            add(outKV, 'aggregatedStats', sumStats);
        end

        % function stats = computeSufficientStats(Xc, Zc, W)
        %     % Xc: [D x Nc] data chunk
        %     % Zc: GaussianContainer object, size: Nc, dim: Kc
        %     % W:  [D x K] global factor loading matrix        
        % 
        % end

        % % Each worker handles its own `Zc` instance
        % function modelMapper(data, ~, intermKVStore, params)
        %     persistent Zc;
        %     persistent Xc;
        % 
        %     W_ = params.W;
        %     mu_ = params.mu;
        %     tau_ = params.tau;
        %     it_ = params.it;
        % 
        %     Kc = size(W_, 2); % Same as obj.K.Val
        %     Nc = size(data, 1);
        % 
        %     if isempty(Zc)
        %         initZMu = randn(Kc, 1);
        %         %                      type, size_, cols, dim,   mu
        %         Zc = GaussianContainer("DS",  Nc,   true,  Kc, initZMu);
        %     end
        % 
        %     if isempty(Xc)
        %         Xc = ViewHandler(data, false);
        %     end
        % 
        % 
        %     % STEP 1: qZUpdate
        %     % TODO: Covariance is shared, don't recompute it! HOW???
        %     tauExp = Utility.ternary(it_ == 1, tau_.getExpInit(), tau_.E);
        %     muExp = Utility.ternary(it_ == 1, mu_.getExpInit(), mu_.E);
        % 
        %     covNew = Utility.matrixInverse(eye(Kc) + tauExp * W_.E_XtX);
        %     muNew = tauExp * covNew * W_.E_Xt * (Xc.X - muExp);
        % 
        %     Zc.updateDistributionsParameters(muNew, covNew);
        % 
        %     % STEP 2: Compute sufficient stats for reducer
        %     stats = computeSufficientStats(Xc, Zc, params);
        %     add(intermKVStore, 'stats', stats);
        % end



        % Updated FIT method using mapreduce and local GaussianContainer for Z
        function [elboVals, it] = fit2(obj, elboIterStep)
            if nargin < 2
                elboIterStep = 1;
            end
        
            elboVals = -Inf(1, obj.maxIter);
            elboIdx = 1;
        
            for it = 1:obj.maxIter
                % Save global params to shared storage (if needed)
                params.mu = obj.mu;
                params.W = obj.W;
                params.tau = obj.tau;
                params.alpha = obj.alpha;
                params.iter = it;
        

                %% Run mapreduce
                % data = tabularTextDatastore('Xdata.csv');
                % data.ReadSize = 100;  % Read 1000 rows at a time


                data = tabularTextDatastore('dataChunks/*.csv');
                data.ReadSize = 'file';  % Now each file is treated as a chunk

                % Here we will use an anonymous function as our mapper. This function
                % definition includes the value of model params (`W`, `mu`, `tau`) computed in the previous iteration.
                mapper = @(data, info, intermKVStore) modelMapper(data, info, intermKVStore, params);
                result = mapreduce(data, mapper, @modelReducer); % 'Display', 'off'
        
        
 

          


                % Extract and aggregate results from reducer output
                tbl = readall(result);
                stats = tbl.Value{1};
                obj.aggregateStats(tbl);
        
                % Update global variables
                obj.qWUpdate(it);
                obj.qMuUpdate();
                obj.qAlphaUpdate();
                obj.qTauUpdate();
        
                % Compute ELBO if needed
                if it ~= 1 && mod(it, elboIterStep) ~= 0
                    continue;
                end
        
                currElbo = obj.computeELBO();
                elboVals(elboIdx) = currElbo;
        
                if RunConfig.getInstance().enableLogging
                    if elboIdx ~= 1
                        disp(['======= ELBO increased by: ', num2str(currElbo - elboVals(elboIdx - 1))]);
                    end
                end
        
                if elboIdx ~= 1 && currElbo < elboVals(elboIdx - 1)
                    fprintf(2, 'ELBO decreased in iteration %d by %f\n!!!', it, abs(currElbo - elboVals(elboIdx - 1)));
                end
        
                if elboIdx ~= 1 && abs(currElbo - elboVals(elboIdx - 1)) / abs(currElbo) < obj.tol
                    disp(['Convergence at iteration: ', num2str(it)]);
                    elboVals = elboVals(1:elboIdx);
                    break;
                end
                elboIdx = elboIdx + 1;
        
                if it == obj.maxIter
                    fprintf(2, 'Model did not converge in %d\n!!!', obj.maxIter);
                end
            end
        end
        

        
        



    end
end
