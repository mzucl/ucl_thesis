% TODO: Remove intermediate files from mr/ or any other
classdef BPCA_mr < handle
    properties 
        % view            % ViewHandler
        ds                  % Datastore
        K               % Number of latent dimensions/principal components

        % Model parameters
        % Z               % [K x N] GaussianContainer      [size: N; for each latent variable zn]
        
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

        stats
    end

    methods
        function obj = BPCA_mr(dataPath, W_init, maxIter, tol)
            CustomError.validateNumberOfParameters(nargin, 1, 4);

            dataPath = 'dataChunks/*.csv'; % TODO: hardcoded for now

            obj.ds = tabularTextDatastore(dataPath);
            obj.ds.ReadSize = 'file';  % Now each file is treated as a chunk

            obj.D = numel(obj.ds.SelectedVariableNames);

            % obj.view = ViewHandler(X, false);
            % obj.N = obj.view.N; -> ?? this is distrubuted now
            % obj.D = obj.view.D;

            % BPCA can infer the right number of components, thus K is not
            % passed in as a parameter
            obj.K = obj.D - 1;


            % Set default values
            if nargin < 4, tol = Utility.getConfigValue('Optimization', 'DEFAULT_TOL'); end
            if nargin < 3, maxIter = Utility.getConfigValue('Optimization', 'DEFAULT_MAX_ITER'); end
            if nargin < 2, W_init = randn(obj.D, obj.K); end

            % TODO: Why transpose above: mean of W is stored in columns even thought cols = false for obj.W
            % Check how we init W from the PPCA model
            
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
            % TODO: Check if mean sixe is validated in the constructor
            
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
            % TODO: don't have obj.N! so this can't be the first update
            % equations that are run
            % obj.qConstantUpdates();
            % obj.alpha.updateAllDistributionsA(obj.alpha.prior.a + obj.D/2);
        end


        
        %% Update methods
        function obj = qConstantUpdates(obj)
            % alpha.a
            obj.alpha.updateAllDistributionsA(obj.alpha.prior.a + obj.D/2);
            % tau.a
            obj.tau.updateA(obj.tau.prior.a + (obj.N * obj.D)/2);
        end


        
        function obj = qWUpdate(obj, it)
            alphaExp = Utility.ternary(it == 1, obj.alpha.getExpInit(), obj.alpha.E);
            tauExp = Utility.ternary(it == 1, obj.tau.getExpInit(), obj.tau.E);
            % muExp = Utility.ternary(it == 1, obj.mu.getExpInit(), obj.mu.E);

            covNew = Utility.matrixInverse(diag(alphaExp) + tauExp * obj.stats.E_ZcZct);
            muNew = tauExp * covNew * obj.stats.E_Zc_times_centered_data_t;
            
            obj.W.updateDistributionsParameters(muNew, covNew);
        end

        % FINE!
        function obj = qAlphaUpdate(obj)
            obj.alpha.updateAllDistributionsB(obj.alpha.prior.b + 1/2 * obj.W.E_SNC);
        end


        % FINE
        function obj = qMuUpdate(obj, it)
            tauExp = Utility.ternary(it == 1, obj.tau.getExpInit(), obj.tau.E);

            covNew = (1/(obj.mu.priorPrec + obj.N * tauExp)) * eye(obj.D);

            muNew = tauExp * covNew * (obj.stats.Xc_col_sum - obj.stats.E_WZc_col_sum);
            
            obj.mu.updateParameters(muNew, covNew);
        end



        % TODO: qTauUpdate and elbo (getExpectationLpPX) have a lot in
        % common
        % FINE!
        function obj = qTauUpdate(obj)
            expWtW_tr = obj.W.E_XtX';
            expZZt = obj.stats.E_ZcZct;
           

            bNew = obj.tau.prior.b + 1/2 * obj.stats.Tr_XctXc + obj.N/2 * obj.mu.E_XtX + ...
                1/2 * sum(expWtW_tr(:) .* expZZt(:)) - obj.stats.Tr_E_WZ_times_Xt + ...
                obj.mu.E' * (obj.stats.E_WZc_col_sum - obj.stats.Xc_col_sum);
        
            obj.tau.updateB(bNew);
        end
        

        %% fit() and ELBO
        % elboIterStep - specifies the interval at which the ELBO should 
        % be computed; e.g. if elboIterStep = 2 elbo will be computed every
        % second iteration.
        

        function elbo = computeELBO(obj)
            elbo = obj.getExpectationLnPX() + obj.stats.E_LnP + obj.getExpectationLnPW() + ... % p(.)
                obj.alpha.E_LnP + obj.mu.E_LnP + obj.tau.E_LnP + ... % p(.)
                obj.stats.H + obj.W.H + obj.alpha.H + obj.mu.H + obj.tau.H; % q(.)
        end


        function value = getExpectationLnPX(obj)
            % Setup
            expWtW_tr = obj.W.E_XtX';
            expZZt = obj.stats.E_ZcZct;

            value = obj.N * obj.D/2 * (obj.tau.E_LnX - log(2 * pi)) - obj.tau.E * ( ...
                1/2 * obj.stats.Tr_XctXc + obj.N/2 * obj.mu.E_XtX + ...
                1/2 * sum(expWtW_tr(:) .* expZZt(:)) - obj.stats.Tr_E_WZ_times_Xt + ...
                obj.mu.E' * (obj.stats.E_WZc_col_sum - obj.stats.Xc_col_sum));
        end


        function value = getExpectationLnPW(obj)
            value = -1/2 * dot(obj.W.E_SNC, obj.alpha.E) + ...
                + obj.D/2 * (obj.alpha.E_LnX - obj.K * log(2*pi));
        end


        % Updated FIT method using mapreduce and local GaussianContainer for Z
        function [elboVals, it] = fit(obj, elboIterStep)
            if nargin < 2
                elboIterStep = 1;
            end
        
            elboVals = -Inf(1, obj.maxIter);
            elboIdx = 1;
        
            obj.maxIter = 1200;
            for it = 1:obj.maxIter
                % Save global params to shared storage (if needed)
                params.mu = obj.mu;
                params.W = obj.W;
                params.tau = obj.tau;
                params.alpha = obj.alpha;
                params.it = it;
                params.K = obj.K;
                
                % Covariance is shared, don't recompute it!
                tauExp = Utility.ternary(it == 1, obj.tau.getExpInit(), obj.tau.E);
                params.covZ = Utility.matrixInverse(eye(obj.K) + tauExp * obj.W.E_XtX);
                
                
                % This is basically qZ update! In this setup as it is now
                % it needs to be run first!


                % Clean storage
                % TODO: Add the name as a const!!!
                % if isfolder("ZcStorage")
                %     delete("ZcStorage/Zc_*.mat");
                % else
                %     mkdir("ZcStorage");
                % end

                %% Run mapreduce
                % data = tabularTextDatastore('Xdata.csv');
                % data.ReadSize = 100;  % Read 1000 rows at a time


                % data = tabularTextDatastore('dataChunks/*.csv');
                % data.ReadSize = 'file';  % Now each file is treated as a chunk

                % Here we will use an anonymous function as our mapper. This function
                % definition includes the value of model params (`W`, `mu`, `tau`) computed in the previous iteration.
                % mapr = mapreducer(0);  % Use in-memory map reducer
                mapper = @(data, info, intermKVStore) modelMapper(data, info, intermKVStore, params);
                result = mapreduce(obj.ds, mapper, @modelReducer, ...
                    'OutputFolder', 'mr', ...
                    'Display', 'off');

                % Extract and aggregate results from reducer output
                tbl = readall(result);
                obj.stats = tbl.Value{1};



                obj.N = obj.stats.N;

                % Do const updates prev done in the constructor
                if it == 1
                    obj.qConstantUpdates();
                end




        
                % Update global variables
                obj.qWUpdate(it);
                obj.qMuUpdate(it);
                obj.qAlphaUpdate();
                obj.qTauUpdate();
        

                % Compute ELBO if needed
                if it ~= 1 && mod(it, elboIterStep) ~= 0
                    continue;
                end

                currElbo = obj.computeELBO();
                elboVals(elboIdx) = currElbo;

                prevElbo = Utility.ternaryOpt(elboIdx == 1, @()nan, @()elboVals(elboIdx - 1));
                
                % The ELBO must increase with each iteration. This is a critical error,
                % so it is logged regardless of the `RunConfig.getInstance().enableLogging` setting.
                if ~isnan(prevElbo)
                    if currElbo < prevElbo
                        fprintf(2, '[ERROR] ELBO decreased in iteration %d by %f!\n', it, abs(currElbo - prevElbo));
                    else
                        if RunConfig.getInstance().enableLogging
                            fprintf('------ ELBO increased by: %.4f\n', currElbo - prevElbo);
                        end
                    end
                end

                
                % Check for convergence
                if elboIdx ~= 1 && abs(currElbo - prevElbo) / abs(currElbo) < obj.tol
                    fprintf('### Convergence at iteration: %d\n', it);
                    elboVals = elboVals(1:elboIdx); % cut the -Inf values at the end
                    break;
                end

                elboIdx = elboIdx + 1;

                if it == obj.maxIter
                    fprintf(2, '[ERROR] Model did not converge in %d iterations!\n', obj.maxIter);
                end
            end
        end
        

        
        



    end
end
