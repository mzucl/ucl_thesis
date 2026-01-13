classdef BaseModel < handle
    properties
        K               % Number of latent dimensions/principal components
        Z               % [K x N] GaussianContainer [size: N; for each latent variable zn]
        views           % An array of views
        bound           % Used for binary views

        % Optimization parameters
        maxIter
        tol

        doRotation
    end


    
    % The next two sections of variables support `dependent properties` that 
    % are initialized in the constructor and should remain immutable afterward
    properties (Dependent, SetAccess = private)
        M               % Number of views/groups
        Mc              % Number of continuous views/groups
        Mb              % Number of binary views/groups
        N               % Number of observations
        D               % Array containing dimensions for each view
    end

    % Backing variables
    properties (Access = private) 
        M_
        Mc_
        Mb_
        N_
        D_
    end

    properties (Access = public, Dependent)
        W               % The `big W` matrix, containing `W` matrix of all views
        alpha           % The `big alpha` matrix, conrainint `alpha` vector of all vierws
    end

    







    % TODO (high): Implement
    % methods(Access = private)
    %     function [M, N] = validateSources(obj, idx)
    % 
    %     end
    % end













    properties(Access = private, Constant)
        SETTINGS = ModelSettings.getInstance();
    end

    methods
        % [NOTE] We need to deal with the form of the datasets here (e.g.
        % are they in table format, or we need to import
        % data from .csv files. Also, we should check if all datasets have the
        % same number of observations.
        %
        % [path, featureInCols: true/false, data] triplet can be used to
        % describe sources, where if path is None data is in 'data' and
        % viceversa

        % [NOTE] For now we have 2 datasets passed in as matrices
        %% Constructors
        function obj = BaseModel(data, Mc, K, bound, maxIter, tol, doRotation) 
            % TODO: Validate Mc somehow?
            % bound = 'J' or 'B'
            % Mc -> number of continous views

            % Optional parameters: maxIter, tol, doRotation
            if nargin < 4
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
            end

            % TODO (very high!): Implementent method below as soon as the
            % toy example with 2 views starts working
            % [obj.M, obj.N] = obj.validateSources(...)
            obj.Mc = Mc;
            obj.M = length(data);
            obj.N = size(data{1}, 2); % Data passed in is DxN

            obj.K = DoubleWrapper(K);

            obj.bound = bound;

            % Default values of optional parameters
            obj.maxIter = BaseModel.SETTINGS.DEFAULT_MAX_ITER;
            obj.tol = BaseModel.SETTINGS.DEFAULT_TOL;
            obj.doRotation = false;

            if nargin > 4
                obj.maxIter = maxIter;
                if nargin > 5
                    obj.tol = tol;
                    if nargin > 6
                        obj.doRotation = doRotation;
                    end
                end
            end

            %% Model setup and initialization
            if bound == 'B'
                %                         type, size_, cols, dim,     mu, cov, priorPrec
                obj.Z = GaussianContainer("DS", obj.N, true, obj.K.Val, zeros(obj.K.Val, 1)); % STEP1
            elseif bound == 'J'
                %                         type, size_, cols,   dim,              mu,            -- cov, priorPrec
                obj.Z = GaussianContainer("DD", obj.N, true, obj.K.Val, randn(obj.K.Val, obj.N)); % STEP1
            end

            % obj.views(1:obj.M) = View(); % Preallocate with dummy `View` objects

            % Preallocate - DONT for now!
            % obj.views = BaseView.empty(obj.M, 0);
            obj.views = {};

            for i = 1:obj.Mc
                obj.views{i} = SGFAGroup(data{i}, obj.Z, obj.K, false); % featuresInCols = false;
            end

            for i = obj.Mc + 1:obj.M
                obj.views{i} = BinaryView(data{i}, obj.Z, obj.K, false, obj.bound); % featuresInCols = false;
            end

            obj.D = [cellfun(@(x) x.D, obj.views)];
            % obj.D = [obj.views.D];
        end



        %% Update methods
        % function obj = qZUpdate(obj)
        %     if obj.bound == 'B'
        %         covNew = zeros(obj.K.Val);
        %         muNew = zeros(obj.K.Val, obj.N);
        % 
        %         for m = 1:obj.Mc
        %             view = obj.views{m};
        %             covNew = covNew + view.tau.E * view.W.E_XtX;
        %             muNew = muNew + view.tau.E * view.W.E_Xt * (view.X.X - view.mu.E);
        %         end
        % 
        %         for m = obj.Mc + 1:obj.M
        %             view = obj.views{m};
        %             covNew = covNew + 1/4 * view.W.E_XtX;
        %             muNew = muNew + view.W.E_Xt * (view.X.X + view.bound.T - 1/4 * view.mu.E);
        %         end
        % 
        %         covNew = LinearAlgebra.inverseLU(eye(obj.K.Val) + covNew);
        %         muNew = covNew * muNew;
        % 
        %         obj.Z.updateDistributionsParameters(muNew, covNew);
        %     elseif obj.bound == 'J'
        %         % TODO
        %     end
        % end

        function obj = qWUpdate(obj, it)
            for i = 1:obj.M
                obj.views{i}.qWUpdate(it);
            end
        end

        function obj = qAlphaUpdate(obj)
            for i = 1:obj.M
                obj.views{i}.qAlphaUpdate();
            end
        end

        function obj = qMuUpdate(obj)
            for i = 1:obj.M
                obj.views{i}.qMuUpdate();
            end
        end

        function obj = qTauUpdate(obj)
            for i = 1:obj.Mc % for continuous views only
                obj.views{i}.qTauUpdate();
            end
        end

        function obj = qXiUpdate(obj)
            for i = obj.Mc + 1:obj.M % for binary views only
                obj.views{i}.qXiUpdate();
            end
        end



        %% fit() and ELBO
        function [elboVals, it] = fit(obj, elboIterStep)
            if nargin < 2
                elboIterStep = 1;
            end

            elboVals = -Inf(1, obj.maxIter);
            % [NOTE] When elboIterStep ~= 1, indexing into elbo array is not done
            % using 'iter'; iter / elboIterStep + 1, but having independent
            % counter is cleaner; '+ 1' because we compute elbo in the
            % first iteration.
            elboIdx = 1;

            for it = 1:obj.maxIter
                obj.qZUpdate();
                obj.qWUpdate(it);
                obj.qMuUpdate();
                obj.qXiUpdate();
                % obj.qZUpdate();
                % if it > 0
                %     obj.updateRotation();
                % end  
                obj.qAlphaUpdate();
                obj.qTauUpdate();

                obj.removeFactors(it);

                if it ~= 1 && mod(it, elboIterStep) ~= 0
                    continue;
                end

                currElbo = obj.computeELBO();
                elboVals(elboIdx) = currElbo;

                if BaseModel.SETTINGS.DEBUG
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
            elbo = 0;
            for m = 1:obj.M
                % p
                view = obj.views{m};
                elbo = elbo + view.getExpectationLnPX() + view.getExpectationLnW() ... % p(.)
                    + view.alpha.E_LnP + view.mu.E_LnP + ... % p(.)
                    + view.W.H + view.alpha.H + view.mu.H; % q(.)
            end

            % 'tau' for continuous views only
            for m = 1:obj.Mc
                view = obj.views{m};
                elbo = elbo + view.tau.E_LnP + view.tau.H;
            end

            elbo = elbo + obj.Z.H + obj.Z.E_LnP;
        end





        %% Additional methods
        function obj = removeFactors(obj, it, threshold)
            if nargin < 3
                threshold = Utility.getConfigValue('Model', 'LATENT_FACTORS_THRESHOLD');
            end
            % Calculate the average of the square of elements for each row of Z
            avgSquare = mean(obj.Z.E.^2, 2);
        
            removeIdx = find(avgSquare < threshold);

            if isempty(removeIdx)
                return;
            end

            % if obj.DEBUG
            disp(['Removed ', num2str(length(removeIdx)), ' factors in iteration ', num2str(it)]);
            % end
            
            % Update number of factors
            obj.K.Val = obj.K.Val - length(removeIdx);
        
            % Remove those rows from Z, corresponding columns from W, and elements from alpha
            obj.Z.removeDimensions(removeIdx);
            for m = 1:obj.M
                obj.views{m}.alpha.removeDimensions(removeIdx);
                obj.views{m}.W.removeDimensions(removeIdx);
            end
        end
    




        %% Getters
        function value = get.M(obj)
            value = obj.M_;
        end

        function value = get.Mc(obj)
            value = obj.Mc_;
        end

        function value = get.Mb(obj)
            value = obj.Mb_;
        end

        function value = get.N(obj)
            value = obj.N_;
        end

        function value = get.D(obj)
            value = obj.D_;
        end











        % Variables with '_' are expectations
        % X_tr and y_tr are used to set the threshold
        function [K_eff, predictions_te] = makePredictions(obj, X_tr, y_tr, X_te)
            Z_ = obj.Z.E;
            K_eff = size(Z_, 1);
            
            W1_ = obj.views{1}.W.E;
            W2_ = obj.views{2}.W.E;
            mu1_ = obj.views{1}.mu.E;
            mu2_ = obj.views{2}.mu.E;
            T1_ = obj.views{1}.tau.E * eye(obj.D(1));
            
            sigma_Z = LinearAlgebra.inverseLU(eye(K_eff) + W1_' * T1_ * W1_);
    
            % Find the best threshold on the train data
            MU_Z = sigma_Z * (W1_' * T1_ * (X_tr' - mu1_));
            
            predictions_tr = Bound.sigma(W2_ * MU_Z + mu2_);
            [fpr, tpr, thresholds, ~] = perfcurve(y_tr', predictions_tr, 1);
    
            % Calculate G-means
            gMeans = sqrt(tpr .* (1 - fpr));
            [~, idx] = max(gMeans);
            train_best_threshold = thresholds(idx);
            
            % Predictions on the test data
            MU_Z = sigma_Z * (W1_' * T1_ * (X_te' - mu1_));
            predictions_te = Bound.sigma(W2_ * MU_Z + mu2_);
            predictions_te = predictions_te >= train_best_threshold;
            predictions_te = double(predictions_te');
        end
    end
end