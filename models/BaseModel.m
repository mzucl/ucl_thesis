classdef (Abstract) BaseModel < handle
    properties
        % Model parameters
        K               % Number of latent dimensions/principal components
        Z               % [K x N] GaussianContainer [size: N; for each latent variable zn]
        views           % An array of views

        % Optimization parameters
        maxIter
        tol

        doRotation
    end

    % The next two sections of variables support `dependent properties` that 
    % are initialized in the constructor and should remain immutable afterward
    properties (Dependent, SetAccess = private)
        M               % Number of views/groups
        N               % Number of observations
        D               % Array containing dimensions for each view
    end

    % Backing variables
    properties (Access = private) 
        M_
        N_
        D_
    end

    properties (Access = public, Dependent)
        W               % The `big W` matrix, containing the `W` matrices of all views
        alpha           % The `big alpha` matrix, containing the `alpha` vectors of all views
    end



    methods (Abstract)
        qZUpdate(obj) 
        computeELBO(obj)
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

        % Closely related to the commnet in the contructor
        % TODO (high): Implement
        % methods(Access = private)
        %     function [M, N] = validateSources(obj, idx)
        % 
        %     end
        % end
        % TODO (very high!): Implementent method below as soon as the
            % toy example with 2 views starts working
            % [obj.M, obj.N] = obj.validateSources(...)


        function obj = BaseModel(data, K, maxIter, tol, doRotation)
            CustomError.validateNumberOfParameters(nargin, 2, 5);

            obj.M_ = numel(data);
            obj.N_ = size(data{1}, 2);

            obj.D_ = zeros(1, obj.M);
            for i = 1:numel(data)
                obj.D_(i) = size(data{i}, 1);
            end

            obj.K = DoubleWrapper(K);
            obj.maxIter = maxIter;
            obj.tol = tol;
            obj.doRotation = doRotation;

            % Initialize views
            obj.views = BaseView.empty(obj.M, 0);
        end



        %% Update methods
        function obj = qWUpdate(obj, it)
            for m = 1:obj.M
                obj.views(m).qWUpdate(it);
            end
        end

        function obj = qAlphaUpdate(obj)
            for m = 1:obj.M
                obj.views(m).qAlphaUpdate();
            end
        end

        function obj = qMuUpdate(obj)
            for m = 1:obj.M
                obj.views(m).qMuUpdate();
            end
        end

        function obj = qTauUpdate(obj)
            for m = 1:obj.M
                obj.views(m).qTauUpdate();
            end
        end



        %% Model training
        function obj = removeFactors(obj, it, threshold)
            CustomError.validateNumberOfParameters(nargin, 2, 3);
            if nargin < 3
                threshold = Utility.getConfigValue('Model', 'LATENT_FACTORS_THRESHOLD');
            end

            % Calculate the average of the square of elements for each row of Z
            avgSquare = mean(obj.Z.E.^2, 2);

            removeIdx = find(avgSquare < threshold);
            if isempty(removeIdx)
                return;
            end

            if (RunConfig.getInstance().enableLogging)
                numRemovedFactors = length(removeIdx);
                factorText = Utility.ternary(numRemovedFactors == 1, 'factor', 'factors');
                
                fprintf('Removed %d %s in iteration %d\n', numRemovedFactors, factorText, it);
            end

            % Update number of factors
            obj.K.Val = obj.K.Val - length(removeIdx);

            % Remove those rows from Z, corresponding columns from W, and elements from alpha
            obj.Z.removeDimensions(removeIdx);
            for m = 1:obj.M
                obj.views(m).alpha.removeDimensions(removeIdx);
                obj.views(m).W.removeDimensions(removeIdx);
            end
        end

        % [IDEA] Make elboRecalcInterval adaptive, based on how close 
        % the model is to convergence.

        % [NOTE] Added because in BGFA there is also `qXiUpdate`
        function stepUpdate(obj, it)
            obj.qZUpdate();
            obj.qWUpdate(it);
            % obj.qZUpdate();
            obj.qMuUpdate();
            % if it > 0
            %     obj.updateRotation();
            % end  
            obj.qAlphaUpdate();
            obj.qTauUpdate();
            obj.removeFactors(it);
        end
        
        function [elboVals, it] = fit(obj, elboRecalcInterval)
            CustomError.validateNumberOfParameters(nargin, 1, 2);
            if nargin < 2
                elboRecalcInterval = RunConfig.getInstance().elboRecalcInterval;
            end

            elboVals = -Inf(1, obj.maxIter);

            % [NOTE] When `elboIterStep` is not equal to 1, the `elbo` array is indexed 
            % differently. Instead of using `iter`, we use `iter / elboIterStep + 1`. 
            % However, using a separate counter for clarity is preferred. The `+ 1` accounts 
            % for the fact that the ELBO is computed in the first iteration.
            elboIdx = 1;
            for it = 1:obj.maxIter
                obj.stepUpdate(it);

                if it ~= 1 && mod(it, elboRecalcInterval) ~= 0
                    continue;
                end

                currElbo = obj.computeELBO();
                elboVals(elboIdx) = currElbo;

                if elboIdx ~= 1
                    prevElbo = elboVals(elboIdx - 1);
                end
                if RunConfig.getInstance().enableLogging
                    if elboIdx ~= 1
                        fprintf('------ ELBO increased by: %.4f\n', currElbo - prevElbo);
                    end
                end

                % The ELBO must increase with each iteration. This is a critical error,
                % so it is logged regardless of the `RunConfig.getInstance().enableLogging` setting.
                if elboIdx ~= 1 && currElbo < elboVals(elboIdx - 1)
                    fprintf(2, '[ERROR] ELBO decreased in iteration %d by %f!\n', it, abs(currElbo - prevElbo));
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



        %% Getters
        function value = get.M(obj)
            value = obj.M_;
        end

        function value = get.N(obj)
            value = obj.N_;
        end

        function value = get.D(obj)
            value = obj.D_;
        end

        function value = get.W(obj)
            totalD = sum(obj.D);
            value = zeros(totalD, obj.K.Val);

            d = 0;
            for m = 1:obj.M
                Dm = obj.views(m).D;
                value(d + 1 : d + Dm, :) = obj.views(m).W.E;
                d = d + Dm;
            end
        end

        function value = get.alpha(obj)
            value = zeros(obj.K.Val, obj.M);
            for m = 1:obj.M
                value(:, m) = obj.views(m).alpha.E;
            end
        end



        % Variables ending with '_' denote expectations.
        % `X_train` and `y_train` are used solely for threshold selection.
        function [K_eff, predictions_test] = makePredictions(obj, X_train, y_train, X_test)
            Z_ = obj.Z.E;
            K_eff = size(Z_, 1);

            W1_ = obj.views(1).W.E;
            W2_ = obj.views(2).W.E;
            mu1_ = obj.views(1).mu.E;
            mu2_ = obj.views(2).mu.E;
            T1_ = obj.views(1).tau.E * eye(obj.D(1));

            sigma_Z = Utility.matrixInverse(eye(K_eff) + W1_' * T1_ * W1_);

            % Find the best threshold on the train data
            MU_Z = sigma_Z * (W1_' * T1_ * (X_train' - mu1_));

            % [NOTE] Although this is a sigmoid function, the output is not 
            % a probability; it's used purely to clip values to the [0, 1] range.
            predictions_tr = Bound.sigma(W2_ * MU_Z + mu2_);
            [fpr, tpr, thresholds, ~] = perfcurve(y_train', predictions_tr, 1);

            % Calculate G-means
            gMeans = sqrt(tpr .* (1 - fpr));
            [~, idx] = max(gMeans);
            train_best_threshold = thresholds(idx);

            % Predictions on the test data
            MU_Z = sigma_Z * (W1_' * T1_ * (X_test' - mu1_));
            predictions_test = Bound.sigma(W2_ * MU_Z + mu2_); 
            predictions_test = predictions_test >= train_best_threshold;
            predictions_test = double(predictions_test');
        end
    end   
end