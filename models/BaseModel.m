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

    properties (Access = private)
        % [NOTE] Not a dependent property because it needs to be sorted at
        % some point
        W               % The `big W` matrix, containing the `W` matrices of all views
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
        alpha           % The `big alpha` matrix, containing the `alpha` vectors of all views
        totalVar
        factorsVar
        varWithin
        nonRelVarWithin
        relVarWithin
        factors
    end



    methods (Abstract)
        qZUpdate(obj, it) 
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
        function obj = removeFactors_v1(obj, it, threshold)
            CustomError.validateNumberOfParameters(nargin, 2, 3);
            if nargin < 3
                threshold = Utility.getConfigValue('Model', 'LATENT_FACTORS_THRESHOLD');
            end

            % Calculate the average of the square of elements for each row of Z
            avgSquare = mean(obj.Z.E.^2, 2);

            if it < 100
                return;
            end

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

        function removeFactors(obj, it, threshold)
            CustomError.validateNumberOfParameters(nargin, 2, 3);
            if nargin < 3
                threshold = Utility.getConfigValue('Model', 'LATENT_FACTORS_THRESHOLD');
            end

            % if it < 50
            %     return;
            % end
            
            W_all = cell(1, obj.M);
            for m = 1:obj.M
                W_all{m} = obj.views(m).W.E;
            end

            % Get column norms for each matrix
            colNorms = cellfun(@(W) sqrt(sum(W.^2, 1)), W_all, 'UniformOutput', false);
            colNormsMat = vertcat(colNorms{:});
            colsAllSmall = all(colNormsMat < threshold, 1);
            
            numRemovedFactors = sum(colsAllSmall);
            if numRemovedFactors == 0
                return;
            end

            if (RunConfig.getInstance().enableLogging)
                factorText = Utility.ternary(numRemovedFactors == 1, 'factor', 'factors');
                fprintf('Removed %d %s in iteration %d\n', numRemovedFactors, factorText, it);
            end

            % Update number of factors
            obj.K.Val = obj.K.Val - numRemovedFactors;

            removeIdx = find(colsAllSmall);

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
            obj.qWUpdate(it);
            obj.qZUpdate(it);
            % obj.qZUpdate(5);
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
                    obj.initW();
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

        function value = getW(obj)
            value = obj.W;
        end

        function value = get.alpha(obj)
            value = zeros(obj.K.Val, obj.M);
            for m = 1:obj.M
                value(:, m) = obj.views(m).alpha.E;
            end
        end

        function val = get.totalVar(obj)
            val = 0;
            for m = 1:obj.M
                view = obj.views(m);
                W_ = view.W.E;
                tau_ = view.tau.E;
                val = val + trace(W_ * W_' + (1 / tau_) * eye(size(W_, 1)));
            end
        end
    
        function val = get.factorsVar(obj)
            val = 0;
            for m = 1:obj.M
                view = obj.views(m);
                W_ = view.W.E;
                val = val + trace(W_ * W_');
            end
        end

        function val = get.nonRelVarWithin(obj)
            val = zeros(obj.M, obj.K.Val);
            for m = 1:obj.M
                W_ = obj.views(m).W.E;
                for k = 1:obj.K.Val
                    val(m, k) = sum(W_(:, k).^2);
                end
            end
        end
    
        function val = get.varWithin(obj)
            val = obj.nonRelVarWithin / obj.totalVar * 100;
        end
    
        function val = get.relVarWithin(obj)
            val = zeros(obj.M, obj.K.Val);
            vw = obj.varWithin;
            for m = 1:obj.M
                rowSum = sum(vw(m, :));
                if rowSum > 0
                    val(m, :) = vw(m, :) / rowSum * 100;
                end
            end
        end

        function val = getFactors(obj)
            
        end

        % [NOTE] Variables ending with '_' denote expectations
        function [K_eff, predictionsTest] = makePredictions(obj, views, trainIdx, testIdx)
            outputView = views{end};
            inputViews = views(1:end-1);

            N_total = size(outputView, 2);  % Total number of observations (train + test) 
                                % [NOTE] obj.N refers only to the number of training samples
            K_eff = size(obj.Z.E, 1);
            sigma_Z = eye(K_eff);
            mu_Z = zeros(K_eff, N_total); 

            for m = 1:numel(inputViews)
                view = obj.views(m);

                W_ = view.W.E;
                mu_ = view.mu.E;
                tau_ = view.tau.E;

                sigma_Z = sigma_Z + W_' * tau_ * W_;
                mu_Z = mu_Z + W_' * tau_ * (inputViews{m} - mu_);
            end

            sigma_Z = Utility.matrixInverse(sigma_Z);
            mu_Z = sigma_Z * mu_Z;

            % [NOTE] Although this is a sigmoid function, the output is not 
            % a probability; it's used purely to clip values to the [0, 1] range.
            predictions = Bound.sigma(obj.views(end).W.E * mu_Z + obj.views(end).mu.E);

            % Find the best threshold on the train data
            [fpr, tpr, thresholds, ~] = perfcurve(outputView(trainIdx), predictions(trainIdx), 1);

            % Calculate G-means
            gMeans = sqrt(tpr .* (1 - fpr));
            [~, idx] = max(gMeans);
            bestThreshold = thresholds(idx);

            % Predictions on the test data
            predictionsTest = double(predictions(testIdx) >= bestThreshold);
        end
    end

    methods (Access = private)
        function obj = setW(obj, newW)
            if size(newW) ~= size(obj.W)
                CustomError.raiseError('InputCheck', 'Dims!');
            end
            obj.W = newW;
        end

        function obj = initW(obj)
            totalD = sum(obj.D);
            obj.W = zeros(totalD, obj.K.Val);
        
            d = 0;
            for m = 1:obj.M
                Dm = obj.views(m).D;
                obj.W(d + 1 : d + Dm, :) = obj.views(m).W.E;
                d = d + Dm;
            end
        end
    end
end