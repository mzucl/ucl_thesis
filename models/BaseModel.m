classdef (Abstract) BaseModel < handle
    properties
        % Model parameters
        K               % Number of latent dimensions/principal components for BPCA
        Z               % [K x N] GaussianContainer [size: N; for each latent variable zn]
        views           % An array of views

        % Optimization parameters
        maxIter
        tol

        doRotation
    end

    % TODO: Update the comment below, public dependent are mutable!!! and
    % for them we need backing vars becasue -> they need to be sorted (changed) at
        % some point
    % The next two sections of variables support `dependent properties` that 
    % are initialized in the constructor and should remain immutable afterward
    properties (Dependent, SetAccess = private)
        M               % Number of views/groups
        N               % Number of observations
        D               % Array containing dimensions for each view
        W               % The `big W` matrix, containing the `W` matrices of all views
        tau             % The `big tau` matrix, containing the `tau` matrices (TODO I did matrices becuase of Irina's code) of all views
    end

    % Backing variables
    properties (Access = private) 
        M_
        N_
        D_
        W_   
        tau_
        totalVar_
        factorsVar_
        varWithin_
        nonRelVarWithin_
        relVarWithin_
    end

    properties (Access = public, Dependent)
        alpha           % The `big alpha` matrix, containing the `alpha` vectors of all views
        totalVar
        factorsVar
        varWithin
        nonRelVarWithin
        relVarWithin
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

                prevElbo = Utility.ternaryOpt(elboIdx == 1, @()nan, @()elboVals(elboIdx - 1));
                
                % The ELBO must increase with each iteration. This is a critical error,
                % so it is logged regardless of the `RunConfig.getInstance().enableLogging` setting.
                if ~isnan(prevElbo) && currElbo < prevElbo
                    fprintf(2, '[ERROR] ELBO decreased in iteration %d by %f!\n', it, abs(currElbo - prevElbo));
                end

                if RunConfig.getInstance().enableLogging
                    if ~isnan(prevElbo)
                        fprintf('------ ELBO increased by: %.4f\n', currElbo - prevElbo);
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

        % TODO: Add caching here?
        function value = get.alpha(obj)
            value = zeros(obj.K.Val, obj.M);
            for m = 1:obj.M
                value(:, m) = obj.views(m).alpha.E;
            end
        end

        % TODO: make private
        function value = computeW(obj)
            totalD = sum(obj.D);
            value = zeros(totalD, obj.K.Val);
        
            d = 0;
            for m = 1:obj.M
                Dm = obj.views(m).D;
                value(d + 1 : d + Dm, :) = obj.views(m).W.E;
                d = d + Dm;
            end
        end

        function value = computeTau(obj) 
            tau_values = cell(1, obj.M);

            for m = 1:obj.M
                tau_values{m} = obj.views(m).tau.E * ones(obj.views(m).D, 1);
            end
            
            value = vertcat(tau_values{:});
        end

        % TODO: make private + move the comment to the declaration section
        % Total variance including the noise and ALL views
        function value = computeTotalVar(obj)
            value = 0;
            for m = 1:obj.M
                view = obj.views(m);
                W_exp = view.W.E;
                tau_exp = view.tau.E;
                value = value + trace(W_exp * W_exp' + (1 / tau_exp) * eye(size(W_exp, 1)));
            end
        end

        % Total factor variance (not including noise), for all views and all factors
        function value = computeFactorsVar(obj)
            value = 0;
            for m = 1:obj.M
                view = obj.views(m);
                W_exp = view.W.E;
                value = value + trace(W_exp * W_exp');
            end
        end

        function value = computeNonRelVarWithin(obj)
            value = zeros(obj.M, obj.K.Val);
            for m = 1:obj.M
                W_exp = obj.views(m).W.E;
                for k = 1:obj.K.Val
                    value(m, k) = sum(W_exp(:, k).^2);
                end
            end
        end

        function value = computeRelVarWithin(obj)
            value = zeros(obj.M, obj.K.Val);
            vw = obj.varWithin;
            for m = 1:obj.M
                rowSum = sum(vw(m, :));
                if rowSum > 0
                    value(m, :) = vw(m, :) / rowSum * 100;
                end
            end
        end

        function value = computeVarWithin(obj)
            value = obj.nonRelVarWithin / obj.totalVar * 100;
        end

        function value = get.W(obj)
            if isempty(obj.W_)
                obj.W_ = obj.computeW();
            end
            value = obj.W_;
        end

        function value = get.tau(obj)
            if isempty(obj.tau_)
                obj.tau_ = obj.computeTau();
            end
            value = obj.tau_;
        end

        function set.W(obj, W)
            if size(W) ~= size(obj.W)
                CustomError.raiseError('InputCheck', 'Dims!'); % TODO!
            end
            obj.W_ = W;
        end

        function value = get.totalVar(obj)
            if isempty(obj.totalVar_)
                obj.totalVar_ = obj.computeTotalVar();
            end
            value = obj.totalVar_;
        end

        function set.totalVar(obj, totalVar)
            % TODO: Validate value
            obj.totalVar_ = totalVar;
        end

        function value = get.factorsVar(obj)
            if isempty(obj.factorsVar_)
                obj.factorsVar_ = obj.computeFactorsVar();
            end
            value = obj.factorsVar_;
        end

        function set.factorsVar(obj, factorsVar)
            % TODO: Validate value (dim at least)
            obj.factorsVar_ = factorsVar;
        end

        function value = get.nonRelVarWithin(obj)
            if isempty(obj.nonRelVarWithin_)
                obj.nonRelVarWithin_ = obj.computeNonRelVarWithin();
            end
            value = obj.nonRelVarWithin_;
        end

        function set.nonRelVarWithin(obj, nonRelVarWithin)
            % TODO: Validate value (dim at least)
            obj.nonRelVarWithin_ = nonRelVarWithin;
        end

        function value = get.relVarWithin(obj)
            if isempty(obj.relVarWithin_)
                obj.relVarWithin_ = obj.computeRelVarWithin();
            end
            value = obj.relVarWithin_;
        end

        function set.relVarWithin(obj, relVarWithin)
            % TODO: Validate value (dim at least)
            obj.relVarWithin_ = relVarWithin;
        end

        function value = get.varWithin(obj)
            if isempty(obj.varWithin_)
                obj.varWithin_ = obj.computeVarWithin();
            end
            value = obj.varWithin_;
        end

        function set.varWithin(obj, varWithin)
            % TODO: Validate value (dim at least)
            obj.varWithin_ = varWithin;
        end

        % TODO: Add desc. Sorts factors based on variances in descending order, and stores
        % then within object (factors) for clustering.
        function factors = getFactors(obj)
            factors = cell(1, obj.K.Val);
            varWithinFactors = sum(obj.varWithin, 1);
            [~, sortedOrder] = sort(varWithinFactors, 'descend');
            
            obj.W = obj.W(:, sortedOrder);
            obj.varWithin = obj.varWithin(:, sortedOrder);
            obj.relVarWithin = obj.relVarWithin(:, sortedOrder);
            
            for k = 1:size(obj.W, 2)
                factors{k} = obj.W(:, k);
            end
        end

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

                W_exp = view.W.E;
                mu_exp = view.mu.E;
                tau_exp = view.tau.E;

                sigma_Z = sigma_Z + W_exp' * tau_exp * W_exp;
                mu_Z = mu_Z + W_exp' * tau_exp * (inputViews{m} - mu_exp);
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
end