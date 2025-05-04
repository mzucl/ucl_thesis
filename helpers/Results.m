classdef Results < handle
    properties
        elbo = [];
        bAcc = [];
        % rocAuc = [];
        precision = []; % Positive Predictive Value (PPV)
        recall = []; % Sensitivity
        f1Score = [];
        npv = [];
        specificity = [];
        numOfEffFactors = [];
        numOfIters = [];
    end

    methods(Static, Access=private)
        function str = formatMeanAndStd(meanValue, stdValue, numOfDecimalSpaces)
            if nargin < 3
                numOfDecimalSpaces = 4;% Utility.getConfigValue('Export', 'NUMBER_OF_DECIMAL_PLACES_MEAN_STD');
            end
        
            formatSpec = sprintf('%%.%df ± %%.%df', numOfDecimalSpaces, numOfDecimalSpaces);
            str = sprintf(formatSpec, meanValue, stdValue);
        end
    end
    
    methods
        function obj = Results(numOfFolds)
            % Preallocate
            % [NOTE] Using `nan` to easily detect if something goes wrong during subsequent calculations
            if nargin > 0
                obj.elbo = nan(1, numOfFolds);
                obj.bAcc = nan(1, numOfFolds);
                % obj.rocAuc = nan(1, numOfFolds);
                obj.precision = nan(1, numOfFolds);
                obj.recall = nan(1, numOfFolds);
                obj.f1Score = nan(1, numOfFolds);
                obj.npv = nan(1, numOfFolds);
                obj.specificity = nan(1, numOfFolds);
                obj.numOfEffFactors = nan(1, numOfFolds);
                obj.numOfIters = nan(1, numOfFolds);
            end
        end
        
        function obj = computeAndAppendMetrics(obj, foldIdx, labels, predictions, numOfEffFactors, numOfIters, elbo)
            % Compute metrics
            confusionMatrix = confusionmat(labels, predictions);
            
            TN = confusionMatrix(1,1);
            FP = confusionMatrix(1,2);
            FN = confusionMatrix(2,1);
            TP = confusionMatrix(2,2);
        
            sensitivity_ = TP / (TP + FN);     % TPR / recall
            specificity_ = TN / (TN + FP);     % TNR
        
            bAcc_ = (sensitivity_ + specificity_) / 2;
            
            precision_ = TP / (TP + FP);
            
            recall_ = sensitivity_;
            
            f1Score_ = 2 * (precision_ * recall_) / (precision_ + recall_);
            
            % Negative Predictive Value (NPV)
            npv_ = TN / (TN + FN);

            obj.elbo(foldIdx) = elbo;
            obj.bAcc(foldIdx) = bAcc_;
            obj.precision(foldIdx) = precision_;
            obj.recall(foldIdx) = recall_;
            obj.f1Score(foldIdx) = f1Score_;
            obj.npv(foldIdx) = npv_;
            obj.specificity(foldIdx) = specificity_;
            obj.numOfEffFactors(foldIdx) = numOfEffFactors;
            obj.numOfIters(foldIdx) = numOfIters;
            
        end

        function summaryTable = computeMeanAndStd(obj, modelName, datasetName)
            header = {'Model', 'Dataset', 'No. of LF', 'bAcc', 'Precision', 'Recall', 'F1 Score', 'NPV', 'Specificity', 'No. of Iterations', 'ELBO'};
            resVals = {
                modelName, ...
                datasetName, ...
                Results.formatMeanAndStd(mean(obj.numOfEffFactors), std(obj.numOfEffFactors)), ...
                Results.formatMeanAndStd(mean(obj.bAcc), std(obj.bAcc)), ...
                Results.formatMeanAndStd(mean(obj.precision), std(obj.precision)), ...
                Results.formatMeanAndStd(mean(obj.recall), std(obj.recall)), ...
                Results.formatMeanAndStd(mean(obj.f1Score), std(obj.f1Score)), ...
                Results.formatMeanAndStd(mean(obj.npv), std(obj.npv)), ...
                Results.formatMeanAndStd(mean(obj.specificity), std(obj.specificity)), ...
                Results.formatMeanAndStd(mean(obj.numOfIters), std(obj.numOfIters)), ...
                Results.formatMeanAndStd(mean(obj.elbo), std(obj.elbo))
            };

            summaryTable = array2table(resVals, 'VariableNames', header);
        end

        function resTable = storeResult(obj, numberOfDecimalPlaces)
            if nargin < 2
                numberOfDecimalPlaces = 6;
            end
            modelIdx = 1:numel(obj.bAcc);
            
            formatValue = @(x, dec) sprintf(['%.' num2str(dec) 'f'], x); % Formating function
            
            % Format each numeric column with the specified decimal places
            bAccFormatted = arrayfun(@(x) formatValue(x, numberOfDecimalPlaces), obj.bAcc, 'UniformOutput', false);
            precisionFormatted = arrayfun(@(x) formatValue(x, numberOfDecimalPlaces), obj.precision, 'UniformOutput', false);
            recallFormatted = arrayfun(@(x) formatValue(x, numberOfDecimalPlaces), obj.recall, 'UniformOutput', false);
            f1ScoreFormatted = arrayfun(@(x) formatValue(x, numberOfDecimalPlaces), obj.f1Score, 'UniformOutput', false);
            npvFormatted = arrayfun(@(x) formatValue(x, numberOfDecimalPlaces), obj.npv, 'UniformOutput', false);
            specificityFormatted = arrayfun(@(x) formatValue(x, numberOfDecimalPlaces), obj.specificity, 'UniformOutput', false);
            elboFormatted = arrayfun(@(x) formatValue(x, numberOfDecimalPlaces), obj.elbo, 'UniformOutput', false);
            
            % Create the table with formatted values
            header = {'Model', 'No. of LF', 'bAcc', 'Precision', 'Recall', 'F1 Score', ...
                'NPV', 'Specificity', 'No. of Iterations', 'ELBO'};
            
            resTable = table(modelIdx', obj.numOfEffFactors', bAccFormatted', precisionFormatted', recallFormatted', ...
                f1ScoreFormatted', npvFormatted', specificityFormatted', obj.numOfIters', elboFormatted', ...
                'VariableNames', header);
        end
    end
end
