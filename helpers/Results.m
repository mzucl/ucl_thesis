classdef Results < handle
    properties
        specificity = [];
        bAcc = [];
        precision = [];
        recall = [];
        f1Score = [];
        npv = [];
        numOfEffFactors = [];
        numOfIters = [];
        elbo = [];
    end

    methods(Static, Access=private)
        function str = formatMeanAndStd(meanValue, stdValue, numOfDecimalSpaces)
            if nargin < 3
                numOfDecimalSpaces = 4;
            end
        
            formatSpec = sprintf('%%.%df ± %%.%df', numOfDecimalSpaces, numOfDecimalSpaces);
            str = sprintf(formatSpec, meanValue, stdValue);
        end
    end
    
    methods
        function obj = Results(numOfFolds)
            % Preallocate
            if nargin > 0
                obj.specificity = nan(1, numOfFolds); % Using 'nan' to easily detect is something is wrong when we proceed with calculations
                obj.bAcc = nan(1, numOfFolds);
                obj.precision = nan(1, numOfFolds);
                obj.recall = nan(1, numOfFolds);
                obj.f1Score = nan(1, numOfFolds);
                obj.npv = nan(1, numOfFolds);
                obj.numOfEffFactors = nan(1, numOfFolds);
                obj.numOfIters = nan(1, numOfFolds);
                obj.elbo = nan(1, numOfFolds);
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

            % Method to append values to each property array
            obj.specificity(foldIdx) = specificity_;
            obj.bAcc(foldIdx) = bAcc_;
            obj.precision(foldIdx) = precision_;
            obj.recall(foldIdx) = recall_;
            obj.f1Score(foldIdx) = f1Score_;
            obj.npv(foldIdx) = npv_;

            obj.numOfEffFactors(foldIdx) = numOfEffFactors;
            obj.numOfIters(foldIdx) = numOfIters;
            obj.elbo(foldIdx) = elbo;
        end

        function resTable = computeMeanAndStd(obj, modelName, datasetName)
            header = {'Model', 'Dataset', 'No. of LF', 'bAcc', 'Specificity', 'Precision', 'Recall', 'F1 Score', 'NPV', 'Iterations', 'ELBO'};
            resVals = {
                modelName, ...
                datasetName, ...
                Results.formatMeanAndStd(mean(obj.numOfEffFactors), std(obj.numOfEffFactors)), ...
                Results.formatMeanAndStd(mean(obj.bAcc), std(obj.bAcc)), ...
                Results.formatMeanAndStd(mean(obj.specificity), std(obj.specificity)), ...
                Results.formatMeanAndStd(mean(obj.precision), std(obj.precision)), ...
                Results.formatMeanAndStd(mean(obj.recall), std(obj.recall)), ...
                Results.formatMeanAndStd(mean(obj.f1Score), std(obj.f1Score)), ...
                Results.formatMeanAndStd(mean(obj.npv), std(obj.npv)), ...
                Results.formatMeanAndStd(mean(obj.numOfIters), std(obj.numOfIters)), ...
                Results.formatMeanAndStd(mean(obj.elbo), std(obj.elbo))
            };

            resTable = array2table(resVals, 'VariableNames', header);
        end
    end
end
