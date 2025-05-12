% Clear the workspace
close all; clearvars; clc;

rc = RunConfig.getInstance();
rc.inputValidation = true;
rc.enableLogging = true;

% Logging
logFileName = ['logs/', mfilename, '.txt'];
if ~exist('logs', 'dir')
    mkdir('logs');
end

% diary(logFileName); % start logging

% Start timing
tic;

% Setup
numOfFolds = 2;
numModelSelectionRuns = 5;
K = 100;
baseDir = 'datasets/hermes/';

% Progress bar
numRuns = numOfFolds * numModelSelectionRuns;
h = waitbar(0, 'Progress...');

res = Results(numOfFolds); % obj

% [NOTE] For prediction tasks, list input views before output views.
viewFilenames = {
    'rs_fMRI_ALFF.csv', ...
    'rs_fMRI_REHO.csv', ...
    'demographics.csv', ...
    'MRS.csv', ...
    'clinical.csv', ... 
    'cognitive.csv', ...
    'sMRI.csv', ...
    'nutrition.csv', ...
    'microbiome.csv', ...
    'classification.csv'
};

% [NOTE] For now this is always set to 1!
outputViewsCnt = 1;

% Initialize cell array to hold the data and read each file
views = cell(size(viewFilenames));
M = numel(views);
for m = 1:M
    filePath = fullfile(baseDir, viewFilenames{m});
    views{m} = readmatrix(filePath, 'FileType', 'text')'; % TODO: set `featuresInCols` property
end

% Labels
y = views{end};
% cv = cvpartition(y, 'KFold', numOfFolds);

% train_idx = readmatrix(fullfile(baseDir, ''), 'FileType', 'text');
% test_idx = readmatrix(filePath, 'FileType', 'text');

totalVars = nan(1, numOfFolds);
factorsVars = nan(1, numOfFolds);
varsWithin = cell(1, numOfFolds);
relVarsWithin = cell(1, numOfFolds);
factors = {};
Ks = {};
elbos = {};
Ws = {};
taus = {};

for foldIdx = 1:numOfFolds
    % Get training and testing indices
    % trainIdx = cv.training(foldIdx);
    % testIdx = cv.test(foldIdx);

    % Import train and test indices
    filePath = fullfile(baseDir, sprintf('FoldsIdx/fold_%d', foldIdx - 1), 'train.csv');
    trainIdx = readmatrix(filePath, 'FileType', 'text');

    filePath = fullfile(baseDir, sprintf('FoldsIdx/fold_%d', foldIdx - 1), 'test.csv');
    testIdx = readmatrix(filePath, 'FileType', 'text');

    % fprintf('Fold %d:\n', foldIdx);
    % fprintf('  Training label distribution:\n');
    % tabulate(y(trainIdx))
    % 
    % fprintf('  Test label distribution:\n');
    % tabulate(y(testIdx))
    % fprintf('\n');

    trainViews = cell(size(views));

    for v = 1:numel(views)
        trainViewData = views{v}(:, trainIdx);
        scaler = StandardScaler();

        scaler = scaler.fit(trainViewData);

        % Don't scale the last view
        if v ~= M
            trainViews{v} = scaler.transform(trainViewData);
            views{v} = scaler.transform(views{v});
        else 
            trainViews{v} = trainViewData;
        end
    end

    bestModel = NaN;
    bestElbo = -inf;
    bestIter = 0;

    for s = 1:numModelSelectionRuns
        % model = BGFA(views_train, M, K, 'B', 10000, 1e-4); % Views are expected in DxN
        model = SGFA(trainViews, K, 10000, 1e-4); % Views are expected in DxN
        [elboVals, iter] = model.fit(10);

        if elboVals(end) > bestElbo
            bestModel = model;
            bestElbo = elboVals(end);
            bestIter = iter;
        end
        % Progress bar
        k = (foldIdx - 1) * numModelSelectionRuns + s;
        waitbar(k/numRuns, h, sprintf('Progress: %d%%', round(100 * k/numRuns)));
        fprintf('\rProgress: %3.0f%%\n', 100 * k/numRuns);
    end

    [K_eff, predictionsTest] = bestModel.makePredictions(views, trainIdx, testIdx);

    % Compute variances
    totalVars(foldIdx) = bestModel.totalVar;
    factorsVars(foldIdx) = bestModel.factorsVar;
    varsWithin{foldIdx} = bestModel.varWithin;
    relVarsWithin{foldIdx} = bestModel.relVarWithin;

    factors = [factors, bestModel.getFactors()];
    Ks = [Ks, K_eff];
    elbos = [elbos, bestElbo];
    Ws = [Ws, bestModel.W];
    taus = [taus, bestModel.tau];

    res.computeAndAppendMetrics(foldIdx, y(testIdx), predictionsTest, K_eff, bestIter, bestElbo);
end

%% Export results
summaryTable = res.computeMeanAndStd('BGFA', 'HERMES');  
writetable(summaryTable, [mfilename, 'summary', '.csv']);

resTable = res.storeResult();
writetable(resTable, [mfilename, '.csv']);

%% Print summary 
close(h);
elapsedTime = toc;
fprintf('The experiment took: %.4f seconds\n', elapsedTime);

%% Cluster factors
function clusteringLabels = clusterFactors(factors, cosineSimilarityThreshold)
    factorsMatrix = cat(2, factors{:})';

    % Compute cosine distance matrix (1 - abs cosine similarity)
    cosSim = abs(factorsMatrix * factorsMatrix');  % Cosine similarity via dot product
    norms = sqrt(sum(factorsMatrix.^2, 2));
    cosSim = cosSim ./ (norms * norms');
    cosDist = 1 - cosSim;

    % Make sure diagonal is zero
    cosDist(1:size(cosDist,1)+1:end) = 0;
    
    % Run DBSCAN
    epsilon = 1 - cosineSimilarityThreshold;
    minPts = 1;
    clusteringLabels = dbscan(squareform(cosDist), epsilon, minPts, 'Distance', 'precomputed');
end

clusteringLabels = clusterFactors(factors, 0.8);


%%
factorsClusterLabels = clusteringLabels;

%%
function factorInfo = getFactorInfo(factorLabelsPerModel, uniqueFactorsLabelsPerModel, label, elbos, Ws, varsWithin, relVarsWithin)
    % Compiles information about a stable latent factor and returns a structure
    % with it, including the factor weights from the best model, all factor weights 
    % from all models that identified it, and the factors variances. 
    % All information is needed for visualizations.

    % Input:
    % - factorLabelsPerModel: Cell array of clustering labels for each factor identified by each model.
    % - uniqueFactorsLabelsPerModel: Cell array of unique clustering factor labels identified by each model.
    % - label: Label of the factor for which the information is compiled.

    % Output:
    % - factorInfo: Structure with factor information, including its weights from the best model, 
    %               all its weights from all models that identified it, and its variances.

    factorInfo = struct();
    
    % Get factor weight
    models = find(cellfun(@(x) any(x == label), uniqueFactorsLabelsPerModel)); % Get models that identified the factor
    elbos = cellfun(@(i) elbos{i}, num2cell(models)); % Extract ELBOs for these models
    [~, bestModelIdx] = max(elbos); % Find the best model (max ELBO)
    bestModelIdx = models(bestModelIdx); % Best model index
    factorIdx = find(factorLabelsPerModel{bestModelIdx} == label, 1); % Find factor index in the best model
    factorInfo.W = Ws{bestModelIdx}(:, factorIdx); % Factor weight from the best model
    
    % Get factor variances
    factorInfo.VarsWithin = varsWithin{bestModelIdx}(:, factorIdx);
    totalVar = sum(factorInfo.VarsWithin);
    factorInfo.VarsPercentage = floor((factorInfo.VarsWithin / totalVar * 100) * 100) / 100; % Percentage of variance explained by the factor
    factorInfo.VarsRelative = relVarsWithin{bestModelIdx}(:, factorIdx);
    
    % Get weights from all models
    allWs = cellfun(@(model) Ws{model}(:, find(factorLabelsPerModel{model} == label, 1)), ...
        num2cell(models), 'UniformOutput', false);

    % % Initialize a matrix to store the transformed data
    % allWs_matrix = zeros(10, 738);
    % 
    % % Loop through each cell in allWs and fill the rows of the matrix
    % for i = 1:numel(allWs)
    %     allWs_matrix(i, :) = allWs{i};  % Each cell contains a 738x1 column vector
    % end
    % 

    factorInfo.allWs = cell2mat(allWs)'; % Weights from all models that identified the factor
    
end



%% FROM HERE it is get_stable_factors(stabilityThreshold)
stabilityThreshold = 1;
% Initialize the factor model label array
factorModelLabel = [];
% TODO: Q: Why is Ks a cell array??
Ks = cell2mat(Ks);
for i = 1:length(Ks)
    k = Ks(i);
    factorModelLabel = [factorModelLabel, repmat(i, 1, k)];
end

% Convert factor_model_label to array
factorModelLabel = factorModelLabel(:);  % Make sure it's a column vector

% Split the factors cluster labels by model
factorLabelsPerModel = cell(1, length(Ks));
cumulativeSum = [0, cumsum(Ks)];
for i = 1:length(Ks)
    factorLabelsPerModel{i} = factorsClusterLabels(cumulativeSum(i)+1 : cumulativeSum(i+1));
end

% Get unique factor labels per model
uniqueFactorsLabelsPerModel = cellfun(@unique, factorLabelsPerModel, 'UniformOutput', false);

% Get stable factors labels
factorsPrevalence = histcounts(factorsClusterLabels, 'BinMethod', 'integers');
minCount = stabilityThreshold * numOfFolds;
stableFactorsLabels = find(factorsPrevalence >= minCount);

% Save stable factors weights and variances for visualizations
stableFactors = cell(1, length(stableFactorsLabels));
for i = 1:length(stableFactorsLabels)
    label = stableFactorsLabels(i);
    % Get stable factor information for detailed plots
    factorInfo = getFactorInfo(factorLabelsPerModel, uniqueFactorsLabelsPerModel, label, elbos, Ws, varsWithin, relVarsWithin);
    stableFactors{i} = factorInfo;
end

%% UNTIL HERE

% visualisation_instance.get_factors_matrix_plot_params()


%% Big figure

% % factors_Ws = cellfun(@(f) round(f.allWs, 5), stableFactors, 'UniformOutput', false);
% 
% 
% factors_Ws = cellfun(@(factor) round(factor.allWs, 5), stableFactors, 'UniformOutput', false);
% d = bestModel.D; % All models have the same D, maybe acccess this not from the model???
% featureLabels = prefixedLabels;
% 
% FVisualization.plotStableFactorsMatrix(factors_Ws, d, featureLabels, true);


%%
stableFactors(end) = []; % TODO: Delete, check what is wrong with last factor
factors_Ws = cell(1, numel(stableFactors));

for i = 1:numel(stableFactors)
    tmp = stableFactors{i}.allWs;
    if iscell(tmp)
        tmp = cell2mat(tmp')';  % Convert to 10x738 matrix
    end
    factors_Ws{i} = round(tmp, 5);
end


py_factors_Ws = py.list();
for i = 1:numel(factors_Ws)
    py_factors_Ws.append(py.numpy.array(factors_Ws{i}));
end


visualisation_instance.plot_stable_factors_matrix( ...
    py_factors_Ws, ...
    py.list(int32(bestModel.D)), ...
    py.list(prefixedLabels), ...
    true);
