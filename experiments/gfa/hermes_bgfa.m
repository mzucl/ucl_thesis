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
numOfFolds = 1;
numModelSelectionRuns = 2;
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

clusteringLabels = clusterFactors(factors, 0.99);