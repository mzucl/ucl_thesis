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
numOfFolds = 5;
numModelSelectionRuns = 10;
K = 100;
baseDir = 'datasets/aida/';

res = Results(numOfFolds); % obj
resTable = table(); % table

% [NOTE] For prediction tasks, list input views before output views.
viewFilenames = {
    'eeg_avg_neutral.csv', ...
    'eeg_avg_sad', ...
    'questionnaires.csv', ...
    'gender.csv'
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
cv = cvpartition(y, 'KFold', numOfFolds);

for foldIdx = 1:numOfFolds
    % Get training and testing indices
    trainIdx = cv.training(foldIdx);
    testIdx = cv.test(foldIdx);

    fprintf('Fold %d:\n', foldIdx);
    fprintf('  Training label distribution:\n');
    tabulate(y(trainIdx))

    fprintf('  Test label distribution:\n');
    tabulate(y(testIdx))
    fprintf('\n');

    % Split each view
    trainViews = cell(size(views));
    testViews = cell(size(views));
    
    for v = 1:numel(views)
        trainViews{v} = views{v}(:, trainIdx);
        testViews{v} = views{v}(:, testIdx);
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
    end

    [K_eff, predictionsTest] = bestModel.makePredictions(views, trainIdx, testIdx);

    res.computeAndAppendMetrics(foldIdx, y(testIdx), predictionsTest, K_eff, bestIter, bestElbo);
end

resTable = [resTable; res.computeMeanAndStd('bGFA', viewFilenames{m})];  


writetable(resTable, [mfilename, '.csv']);

elapsedTime = toc;
fprintf('The experiment took: %.4f seconds\n', elapsedTime);