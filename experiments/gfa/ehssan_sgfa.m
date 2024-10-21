% Clear the workspace
close all; clear all; clc;

% Start timing
tic;

% Logging
logFileName = ['logs/', mfilename, '.txt'];
if ~exist('logs', 'dir')
    mkdir('logs');
end
% diary(logFileName);

% Model settings
settings = ModelSettings.getInstance();
settings.VALIDATE = false;
settings.DEBUG = true;

% Setup
figsSubfolder = 'ehssan';
numOfFolds = 5;
stabilityRuns = 5;
K = 200;

res = Results(numOfFolds); % obj
resTable = table(); % table

% X
viewFileNames = {'connectivity__.csv'};

% Labels
y = readmatrix('datasets/ehssan/demographics__.csv', 'FileType', 'text');

for i = 1:length(viewFileNames)
    X = readmatrix(['datasets/ehssan/', viewFileNames{i}], 'FileType', 'text');

    cv = cvpartition(y, 'KFold', numOfFolds);

    for foldIdx = 1:numOfFolds
        % [X_tr, y_tr, X_te, y_te] = Datasets.trainTestSplit(X, y, true);

        % Get training and testing indices for fold 'foldIdx'
        trainIdx = cv.training(foldIdx);
        testIdx = cv.test(foldIdx);

        % Train
        X_tr = X(trainIdx, :); y_tr = y(trainIdx); 
        % Test
        X_te = X(testIdx, :); y_te = y(testIdx); 


        bestModel = NaN;
        bestElbo = -inf;
        bestIter = 0;

        for s = 1:stabilityRuns
            model = SGFA({X_tr', y_tr'}, K, 10000, 1e-4); % Views are expected in DxN
            [elboVals, iter] = model.fit(10);

            if elboVals(end) > bestElbo
                bestModel = model;
                bestElbo = elboVals(end);
                bestIter = iter;
            end
        end

        [K_eff, predictions_te] = bestModel.makePredictions(X_tr, y_tr, X_te);

        res.computeAndAppendMetrics(foldIdx, y_te, predictions_te, K_eff, bestIter, bestElbo);
    end

    resTable = [resTable; res.computeMeanAndStd('cGFA', viewFileNames{i})];
    
end
writetable(resTable, 'lala_v2.csv');

elapsedTime = toc;
fprintf('The experiment took: %.4f seconds\n', elapsedTime);