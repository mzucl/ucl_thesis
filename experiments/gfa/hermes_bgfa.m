
clear all
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
stabilityRuns = 2;
K = 100;

res = Results(numOfFolds); % obj
resTable = table(); % table

% viewFileNames = {
%     % 'microbiome.csv', ...
%     'rs_fMRI_ALFF.csv', ...
%     'rs_fMRI_REHO.csv'
%     };


viewFileNames = {
    
   
    
    
    
    

    'rs_fMRI_ALFF.csv', ...
    'rs_fMRI_REHO.csv', ...
    'demographics.csv', ...
    'MRS.csv', ...
    'clinical.csv', ... 
    'cognitive.csv', ...
    'sMRI.csv', ...
    'nutrition.csv', ...
    'microbiome_1.csv', ...
    'microbiome_2.csv', ...
    'microbiome_3.csv', ...
    'microbiome_4.csv', ...

    % 'microbiome_5.csv', ...
    % 'microbiome_6.csv', ...
    % 'microbiome_7.csv', ...
    
    };

% Base directory
baseDir = 'datasets/hermes/';

% Number of input views
M_input = numel(viewFileNames);


% Initialize cell array to hold the data
views = cell(1, size(viewFileNames, 2) + 1);

% Read each file
for i = 1:M_input
    filePath = fullfile(baseDir, viewFileNames{i});
    views{i} = readmatrix(filePath, 'FileType', 'text')';
end

% Labels
y = readmatrix(fullfile(baseDir, 'classification.csv'), 'FileType', 'text')';

cv = cvpartition(y, 'KFold', numOfFolds);

for foldIdx = 1:numOfFolds
    % Get training and testing indices
    trainIdx = cv.training(foldIdx);
    testIdx = cv.test(foldIdx);

    % Split labels
    y_train = y(trainIdx);
    y_test = y(testIdx);

    % Split each view
    views_train = cell(size(views));
    views_test = cell(size(views));
    
    for v = 1:M_input
        views_train{v} = views{v}(:, trainIdx);
        views_test{v} = views{v}(:, testIdx);
    end

    % Add label view (to the train data only)
    views_train{M_input + 1} = y_train;

    bestModel = NaN;
    bestElbo = -inf;
    bestIter = 0;

    for s = 1:stabilityRuns
        % model = BGFA(views_train, M_input, K, 'B', 10000, 1e-4); % Views are expected in DxN
        model = SGFA(views_train, K, 10000, 1e-4); % Views are expected in DxN
        [elboVals, iter] = model.fit(10);

        if elboVals(end) > bestElbo
            bestModel = model;
            bestElbo = elboVals(end);
            bestIter = iter;
        end
    end

    % [K_eff, predictions_te] = bestModel.makePredictions(X_tr, y_tr, X_te);
    [K_eff, predictions_test] = bestModel.makePredictions(views_train, views_test);

    res.computeAndAppendMetrics(foldIdx, y_test, predictions_test, K_eff, bestIter, bestElbo);
end

resTable = [resTable; res.computeMeanAndStd('bGFA', viewFileNames{i})];  


writetable(resTable, [mfilename, '.csv']);

elapsedTime = toc;
fprintf('The experiment took: %.4f seconds\n', elapsedTime);







% 
% 
% 
% 
% % Clear the workspace
% close all; clearvars; clc;
% 
% rc = RunConfig.getInstance();
% rc.inputValidation = true;
% rc.enableLogging = true;
% 
% % Logging
% logFileName = ['logs/', mfilename, '.txt'];
% if ~exist('logs', 'dir')
%     mkdir('logs');
% end
% 
% % diary(logFileName); % start logging
% 
% % Start timing
% tic;
% 
% % Setup
% numOfFolds = 5;
% stabilityRuns = 10;
% K = 100;
% 
% res = Results(numOfFolds); % obj
% resTable = table(); % table
% 
% % X
% viewFileNames = {
%     'clinical.csv', ... 
%     'cognitive.csv', ...
%     'demographics.csv', ...
%     'microbiome.csv'
%     % 'MRS.csv', ...
%     % 'nutrition.csv', ...
%     % 'rs_fMRI_ALFF.csv', ...
%     % 'rs_fMRI_REHO.csv', ...
%     % 'sMRI.csv'
%     };
% 
% % Base directory
% baseDir = 'datasets/hermes/';
% 
% % Initialize cell array to hold the data
% views = cell(size(viewFileNames));
% 
% % Read each file
% for i = 1:numel(viewFileNames)
%     filePath = fullfile(baseDir, viewFileNames{i});
%     views{i} = readmatrix(filePath, 'FileType', 'text');
% end
% 
% % Labels
% y = readmatrix(fullfile(baseDir, 'classification.csv'), 'FileType', 'text');
% 
% cv = cvpartition(y, 'KFold', numOfFolds);
% 
% for foldIdx = 1:numOfFolds
%     % Get training and testing indices
%     trainIdx = cv.training(foldIdx);
%     testIdx = cv.test(foldIdx);
% 
%     % Split labels
%     y_train = y(trainIdx);
%     y_test = y(testIdx);
% 
%     % Split each view
%     views_train = cell(size(views));
%     views_test = cell(size(views));
% 
%     for v = 1:numel(views)
%         views_train{v} = views{v}(trainIdx, :);
%         views_test{v} = views{v}(testIdx, :);
%     end
% 
%     bestModel = NaN;
%     bestElbo = -inf;
%     bestIter = 0;
% 
%     for s = 1:stabilityRuns
%         % model = BGFA({views_train{1}', views_train{2}', views_train{3}', ...
%         %     views_train{4}', views_train{5}', views_train{6}', views_train{7}', ...
%         %     views_train{8}', views_train{9}', y_train'}, 1, K, 'B', 10000, 1e-4); % Views are expected in DxN
%         % 
%         model = SGFA({views_train{1}', views_train{2}', views_train{3}', ...
%             views_train{4}', y_train'}, K, 10000, 1e-4); % Views are expected in DxN
% 
% 
%         % model = BGFA({views_train{1}', views_train{2}', views_train{3}', ...
%         %     views_train{4}', y_train'}, 4, K, 'B', 10000, 1e-4); % Views are expected in DxN
% 
%         % model = SGFA({X_tr', y_tr'}, K, 10000, 1e-4); % Views are expected in DxN
%         [elboVals, iter] = model.fit(10);
% 
%         if elboVals(end) > bestElbo
%             bestModel = model;
%             bestElbo = elboVals(end);
%             bestIter = iter;
%         end
%     end
% 
%     [K_eff, predictions_te] = bestModel.makePredictions(views_train, y_train, views_views_test);
% 
%     res.computeAndAppendMetrics(foldIdx, y_test, predictions_te, K_eff, bestIter, bestElbo);
% end
% 
% resTable = [resTable; res.computeMeanAndStd('bGFA', viewFileNames{i})];  
% 
% writetable(resTable, [mfilename, '.csv']);
% 
% elapsedTime = toc;
% fprintf('The experiment took: %.4f seconds\n', elapsedTime);
