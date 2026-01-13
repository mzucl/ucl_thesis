%% Prerequisites
% [1] Ensure that all required `.csv` files are located in the `datasets/hermes` directory.
% 
% Note: The folder already contains a `FoldsIdx` subdirectory with precomputed train/test indices 
% used to replicate the splits from the Python experiments. While the code to generate folds is 
% provided below in a commented-out block, it is not executed here. In a real experimental setting, 
% cross-validation folds are typically generated dynamically rather than fixed in advance.


%% Clear the workspace
close all; clearvars; clc;


%% Runtime configuration
rc = RunConfig.getInstance();
rc.inputValidation = false;
rc.enableLogging = false;


%% Setup logging
logFileName = ['logs/', mfilename, '.txt'];
if ~exist('logs', 'dir')
    mkdir('logs');
end
diary(logFileName); % start logging


%% Start timer
tic;


%% Experiment setup
numOfFolds = 10;
numModelSelectionRuns = 2;
K = 100;
baseDir = 'datasets/hermes/';


%% Initialize progress bar and results
numRuns = numOfFolds * numModelSelectionRuns;
h = waitbar(0, 'Progress...');
res = Results(numOfFolds); % obj


%% Define view filenames
% [NOTE] Output view always last
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
outputViewsCnt = 1;


%% Load all views
views = cell(size(viewFilenames));
M = numel(views);
for m = 1:M
    filePath = fullfile(baseDir, viewFilenames{m});
    views{m} = readmatrix(filePath, 'FileType', 'text')'; % TODO: set `featuresInCols` property
end
y = views{end}; % Labels

% cv = cvpartition(y, 'KFold', numOfFolds);

%% Initialize storage variables
% totalVars = nan(1, numOfFolds);
% factorsVars = nan(1, numOfFolds);
% varsWithin = cell(1, numOfFolds);
% relVarsWithin = cell(1, numOfFolds);
% factors = {};
% Ks = {};
% elbos = {};
Ws = {};
taus = {};


%% Cross-validation loop
for foldIdx = 1:numOfFolds
    % Get training and testing indices
    % trainIdx = cv.training(foldIdx);
    % testIdx = cv.test(foldIdx);

    % Load train/test indices
    filePath = fullfile(baseDir, sprintf('FoldsIdx/fold_%d', foldIdx - 1), 'train.csv');
    trainIdx = readmatrix(filePath, 'FileType', 'text');

    filePath = fullfile(baseDir, sprintf('FoldsIdx/fold_%d', foldIdx - 1), 'test.csv');
    testIdx = readmatrix(filePath, 'FileType', 'text');

    trainViews = cell(size(views));
    for v = 1:M
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

    % Model selection loop
    bestModel = NaN;
    bestElbo = -inf;
    bestIter = 0;
    for s = 1:numModelSelectionRuns
        model = BGFA(trainViews, M - 1, K, 'B', 10000, 1e-4); % Views are
        % expected in DxN format
        % model = SGFA(trainViews, K, 10000, 1e-4); % Views are expected in DxN format
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

    % Prediction & results
    [K_eff, predictionsTest] = bestModel.makePredictions(views, trainIdx, testIdx);
    % totalVars(foldIdx) = bestModel.totalVar;
    % factorsVars(foldIdx) = bestModel.factorsVar;
    % varsWithin{foldIdx} = bestModel.varWithin;
    % relVarsWithin{foldIdx} = bestModel.relVarWithin;
    % factors = [factors, bestModel.getFactors()];
    % Ks = [Ks, K_eff];
    % elbos = [elbos, bestElbo];
    Ws = [Ws, bestModel.W];
    taus = [taus, bestModel.tau];
    res.computeAndAppendMetrics(foldIdx, y(testIdx), predictionsTest, K_eff, bestIter, bestElbo);
end


%% Export results to CSV
summaryTable = res.computeMeanAndStd(LogicUtils.ternary(isa(model, 'BGFA'), 'BGFA', 'SGFA'), 'HERMES');  
writetable(summaryTable, [mfilename, '_summary', '.csv']);
resTable = res.storeResult();
writetable(resTable, [mfilename, '.csv']);


%% Print elapsed time and close progress bar
close(h);
elapsedTime = toc;
fprintf('The experiment took: %.4f seconds\n', elapsedTime);


%% Export `W` and `tau` to .npy files
if isa(model, 'SGFA')
    if ~exist('W', 'dir')
        mkdir('W');
    end
    if ~exist('tau', 'dir')
        mkdir('tau');
    end
    for foldIdx = 1:numel(Ws)
        py.numpy.save(sprintf('W/W_%d.npy', foldIdx - 1), Ws{foldIdx});
        py.numpy.save(sprintf('tau/tau_%d.npy', foldIdx - 1), taus{foldIdx});
    end
end





return;
%% Python visualisation (Optional)
% Caught unexpected exception of unknown type. in line 185 in
% visualisation.py file

% STEP 1
pyenv('Version', '/opt/anaconda3/envs/matlab_env/bin/python');

% Import the module if it's not automatically recognized
vis_mod = py.importlib.import_module('visualization');

% STEP 2
visualisation_instance = py.visualisation.Visualisation('python_plots');

% ...
py.help('visualisation');

factors_Ws = cell(1, numel(stableFactors));

for i = 1:numel(stableFactors)
    factors_Ws{i} = round(stableFactors{i}.allWs, 5);
end

% Convert to Python list
py_factors_Ws = py.list();
for i = 1:numel(factors_Ws)
    py_factors_Ws.append(py.numpy.array(factors_Ws{i}));
end

visualisation_instance.plot_stable_factors_matrix( ...
    py_factors_Ws, ...
    py.list(int32(bestModel.D)), ...
    py.list(prefixedLabels), ...
    true);