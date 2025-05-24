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
numOfFolds = 10;
numModelSelectionRuns = 10;
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
        % model = BGFA(trainViews, M - 1, K, 'B', 10000, 1e-4); % Views are
        % expected in DxN format
        model = SGFA(trainViews, K, 10000, 1e-4); % Views are expected in DxN format
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
summaryTable = res.computeMeanAndStd('SGFA', 'HERMES');  
writetable(summaryTable, [mfilename, 'summary', '.csv']);

resTable = res.storeResult();
writetable(resTable, [mfilename, '.csv']);

%% Print summary 
close(h);
elapsedTime = toc;
fprintf('The experiment took: %.4f seconds\n', elapsedTime);


%% Export `W` and `tau` to `.npy` format
if ~exist('W', 'dir')
    mkdir('W');
end

if ~exist('tau', 'dir')
    mkdir('tau');
end

for foldIdx = 1:numel(Ws)
    % Save `W`
    filename = sprintf('W/W_%d.npy', foldIdx - 1);
    py.numpy.save(filename, Ws{foldIdx});  % Save the W matrix for the current fold

    % Save `tau`
    filename = sprintf('tau/tau_%d.npy', foldIdx - 1);
    py.numpy.save(filename, taus{foldIdx});  % Save the tau matrix for the current fold
end






return;
%% Execute Python code from MATLAB
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