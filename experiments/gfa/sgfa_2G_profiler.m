% Clear the workspace
close all; clearvars; clc;

rc = RunConfig.getInstance();
rc.inputValidation = false;
rc.enableLogging = false;

% Logging
logFileName = ['logs/', mfilename, '.txt'];
if ~exist('logs', 'dir')
    mkdir('logs');
end

% Start logging: only potential errors will be logged, since in the 
% profiler scripts we prefer to skip validation and logging to focus 
% solely on profiling the actual computational code.
diary(logFileName);

%% Generate data and train the model
data = Datasets.generateSyntheticGFAData(2);

profile on;

K = 10;
elboRecalcInterval = 10;
sgfaModel = SGFA(data.X_train, K);
[elboVals, convIt] = sgfaModel.fit(elboRecalcInterval);

profile off;
profile viewer;

%% Visualize true and inferred latent factors
Visualization.plotLatentFactors(data.Z, 'True latent factors', '', mfilename);
Visualization.plotLatentFactors(sgfaModel.Z.E, 'Inferred latent factors', '', mfilename);


%% Visualize loadings and alpha
Visualization.plotFactorLoadingsAndAlpha(data.W, data.D, data.alpha, 'bottom', '', 2.5, ...
    'True $\mathbf{W}^\top$ and \boldmath{$\alpha$}$^\top$');
Visualization.plotFactorLoadingsAndAlpha(sgfaModel.W, sgfaModel.D, sgfaModel.alpha, 'bottom', '', 2.5, ...
    'Inferred $\mathbf{W}^\top$ and \boldmath{$\alpha$}$^\top$');
