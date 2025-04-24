% Clear the workspace
close all; clearvars; clc;

% Logging
logFileName = 'logs/sgfa_2G.txt';
if ~exist('logs', 'dir')
    mkdir('logs');
end

% Figures folder
figsSubfolder = 'sgfa';

diary(logFileName); % start logging

% Model settings
settings = ModelSettings.getInstance();
settings.VALIDATE = false;
settings.DEBUG = false;


%% Generate data and train the model
data = Datasets.generateSyntheticGFAData(2);

X1 = data.X_train{1}; % [D1 x N]
X2 = data.X_train{2}; % [D2 x N]

profile on;

K = 10;
sgfaModel = SGFA({X1, X2}, K);
[elboVals, convIt] = sgfaModel.fit(10);


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
