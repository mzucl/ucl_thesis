% Clear the workspace
close all; clear all; clc;

% Logging
logFileName = 'logs/sgfa_3G.txt';
if ~exist('logs', 'dir')
    mkdir('logs');
end

% Figures folder
figsSubfolder = 'sgfa3G';

diary(logFileName); % start logging

% Model settings
settings = ModelSettings.getInstance();
settings.VALIDATE = false;
settings.DEBUG = false;


%% Generate data and train the model
data = Datasets.generateSyntheticGFAData(3);

X1 = data.X_train{1}; % [D1 x N]
X2 = data.X_train{2}; % [D2 x N]
X3 = data.X_train{3}; % [D3 x N]

% profile on;

K = 10;
sgfaModel = SGFA({X1, X2, X3}, K);
[elboVals, convIt] = sgfaModel.fit(10);

% 
% profile off;
% profile viewer;

%% Visualize latent factors
%                             Z,         figureTitle,        figName,      folderName)
Visualization.plotLatentFactors(data.Z, 'True latent factors', 'trueLatentFactors', figsSubfolder);
Visualization.plotLatentFactors(sgfaModel.Z.E, 'Inferred latent factors', 'latentFactors', figsSubfolder);

%% Visualize loadings and alpha
Visualization.plotFactorLoadingsAndAlpha(data.W, data.D, data.alpha, 'bottom', '', 3.5, 'True $\mathbf{W}^\top$ and \boldmath{$\alpha$}$^\top$');
Visualization.plotFactorLoadingsAndAlpha(sgfaModel.W, sgfaModel.D, sgfaModel.alpha, 'bottom', '', 3.5, 'Inferred $\mathbf{W}^\top$ and \boldmath{$\alpha$}$^\top$');
