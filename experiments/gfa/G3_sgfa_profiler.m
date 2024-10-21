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
data = Datasets.generateGFA_3G();

X1 = data.X_tr{1}; % [D1 x N];
X2 = data.X_tr{2}; % [D2 x N]
X3 = data.X_tr{3}; % [D3 x N]

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
totalD = sum(sgfaModel.D); % Total number of dimensions

trueW = zeros(totalD, data.K); % True K
estW = zeros(totalD, sgfaModel.K.Val);
estAlpha = zeros(sgfaModel.K.Val, sgfaModel.M);

d = 0;
for m = 1:sgfaModel.M
    Dm = sgfaModel.views(m).D;
    trueW(d + 1 : d + Dm, :) = data.W{m};
    estW(d + 1 : d + Dm, :) = sgfaModel.views(m).W.E;
    d = d + Dm;

    estAlpha(:, m) = sgfaModel.views(m).alpha.E;
end
% 
% %                            W,     dimList,   figTitle, figName, subfolderName
% Visualization.plotLoadings(trueW, sgfaModel.D, 'True $\mathbf{W}$', 'trueW', figsSubfolder);
% Visualization.plotLoadings(estW, sgfaModel.D, 'Estimated $\mathbf{W}$', 'estimatedW', figsSubfolder);
% 
% % alpha
% Visualization.plotHintonDiagrams({-data.alpha, -estAlpha}, 'alpha');


%% W and alpha
% Visualization.plotLoadingsAndAlpha(trueW, sgfaModel.D, -data.alpha, 'True $\mathbf{W}$ and \boldmath$\alpha$', 'trueW', figsSubfolder);
% Visualization.plotLoadingsAndAlpha(estW, sgfaModel.D, -estAlpha, 'Estimated $\mathbf{W}$ and \boldmath$\alpha$', 'trueW', figsSubfolder);

Visualization.plotLoadingsAndAlpha(trueW, sgfaModel.D, data.alpha, 'True $\mathbf{W}^\top$ and \boldmath{$\alpha$}$^\top$', 'bottom');
Visualization.plotLoadingsAndAlpha(estW, sgfaModel.D, estAlpha, 'Inferred $\mathbf{W}^\top$ and \boldmath{$\alpha$}$^\top$', 'top');