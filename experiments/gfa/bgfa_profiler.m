% Clear the workspace
close all; clear all; clc;

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
data = generateTwoViews();

X1 = data.X_tr{1}; % [N x D1];
X2 = data.X_tr{2}; % [N x D2]

profile on;

K = 10;
%                   data,   Mc, K, bound, maxIter, tol, doRotation 
bgfaModel = BGFA({X1', X2'}, 1, K, 'B');

return;
[elboVals, convIt] = bgfaModel.fit(10);


profile off;
profile viewer;

%%
%                             Z,         figureTitle,        figName,      folderName)
Visualization.plotLatentFactors(data.Z, 'True latent factors', 'trueLatentFactors', figsSubfolder);
Visualization.plotLatentFactors(bgfaModel.Z.E', 'Latent factors', 'latentFactors', figsSubfolder);
% TODO: reorder these!!! no '

%% Visualize loadings and alpha
totalD = sum(bgfaModel.D); % Total number of dimensions

trueW = zeros(totalD, data.trueK); % True K
estW = zeros(totalD, bgfaModel.K.Val);
estAlpha = zeros(bgfaModel.K.Val, bgfaModel.M);

d = 0;
for m = 1:bgfaModel.M
    Dm = bgfaModel.views(m).D;
    trueW(d + 1 : d + Dm, :) = data.W{m};
    estW(d + 1 : d + Dm, :) = bgfaModel.views(m).W.E;
    d = d + Dm;

    estAlpha(:, m) = bgfaModel.views(m).alpha.E;
end

%                            W,     dimList,   figTitle, figName, subfolderName
Visualization.plotLoadings(trueW, bgfaModel.D, 'True $\mathbf{W}$', 'trueW', figsSubfolder);
Visualization.plotLoadings(estW, bgfaModel.D, 'Estimated $\mathbf{W}$', 'estimatedW', figsSubfolder);

% alpha
Visualization.plotHintonDiagrams({-data.alpha', -estAlpha}, 'Lala');