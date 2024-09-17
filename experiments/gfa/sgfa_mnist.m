% Clear the workspace
close all; clear all; clc;

% Logging
logFileName = 'logs/sgfa_mnist.txt';
if ~exist('logs', 'dir')
    mkdir('logs');
end

% Figures folder
figsSubfolder = 'sgfa_mnist';

diary(logFileName); % start logging

% Model settings
settings = ModelSettings.getInstance();
settings.VALIDATE = true;
settings.DEBUG = true;

folderName = 'mnist18'; % mnist38

%% Import data and train the model
csvread('datasets/mnist38/pixels.tsv');

X1 = csvread(['datasets/', folderName, '/pixels.tsv']); % [N x D1];
X2 = csvread(['datasets/', folderName, '/continuousLabels.tsv']); % [N x D2]; D2 = 1;

% profile on;


K = 10;
sgfaModel = SGFA({X1', X2'}, K);
[elboVals, convIt] = sgfaModel.fit(10);
return;

% profile off;
% profile viewer;

%% Predictive distribution
Z = sgfaModel.Z.E;
K = size(Z, 1);
W1 = sgfaModel.views(1).W.E;
W2 = sgfaModel.views(2).W.E;
mu1 = sgfaModel.views(1).mu.E;
mu2 = sgfaModel.views(2).mu.E;
T1 = sgfaModel.views(1).tau.E * eye(sgfaModel.D(1));
T2 = sgfaModel.views(2).tau.E * eye(sgfaModel.D(2));

sigma_Z = Utility.matrixInverse(eye(K) + W1' * T1 * W1);

idx = 13650;
xn1 = X1(idx, :)';
 

mu_zn = sigma_Z * (W1' * T1 * (xn1 - mu1));

predictedLabel = W2 * mu_zn + mu2;
disp(predictedLabel)

%%
%                             Z,         figureTitle,        figName,      folderName)
Visualization.plotLatentFactors(data.Z, 'True latent factors', 'trueLatentFactors', figsSubfolder);
Visualization.plotLatentFactors(sgfaModel.Z.E', 'Latent factors', 'latentFactors', figsSubfolder);
% TODO: reorder these!!! no '

%% Visualize loadings and alpha
totalD = sum(sgfaModel.D); % Total number of dimensions

trueW = zeros(totalD, data.trueK); % True K
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

%                            W,     dimList,   figTitle, figName, subfolderName
Visualization.plotLoadings(trueW, sgfaModel.D, 'True $\mathbf{W}$', 'trueW', figsSubfolder);
Visualization.plotLoadings(estW, sgfaModel.D, 'Estimated $\mathbf{W}$', 'estimatedW', figsSubfolder);

% alpha
Visualization.plotHintonDiagrams({-data.alpha', -estAlpha}, 'Lala');