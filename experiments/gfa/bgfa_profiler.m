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

diary(logFileName); % start logging

%% Import data and train the model
folderName = 'mnist38'; 
%                                                                             binaryLabels
[X1_train, X2_train, X1_test, X2_test] = Datasets.trainTestSplitMNIST(folderName, true); 

% Speed up
% X1_train = X1_train(1:500, :);
% X2_train = X2_train(1:500, :);

% profile on;

K = 10;
%                   data,   Mc, K, bound, maxIter, tol, doRotation 
model = BGFA({X1_train', X2_train'}, 1, K, 'B');

%%
% return;
[elboVals, convIt] = model.fit();

% profile off;
% profile viewer;

%%
Z = model.Z.E;
W1 = model.views(1).W.E;
W2 = model.views(2).W.E;
mu1 = model.views(1).mu.E;
mu2 = model.views(2).mu.E;

%%
Z = model.Z.E;
K = size(Z, 1);
W1 = model.views(1).W.E;
W2 = model.views(2).W.E;
mu1 = model.views(1).mu.E;
mu2 = model.views(2).mu.E;
T1 = model.views(1).tau.E * eye(model.D(1));


sigma_Z = Utility.matrixInverse(eye(K) + W1' * T1 * W1);

MU_Z = sigma_Z * (W1' * T1 * (X1_test' - mu1));

predictedLabels = Bound.sigma(W2 * MU_Z + mu2);
predictedLabels = predictedLabels';
%%

predictedLabels(predictedLabels >= 0.5) = 1;  % Set elements greater than 0 to 1
predictedLabels(predictedLabels < 0.5) = 0; % Set elements less than 0 to -1