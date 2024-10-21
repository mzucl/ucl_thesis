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

folderName = 'mnist38'; % mnist38

%% Import data and train the model
[X, y] = Datasets.getMNISTData(folderName, false, 1000);












[X1_train, X2_train, X1_test, X2_test] = Datasets.trainTestSplitMNIST(folderName, false);

% profile on;

K = 50;
model = SGFA({X1_train', X2_train'}, K);
[elboVals, convIt] = model.fit(10);

% profile off;
% profile viewer;

%% Predictive distribution for the whole test dataset
Z = model.Z.E;
K = size(Z, 1);
W1 = model.views(1).W.E;
W2 = model.views(2).W.E;
mu1 = model.views(1).mu.E;
mu2 = model.views(2).mu.E;
T1 = model.views(1).tau.E * eye(model.D(1));
T2 = model.views(2).tau.E * eye(model.D(2));

sigma_Z = Utility.matrixInverse(eye(K) + W1' * T1 * W1);

MU_Z = sigma_Z * (W1' * T1 * (X1_test' - mu1));

predictedLabels = W2 * MU_Z + mu2;
predictedLabels = predictedLabels';

%%
predictedLabels(predictedLabels >= 0) = 1;  % Set elements greater than 0 to 1
predictedLabels(predictedLabels < 0) = -1; % Set elements less than 0 to -1

same_elements = predictedLabels == X2_test;
    
% Count the number of same elements
num_same_elements = sum(same_elements);
disp(num_same_elements * 100/size(X2_test, 1));

% Create the confusion matrix
confMatrix = confusionmat(X2_test, predictedLabels);
disp(confMatrix);