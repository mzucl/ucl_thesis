% Clear the workspace
close all; clear; clc;

% Add folders to MATLAB path
addpath('src');
addpath('models');
addpath('helpers');
addpath('tests');

% Uncomment to run tests
% testResults = runtests('tests');

%% Generate the toy dataset
N = 100;
D = 10;
K = 3;
noiseVariance = 0.05;

[X, Z, W] = Utility.generateToyDataset(N, D, K, noiseVariance); % X is DxN

% X = X - mean(X, 2);

%% PPCA
W_PPCA = PPCA(X', D - 1); % PPCA expects X in NxD format



%% BPCA
K = D - 1; % dim - 1, because BPCA can infer effective number of components
numIter = 50;
obj = BayesianPCA(X, K, numIter); % BayesianPCA constructor expects data in DxN format
[elboVals, convIt, resArr] = obj.fit();



%% Visualization of the results
figure
subplot(1, 3, 1);
Utility.hintonDiagram(W, 'Ground truth');
hold on

subplot(1, 3, 2);
Utility.hintonDiagram(W_PPCA, "PPCA");
hold on;

subplot(1, 3, 3);
Utility.hintonDiagram(obj.W.EC, 'BPCA');
hold on;


Utility.plotStructVariables(resArr);