% % Parameters
% N = 1000;  % number of samples
% D0 = 55;   % input features
% D1 = 3;    % output features
% 
% K = 2;     % common latent variables
% K0 = 3;    % first view's latent variables
% K1 = 3;    % second view's latent variables
% Kc = K + K0 + K1;  % total latent variables
% 
% % Generation of matrix W
% A0 = randn(D0, K);  % D0 x K matrix
% A1 = randn(D1, K);  % D1 x K matrix
% 
% B0 = randn(D0, K0); % D0 x K0 matrix
% B1 = randn(D1, K1); % D1 x K1 matrix
% 
% W0 = [A0, B0, zeros(D0, K1)];  % D0 x Kc matrix
% W1 = [A1, zeros(D1, K0), B1];  % D1 x Kc matrix
% W_tot = [W0; W1];              % (D0 + D1) x Kc matrix
% 
% % Generation of matrix Z
% Z = randn(N, Kc);  % N x Kc matrix
% 
% % Generation of matrix X
% X0 = Z * W0' + randn(N, D0) * 0.1;  % N x D0 matrix
% X1 = Z * W1' + randn(N, D1) * 0.1;  % N x D1 matrix




% Clear the workspace
% close all; clear; clc;

% Add folders to MATLAB path
% addpath('src');
% addpath('models');
% addpath('helpers');
% addpath('tests');

% Uncomment to run tests
% testResults = runtests('tests');

%% Generate the toy dataset
N = 100;
D = 10;
K = 10;
noiseVariance = 0.05;
% 
% [X, Z, W] = Utility.generateToyDataset(N, D, K, noiseVariance); % X is DxN
% [X2, Z2, W2] = Utility.generateToyDataset(N, D, K, noiseVariance); % X is DxN

profile on;
maxIter = 100;
obj =  GFA({X0', X1'}, K, maxIter);
[elboVals, convIt, resArr] = obj.fit();
profile off;

% View profiling results
profile viewer;

% %% PPCA
% W_PPCA = PPCA(X', D - 1); % PPCA expects X in NxD format
% 
% 
% 
% %% BPCA
% K = D - 1; % dim - 1, because BPCA can infer effective number of components
% numIter = 50;
% obj = BayesianPCA(X, K, numIter); % BayesianPCA constructor expects data in DxN format
% [elboVals, convIt, resArr] = obj.fit();
% 
% 
% 
% %% Visualization of the results
% figure
% subplot(1, 3, 1);
% Utility.hintonDiagram(W, 'Ground truth');
% hold on
% 
% subplot(1, 3, 2);
% Utility.hintonDiagram(W_PPCA, "PPCA");
% hold on;
% 
% subplot(1, 3, 3);
% Utility.hintonDiagram(obj.W.EC, 'BPCA');
% hold on;
% 
Utility.plotStructVariables(resArr);