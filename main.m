% Clear the workspace
close all; clear; clc;

% Add folders to MATLAB path
addpath('src');
addpath('models');
addpath('helpers');
addpath('tests');

% Uncomment to run tests
% testResults = runtests('tests');

% Generate toy dataset
N = 100;
D = 10;
K = 3;
noiseVariance = 0.1;

% Generate the toy dataset
[X, Z, W] = Utility.generateToyDataset(N, D, K, noiseVariance);

W_PPCA = PPCA(X', 9); % PPCA expects X in NxD format
% -------------------------------------------------------- %

% BPCA
% K = D - 1; % dim - 1, because BPCA can infer effective number of components
% numIter = 20;
% obj = BayesianPCA(X, K, numIter); % BayesianPCA constructor expects data in DxN format
% [elboVals, it, resArr] = obj.fit();

% figure;
% subplot(1, 3, 1);
% hintonDiagram(W, 'Ground truth');
% disp(size(W));
% hold on

figure
% subplot(1, 3, 2);
Utility.hintonDiagram(W_PPCA, "PPCA");
% hold on;
% disp(size(W_PPCA));
% 
% subplot(1, 3, 3);
% hintonDiagram(obj.W.ExpectationC, 'BPCA');
% hold on;
% 
% disp(size(obj.W.ExpectationC));















% % Generate data
% numPoints = 300;
% dim = 10;
% stdDevs = [0.5, 2, 0.5, 2, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
% data = generateSyntheticData(numPoints, dim, stdDevs);
% 
% % PPCA
% % Uncomment to compare with PPCA
% % W = PPCA(data, 9);
% % hintonDiagram(W, "PPCA");
% % 
% % BPCA
% K = dim - 1; % dim - 1, because BPCA can infer effective number of components
% numIter = 10;
% obj = BayesianPCA(data', K, numIter); % BayesianPCA constructor expects data in DxN format
% [elboVals, it, resArr] = obj.fit();
% 
% % Create a figure
% Utility.plotStructVariables(resArr);
% 
% 
% 










% 
% figure
% hintonDiagram(obj.W.ExpectationC, "BPCA");

% X = randn(100, 5);
% K = 2;
% model = BayesianPCA(X, K);
% model = model.fit();
% 
% numPoints = 100;
% d = 10;
% stdDevs = [0.5 20 0.5 30 0.5 0.5 0.5 0.5 50 0.5];
% 
% try
%     data = generateSyntheticData(numPoints, d, stdDevs);
% 
%     W = PPCA(data, 9);
% 
%     hintonDiagram(W);
% catch e
%     % Handle errors
%     disp('An error occurred:');
%     disp(e.message);
% end
% 
% 
% 
% 
% % 
% % % Example usage
% % X = randn(100, 10); % Replace with your dataset
% % principalComponentsAndHinton(X);




%% Code optimization TODOs
% memory allocation (e.g. preallocating arrays)
% matrix inverses (Choleskey)
% Value vs Handle classes (avoid copying)