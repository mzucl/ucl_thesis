% Clear the workspace
close all; clear; clc;

% Add folders to MATLAB path
addpath('src');
addpath('models');
addpath('helpers');
addpath('tests');

% Uncomment to run tests
% testResults = runtests('tests');

% Generate data
numPoints = 300;
dim = 10;
stdDevs = [0.5, 2, 0.5, 2, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
data = generateSyntheticData(numPoints, dim, stdDevs);

% PPCA
% Uncomment to compare with PPCA
% W = PPCA(data, 9);
% hintonDiagram(W, "PPCA");
% 
% BPCA
K = dim - 1; % dim - 1, because BPCA can infer effective number of components
numIter = 10;
obj = BayesianPCA(data', K, numIter); % BayesianPCA constructor expects data in DxN format
[elboVals, convIter] = obj.fit();
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