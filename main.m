% Clear the workspace
close all; clear; clc;

% Add folders to MATLAB path
addpath('src');
addpath('models');
addpath('helpers');
addpath('tests');

% Uncomment to run tests
testResults = runtests('tests');


% % Generate data
% numPoints = 100;
% dim = 10;
% stdDevs = [0.5 2 0.5 2 0.5 0.5 0.5 0.5 2 0.5];
% 
% data = generateSyntheticData(numPoints, dim, stdDevs);
% 
% K = 2; % dim - 1, because BPCA can infer the number of components
% obj = BayesianPCA(data', K); % Expects data in DxN format



% Usage:
% % X = randn(100, 5); % Example data
% % K = 2; % Number of principal components
% % model = BayesianPCA(X, K);
% % model = model.fit();
% 
% % numPoints = 100;
% % d = 10;
% % stdDevs = [0.5 2 0.5 2 0.5 0.5 0.5 0.5 2 0.5];
% % 
% % try
% %     data = generateSyntheticData(numPoints, d, stdDevs);
% % 
% %     W = PPCA(data, 9);
% % 
% %     hintonDiagram(W);
% % catch e
% %     % Handle errors
% %     disp('An error occurred:');
% %     disp(e.message);
% % end
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