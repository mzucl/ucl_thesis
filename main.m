% 
% clear; close all;
% 
% d = 100;
% beta = 1e-1;
% X = rand(1,d);
% w = randn;
% b = randn;
% t = w'*X+b+beta*randn(1,d);
% x = linspace(min(X),max(X),d);   % test data
% 
% [model,llh] = linRegVb(X,t);
% % [model,llh] = rvmRegVb(X,t);
% % plot(llh);
% [y, sigma] = linRegPred(model,x,t);
% figure
% plotCurveBar(x,y,sigma);
% hold on;
% plot(X,t,'o');
% hold off



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
noiseVariance = 0.05;

% Generate the toy dataset
[X, Z, W] = Utility.generateToyDataset(N, D, K, noiseVariance); % X is DxN

W_PPCA = PPCA(X', D - 1); % PPCA expects X in NxD format
% -------------------------------------------------------- %

% BPCA
K = D - 1; % dim - 1, because BPCA can infer effective number of components
numIter = 5;
obj = BayesianPCA(X, K, numIter); % BayesianPCA constructor expects data in DxN format
[elboVals, convIt, resArr] = obj.fit();

figure;
subplot(1, 3, 1);
Utility.hintonDiagram(W, 'Ground truth');
hold on

subplot(1, 3, 2);
Utility.hintonDiagram(W_PPCA, "PPCA");
hold on;

subplot(1, 3, 3);
Utility.hintonDiagram(obj.W.ExpectationC, 'BPCA');
hold on;


Utility.plotStructVariables(resArr, 3);







% %% Code optimization TODOs
% % memory allocation (e.g. preallocating arrays)
% % matrix inverses (Choleskey)
% Value vs Handle classes (avoid copying)