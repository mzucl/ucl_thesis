% % Experiment 1
% pathTrueW1 = fullfile(pwd, 'figures', 'trueW1.png');
% pathTrueZ = fullfile(pwd, 'figures', 'trueZ.png');
% 
% pathW1 = fullfile(pwd, 'figures', 'W1.png');
% pathZ = fullfile(pwd, 'figures', 'Z.png');
% 
% % Ensure the 'figures' folder exists
% if ~exist(fullfile(pwd, 'figures'), 'dir')
%     mkdir(fullfile(pwd, 'figures'));
% end
% 
% addpath('figures');

data = get_data_2g();
D1 = size(data.X_tr{1}, 2);
% plot_loadings(data.W{1}, D1, pathTrueW1);

maxIter = 50;
K = 10;
X1 = data.X_tr{1};
X2 = data.X_tr{2};
obj =  GFA({X1', X2'}, K, maxIter);
[elboVals, convIt, resArr] = obj.fit();
Utility.plotStructVariables(resArr);

% plot_loadings(obj.views(1).W.EC, D1, pathW1);


% % Specify the path to save the figure
% W_path = fullfile(pwd, 'figures', 'loading_plot.png');
% 
% % Ensure the 'figures' folder exists
% if ~exist(fullfile(pwd, 'figures'), 'dir')
%     mkdir(fullfile(pwd, 'figures'));
% end
% 
% addpath('figures');
% 
% % Call the function
% % plot_loadings(res.W{1}, size(res.X_tr{1}, 2), W_path);
% Z_path = fullfile(pwd, 'figures', 'Z_path.png');
% plot_Z(res.Z, Z_path);


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
% N = 100;
% D = 10;
% K = 10;
% noiseVariance = 0.05;
% % 
% % [X, Z, W] = Utility.generateToyDataset(N, D, K, noiseVariance); % X is DxN
% % [X2, Z2, W2] = Utility.generateToyDataset(N, D, K, noiseVariance); % X is DxN

% profile on;
% maxIter = 100;
% X1 = res.X_tr{1};
% X2 = res.X_tr{2};
% obj =  GFA({X1', X2'}, K, maxIter);
% [elboVals, convIt, resArr] = obj.fit();
% profile off;

% View profiling results
% profile viewer;

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
% Utility.plotStructVariables(resArr);