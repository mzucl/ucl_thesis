% profile on;
% CODE
% profile off;

% profile viewer;

% % Clear the workspace
close all; clear; clc;
% 
% % Add folders to MATLAB path
addpath('src');
addpath('models');
addpath('helpers');
addpath('tests');
addpath('experiments');
addpath('figures');

% Uncomment to run tests
% testResults = runtests('tests');

%% 


% % % Experiment 1
% % pathTrueW1 = fullfile(pwd, 'figures', 'trueW1.png');
% % pathTrueZ = fullfile(pwd, 'figures', 'trueZ.png');
% % 
% % pathW1 = fullfile(pwd, 'figures', 'W1.png');
% % pathZ = fullfile(pwd, 'figures', 'Z.png');
% % 
% % % Ensure the 'figures' folder exists
% % if ~exist(fullfile(pwd, 'figures'), 'dir')
% %     mkdir(fullfile(pwd, 'figures'));
% % end
% % 
% % addpath('figures');
% 
% data = get_data_2g();
% D1 = size(data.X_tr{1}, 2);
% % plot_loadings(data.W{1}, D1, pathTrueW1);
% 
% maxIter = 500;
% K = 10;
% X1 = data.X_tr{1};
% X2 = data.X_tr{2};
% obj =  GFA({X1', X2'}, K, maxIter);
% [elboVals, convIt, resArr] = obj.fit();
% Utility.plotStructVariables(resArr);
% 
% % plot_loadings(obj.views(1).W.EC, D1, pathW1);
% 
% 
% % % Specify the path to save the figure
% % W_path = fullfile(pwd, 'figures', 'loading_plot.png');
% % 
% % % Ensure the 'figures' folder exists
% % if ~exist(fullfile(pwd, 'figures'), 'dir')
% %     mkdir(fullfile(pwd, 'figures'));
% % end
% % 
% % addpath('figures');
% % 
% % % Call the function
% % % plot_loadings(res.W{1}, size(res.X_tr{1}, 2), W_path);
% % Z_path = fullfile(pwd, 'figures', 'Z_path.png');
% % plot_Z(res.Z, Z_path);
% 
% 










% maxIter = 100;
% X1 = res.X_tr{1};
% X2 = res.X_tr{2};
% obj =  GFA({X1', X2'}, K, maxIter);
% [elboVals, convIt, resArr] = obj.fit();

