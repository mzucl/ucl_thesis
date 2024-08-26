lala = true;
tic;
mu = Utility.generateRandomIntMatrix(50, 1);

for i = 1:100000
    E_Xt = Utility.ternary(lala, mu', mu);
end
elapsedTime = toc;
fprintf('Elapsed time : %.4f seconds\n', elapsedTime);


for i = 1:100000
    if lala
        E_Xt = mu';
    else
        E_Xt = mu;
    end
end
elapsedTime = toc;
fprintf('Elapsed time jkjkj : %.4f seconds\n', elapsedTime);


return;

% profile on;
% CODE
% profile off;

% profile viewer;

% % % Clear the workspace
% close all; clear; clc;
% % 
% % Add folders to MATLAB path
% addpath('src');
% addpath('models');
% addpath('helpers');
% addpath('tests');
% addpath('tests/src');
% addpath('experiments');
% addpath('figures');

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

% M = 9;
% N = 9;
% K = 9;
% 
% A = Utility.generateRandomIntMatrix(M, N);
% B = Utility.generateRandomIntMatrix(N, K);
% C = Utility.generateRandomIntMatrix(K, M);
% 
% tic;
% for i = 1:100000
%     trace(A * B);
% end
% elapsedTime = toc;
% fprintf('Elapsed time : %.4f seconds\n', elapsedTime);
% 
% tic;
% for i = 1:100000
%     A = A';
%     dot(A(:), B(:));
% end
% elapsedTime = toc;
% fprintf('Elapsed time : %.4f seconds\n', elapsedTime);
% 
% 
% % tic;
% b =  A * (B * ones(K, 1));
% elapsedTime = toc;
% fprintf('Elapsed time : %.4f seconds\n', elapsedTime);
% disp(b);

% disp(a == b);

% 
% tic;
% AB = (A * B)';
% a = dot(AB(:), C(:));
% elapsedTime = toc;
% fprintf('Elapsed time : %.4f seconds\n', elapsedTime);
% disp(a);
% 
% tic;
% BC = (B * C)';
% b = dot(BC(:), A(:));
% elapsedTime = toc;
% fprintf('Elapsed time : %.4f seconds\n', elapsedTime);
% disp(b);
% 
% tic;
% CA = (C * A)';
% c = dot(CA(:), B(:));
% elapsedTime = toc;
% fprintf('Elapsed time : %.4f seconds\n', elapsedTime);
% disp(c);
% 
% 
% tic;
% BC = B * C;
% A = A';
% b = dot(A(:), BC(:));
% elapsedTime = toc;
% fprintf('Elapsedvvv time : %.4f seconds\n', elapsedTime);
% disp(b);







