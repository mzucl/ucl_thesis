%% Clear the workspace
% Use clearvars when you just want to reset variables without wiping the whole environment. 
% Use clear all when you want a completely clean slate.

close all; clearvars; clc;



%% Add folders to MATLAB path
addpath(genpath('datasets'))
addpath('demos');
addpath(genpath('experiments'))
addpath(genpath('figures'))
addpath(genpath('helpers'));
addpath('logs');
addpath(genpath('models'));
addpath('plots');
addpath('src');
addpath(genpath('tests'))



%% Run tests
testResults = runtests('tests');



%% 
X1 = Utility.importCSV('datasets/mnist38/pixels.csv');
X2 = Utility.importCSV('datasets/mnist38/continuousLabels.csv');


%%

hfig = confusionchart(X2_test, predictedLabel, ...
                            'Title', 'Confusion Matrix', ...
                            'DiagonalColor', Constants.GREEN, ...
                            'OffDiagonalColor', Constants.RED, ...
                            'RowSummary', 'row-normalized', ...
                            'ColumnSummary', 'column-normalized');
Visualization.formatFigure(hfig)
%%





return;


A = Utility.generateRandomIntMatrix(50, 6000);
B = Utility.generateRandomIntMatrix(6000, 50);
C = Utility.generateRandomIntMatrix(50, 50);

tic;
for i = 1:1000
    trA = trace(A * B * C);
end
elapsedTime = toc;
fprintf('Elapsed time ABC: %.4f seconds\n', elapsedTime);


tic;
for i = 1:1
    trB = trace(B * C * A);
end
elapsedTime = toc;
fprintf('Elapsed time BCA (skipped): %.4f seconds\n', elapsedTime);

tic;
for i = 1:1000
    trC = trace(C * A * B);
end
elapsedTime = toc;
fprintf('Elapsed time CAB: %.4f seconds\n', elapsedTime);


tic;
for i = 1:1000
    D = A * B;
    trA = (D(:))' * C(:);
end
elapsedTime = toc;
fprintf('Elapsed time (AB)C: %.4f seconds\n', elapsedTime);

tic;
for i = 1:10000
    D = A * B;
    trA = C(:)' * D(:);
end
elapsedTime = toc;
fprintf('Elapsed time (AB)C: %.4f seconds\n', elapsedTime);

tic;
for i = 1:1000
    D = C*A;
    trB = (D(:))' * B(:);
end
elapsedTime = toc;
fprintf('Elapsed time (CA)B: %.4f seconds\n', elapsedTime);


tic;
for i = 1:1000
    D = C*A;
    trB = B(:)' * (D(:));
end
elapsedTime = toc;
fprintf('Elapsed time (CA)B: %.4f seconds\n', elapsedTime);




return;


% Create a figure
figure;

% Set the figure background color
set(gcf, 'Color', [1 1 1]); % White background

% Create some sample data and plot
x = 1:10;
y = rand(1, 10);
plot(x, y);

% Get current axes handle
ax = gca;


% Set the axes background color
set(ax, 'Color', Visualization.hexToRGB(Constants.BLUE)); % Light grey background
set(gcf, 'PaperPositionMode', 'auto');
% % Optionally, set the axes box to 'on' to see the background color clearly
% set(ax, 'Box', 'on');

% Ensure the figure's background color is preserved when saving
set(gcf, 'PaperPositionMode', 'auto');

% % Save the figure to a file (example: PNG format)
% print(gcf, 'my_plot.png', '-dpng', '-r300'); % Save with 300 dpi resolution

exportgraphics(gca, 'results.pdf', 'ContentType', 'vector');

exportgraphics(gca, 'my_plot.png');


return;

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







