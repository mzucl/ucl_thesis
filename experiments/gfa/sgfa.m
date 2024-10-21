% Clear the workspace
close all; clear all; clc;

% Logging
logFileName = 'logs/sgfa_2G.txt';
if ~exist('logs', 'dir')
    mkdir('logs');
end

diary(logFileName); % start logging

% Model settings
settings = ModelSettings.getInstance();
% settings.VALIDATE = false;
% settings.DEBUG = false;



%% Generate data and train the model
data = generateTwoViews();

X1 = data.X_tr{1}; % [N x D1];
X2 = data.X_tr{2}; % [N x D2]

% Scale datasets
% X1 = Datasets.standardScaler(X1);
% X2 = Datasets.standardScaler(X2);

K = 10;

stabilityRun = 20;
modelSelectionIter = 10;
convItAvg = 0;

tic;

for s = 1:stabilityRun 
    maxElbo = -Inf;
    bestW = NaN;
    convIt = NaN;

    for i = 1:modelSelectionIter
        sgfaModel = SGFA({X1', X2'}, K);
        [elboVals, it] = sgfaModel.fit(10);
    
        if elboVals(end) > maxElbo
            maxElbo = elboVals(end);
            sgfaModel = sgfaModel;
            convIt = it;
        end
    end

    convItAvg = convItAvg + convIt;
    disp(['The best model converged in ', num2str(convIt), 'iterations.\n']);
end

elapsedTime = toc;
fprintf('\n\n\nElapsed time: %.4f [s]\n', elapsedTime);
fprintf('Average number of iterations: %.4f\n', convItAvg / stabilityRun);

diary off; 

% return;
%% Visualize true and recovered latent factors
% True factors
trueNumOfFactors = size(data.Z, 2);

figure;
for i = 1:trueNumOfFactors
    subplot(trueNumOfFactors, 1, i);
    factor = data.Z(:, i);
    plot(factor, '.', 'MarkerSize', 4);
    hold on;
end
sgtitle('True latent factors');

% Recovered latent factors
expZ = sgfaModel.Z.E;

% Effective number of factors
numEffFactors = size(expZ, 1);

figure;
for i = 1:numEffFactors
    subplot(numEffFactors, 1, i);
    factor = expZ(i, :);
    plot(factor, '.', 'MarkerSize', 4);
    hold on;
end
sgtitle('Latent factors');



%% Visualize loadings and alpha
totalD = sum(sgfaModel.D); % Total number of dimensions

trueW = zeros(totalD, data.trueK); % True K
estW = zeros(totalD, sgfaModel.K.Val);
estAlpha = zeros(sgfaModel.K.Val, sgfaModel.M);

d = 0;
for m = 1:sgfaModel.M
    Dm = sgfaModel.views(m).D;
    trueW(d + 1 : d + Dm, :) = data.W{m};
    estW(d + 1 : d + Dm, :) = sgfaModel.views(m).W.E;
    d = d + Dm;

    estAlpha(:, m) = sgfaModel.views(m).alpha.E;
end

Visualization.plotLoadings(trueW, sgfaModel.D, 'True W');
Visualization.plotLoadings(estW, sgfaModel.D, 'Estimated W');

figure;
ax1 = subplot(1, 2, 1);
ax2 = subplot(1, 2, 2);

Visualization.hintonDiagram(-data.alpha, ax1, 'True -1 * alpha');
Visualization.hintonDiagram(-estAlpha, ax2, 'Estimated -1 * alpha');