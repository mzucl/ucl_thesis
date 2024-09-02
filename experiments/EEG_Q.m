%% Import data
data = readtable('control_eeg_questionnaire.csv');

colsToDrop = {'ppid', 'isControl', 'gender', 'age'};

data(:, colsToDrop) = [];

dataArray = table2array(data);

colNumEEG = 155;

X1 = dataArray(:, 1:colNumEEG); % [N x D1]
X2 = dataArray(:, colNumEEG+1:end); % [N x D2]

% Scale datasets
% X1 = Datasets.standardScaler(X1);
% X2 = Datasets.standardScaler(X2);

K = 50;

profile on;

gfaModel =  GFA({X1', X2'}, K);
[elboVals, convIt] = gfaModel.fit(10);

profile off;
profile viewer;


return;

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
expZ = bestModel.Z.E;

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
totalD = sum(bestModel.D); % Total number of dimensions

trueW = zeros(totalD, data.trueK); % True K
estW = zeros(totalD, bestModel.K.Val);
estAlpha = zeros(bestModel.K.Val, bestModel.M);

d = 0;
for m = 1:bestModel.M
    Dm = bestModel.views(m).D;
    trueW(d + 1 : d + Dm, :) = data.W{m};
    estW(d + 1 : d + Dm, :) = bestModel.views(m).W.E;
    d = d + Dm;

    estAlpha(:, m) = bestModel.views(m).alpha.E;
end

Visualization.plotLoadings(trueW, bestModel.D, 'True W');
Visualization.plotLoadings(estW, bestModel.D, 'Estimated W');

figure;
ax1 = subplot(1, 2, 1);
ax2 = subplot(1, 2, 2);

Visualization.hintonDiagram(-data.alpha, ax1, 'True -1 * alpha');
Visualization.hintonDiagram(-estAlpha, ax2, 'Estimated -1 * alpha');