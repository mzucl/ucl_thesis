%% Import data
data = readtable('control_eeg_questionnaire.csv');
colsToDrop = {'ppid', 'gender', 'age'};
data(:, colsToDrop) = [];
data = table2array(data);

% Extract views
colNumEEG = 155; % Number of features in EEG dataset
X1 = data(:, 2:colNumEEG + 1); % [N x D1]
X2 = data(:, colNumEEG + 2:end); % [N x D2]
X3 = data(:, 1); % 'isControl'; [N x 1]

% Train test split
% ...

cv = cvpartition(y, 'HoldOut', 0.3);

% Train and test split
X_train = X(training(cv), :);
y_train = y(training(cv), :);
X_test = X(test(cv), :);
y_test = y(test(cv), :);


% Train model on X1 and X2
bestModel = NaN;
bestElbo = -inf;

for i = 1:10
    K = 100;
    gfaModel =  SGFA({X1', X2'}, K, 15000);
    [elboVals, convIt] = gfaModel.fit(10);
    if elboVals(end) > bestElbo
        bestModel = gfaModel;
    end
end



return;

%% Visualize true and recovered latent factors

% Recovered latent factors
expZ = gfaModel.Z.E;

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
totalD = sum(gfaModel.D); % Total number of dimensions

% trueW = zeros(totalD, data.trueK); % True K
estW = zeros(totalD, gfaModel.K.Val);
estAlpha = zeros(gfaModel.K.Val, gfaModel.M);

d = 0;
for m = 1:gfaModel.M
    Dm = gfaModel.views(m).D;
    estW(d + 1 : d + Dm, :) = gfaModel.views(m).W.E;
    d = d + Dm;

    estAlpha(:, m) = gfaModel.views(m).alpha.E;
end

imshow(estW);

return; 
% Visualization.plotLoadings(trueW, gfaModel.D, 'True W');
Visualization.plotLoadings(estW, gfaModel.D, 'Estimated W', '');

figure;
ax1 = subplot(1, 2, 1);
ax2 = subplot(1, 2, 2);

Visualization.hintonDiagram(-data.alpha, ax1, 'True -1 * alpha');
Visualization.hintonDiagram(-estAlpha, ax2, 'Estimated -1 * alpha');