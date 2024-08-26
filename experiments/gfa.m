%% Generate data and train the model
data = generateTwoViews();

X1 = data.X_tr{1}; % [N x D1];
X2 = data.X_tr{2}; % [N x D2]

% Scale datasets
% X1 = Datasets.standardScaler(X1);
% X2 = Datasets.standardScaler(X2);

maxIter = 1000;
K = 10;

profile on;

gfaModel =  GFA({X1', X2'}, K, maxIter);
[elboVals, convIt, resArr] = gfaModel.fit();

profile off;
profile viewer;

%% Visualize model training
% Visualization.plotStructVariables(resArr);



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

trueW = zeros(totalD, data.trueK); % True K
estW = zeros(totalD, gfaModel.K.Val);
estAlpha = zeros(gfaModel.K.Val, gfaModel.M);

d = 0;
for m = 1:gfaModel.M
    Dm = gfaModel.views(m).D;
    trueW(d + 1 : d + Dm, :) = data.W{m};
    estW(d + 1 : d + Dm, :) = gfaModel.views(m).W.E;
    d = d + Dm;

    estAlpha(:, m) = gfaModel.views(m).alpha.E;
end

Visualization.plotLoadings(trueW, gfaModel.D, 'True W');
Visualization.plotLoadings(estW, gfaModel.D, 'Estimated W');

figure;
ax1 = subplot(1, 2, 1);
ax2 = subplot(1, 2, 2);

Visualization.hintonDiagram(-data.alpha, ax1, 'True -1 * alpha');
Visualization.hintonDiagram(-estAlpha, ax2, 'Estimated -1 * alpha');