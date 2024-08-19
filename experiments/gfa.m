data = generateTwoViews();

X1 = data.X_tr{1}; % [N x D1];
X2 = data.X_tr{2}; % [N x D2]

% Scale datasets
% X1 = Datasets.standardScaler(X1);
% X2 = Datasets.standardScaler(X2);

maxIter = 500;
K = 10;

gfaModel =  GFA({X1', X2'}, K, maxIter);
[elboVals, convIt, resArr] = gfaModel.fit();
Visualization.plotStructVariables(resArr);

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
expZ = gfaModel.Z.EC;

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