[X, D] = Datasets.generateBPCA();

%% PPCA
[W_PPCA, sigmaSq] = PPCA(X, D - 1); % PPCA expects X in [N x D] format


%% BPCA
numIter = 500;

profile on;

% BayesianPCA constructor expects data in [D x N] format
bpcaModel = BayesianPCA(X', numIter); 
[elboVals, convIt, resArr] = bpcaModel.fit();

profile off;
profile viewer;



%% Visualize
figure;

ax1 = subplot(1, 2, 1);
Visualization.hintonDiagram(W_PPCA, ax1, 'PPCA W matrix');

ax2 = subplot(1, 2, 2);
Visualization.hintonDiagram(bpcaModel.W.EC, ax2, 'BPCA W matrix');