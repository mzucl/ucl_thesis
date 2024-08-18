[X, D] = Datasets.generateBPCA();

%% PPCA
W_PPCA = PPCA(X, D - 1); % PPCA expects X in NxD format


%% BPCA
K = D - 1; % dim - 1, because BPCA can infer effective number of components
numIter = 500;

% BayesianPCA constructor expects data in DxN format
obj = BayesianPCA(X', K, numIter); 
[elboVals, convIt, resArr] = obj.fit();

figure;

ax1 = subplot(1, 2, 1);
Visualization.hintonDiagram(W_PPCA, ax1, 'PPCA W matrix');

ax2 = subplot(1, 2, 2);
Visualization.hintonDiagram(obj.W.EC, ax2, 'BPCA W matrix');