[X, D] = Datasets.generateBPCA();

%% PPCA
[W_PPCA, sigmaSq] = PPCA(X, D - 1); % PPCA expects X in [N x D] format


%% BPCA
numIter = 500;

profile on;

% BPCA constructor expects data in [D x N] format
obj = BPCA(X', numIter); 

[elboVals, it] = obj.fit(5);

profile off;
profile viewer;

%% Visualize
Visualization.hintonDiagramPlot({W_PPCA, obj.W.E}, 'figures/pca', 'W_ppca_bpca');
