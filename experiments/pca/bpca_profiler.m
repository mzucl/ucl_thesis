% Clear the workspace
close all; clearvars; clc;

rc = RunConfig.getInstance();
rc.inputValidation = false;
rc.enableLogging = false;

%% Generate dataset
[X, D] = Datasets.generateSyntheticBPCAData(50000, 100);

%% PPCA
[W_PPCA, sigmaSq] = PPCA(X);

%% BPCA
profile on;

obj = BPCA(X);
[elboVals, it] = obj.fit();

profile off;
profile viewer;
    
%% Validate ELBO values
if ~MatrixValidation.isMonotonicIncreasing(elboVals)
    fprintf(2, 'ELBO decreased at some iteration!!!');
end

%% Visualize loadings
Visualization.plotHintonDiagrams({W_PPCA, obj.W.E}, {'PPCA', 'BPCA'}, '', 'W_ppca_bpca_', mfilename);