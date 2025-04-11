%% Clear the workspace
close all; clear all; clc;



%% Generate dataset
[X, D] = Datasets.generateBPCA();



%% Train model
settings = ModelSettings.getInstance();
settings.VALIDATE = false;
settings.DEBUG = false;

% PPCA
[W_PPCA, sigmaSq] = PPCA(X, D - 1); % PPCA expects X in [N x D] format

profile on;

% BPCA
obj = BPCA(X');
[elboVals, it] = obj.fit();

profile off;
profile viewer;
    


%% Validation
if ~Utility.isMonotonicIncreasing(elboVals)
    fprintf(2, 'ELBO decreased at some iteration!!!');
end



%% Visualize
Visualization.plotHintonDiagrams({W_PPCA, obj.W.E}, {'PPCA', 'BPCA'}, '', 'W_ppca_bpca_', mfilename);