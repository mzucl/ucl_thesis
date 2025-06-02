% Clear the workspace
close all; clearvars; clc;

rc = RunConfig.getInstance();
rc.inputValidation = false;
rc.enableLogging = false;


% profile on;

obj = BPCA_mr('fgg');
[elboVals, it] = obj.fit();

% profile off;
% profile viewer;
% 

% PPCA
% Generate dataset
[X, D] = Datasets.generateSyntheticBPCAData();
[W_PPCA, sigmaSq] = PPCA(X, D - 1); % PPCA expects X in [N x D] format


%% Validation
if ~Utility.isMonotonicIncreasing(elboVals)
    fprintf(2, 'ELBO decreased at some iteration!!!');
end


%% Visualize
Visualization.plotHintonDiagrams({W_PPCA, obj.W.E}, {'PPCA', 'BPCA'}, '', 'W_ppca_bpca_', mfilename);