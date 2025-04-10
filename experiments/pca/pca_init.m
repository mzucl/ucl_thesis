% Clear the workspace
close all; clearvars; clc; clear all;

% Logging
logFileName = 'logs/pca_init.txt';
if ~exist('logs', 'dir')
    mkdir('logs');
end

diary(logFileName); % start logging

% Model settings
settings = ModelSettings.getInstance();
% settings.VALIDATE = false;
% settings.DEBUG = false;

% Generate dataset
[X, D] = Datasets.generateBPCA();

% PPCA
[W_PPCA, sigmaSq] = PPCA(X, D - 1); % PPCA expects X in [N x D] format


%% BPCA - no initialization setup
stabilityRun = 20;
modelSelectionIter = 10;
convItAvg = 0;

tic;

% arrW = {}; % Ws of the 'winning' models

for s = 1:stabilityRun 
    maxElbo = -Inf;
    bestW = NaN;
    convIt = NaN;
    for i = 1:modelSelectionIter
        % BPCA constructor expects data in [D x N] format
        obj = BPCA(X', W_PPCA + randn(size(W_PPCA)));
        
        [elboVals, it] = obj.fit();

        % Validation
        if ~Utility.isMonotonicIncreasing(elboVals)
            fprintf(2, 'ELBO decreased at some iteration!!!');
        end

        if elboVals(end) > maxElbo
            maxElbo = elboVals(end);
            bestW = obj.W.E;
            convIt = it;
        end
    end
    convItAvg = convItAvg + convIt;
    disp(['The best model converged in ', num2str(convIt), 'iterations.\n']);

    % arrW{s} = bestW;

    %% Visualize
    Visualization.hintonDiagramPlot({W_PPCA, bestW}, 'figures/pca_init', ['W_ppca_bpca_', num2str(s)], true);
end

elapsedTime = toc;
fprintf('\n\n\nElapsed time: %.4f [s]\n', elapsedTime);
fprintf('Average number of iterations: %.4f\n', convItAvg / stabilityRun);

diary off; 
