% % Clear the workspace
% close all; clearvars; clc;
% 
% RunConfig.resetToDefaults();
% 
% %% Generate data and run the experiment
% data = Datasets.generateSyntheticGFAData(2);
% 
% bestOverallModel = Experiment('sgfa', 10, {data.X_train{1}, data.X_train{2}}, mfilename, '', '', '', false, false).run();
% 
% %% [NOTE] The visualization functionality is deliberately kept outside of 
% % the `Experiment` class to allow for greater flexibility and control over 
% % aspects such as figure title, name, and other presentation settings.
% %% Visualize true and inferred latent factors
% Visualization.plotLatentFactors(data.Z, 'True latent factors', '', mfilename);
% Visualization.plotLatentFactors(bestOverallModel.Z.E, 'Inferred latent factors', '', mfilename);
% 
% %% Visualize loadings and alpha
% Visualization.plotFactorLoadingsAndAlpha(data.W, data.D, data.alpha, 'bottom', '', 2.5, 'True $\mathbf{W}^\top$ and \boldmath{$\alpha$}$^\top$');
% Visualization.plotFactorLoadingsAndAlpha(bestOverallModel.W, bestOverallModel.D, bestOverallModel.alpha, 'bottom', '', 2.5, 'Inferred $\mathbf{W}^\top$ and \boldmath{$\alpha$}$^\top$');
% 
% 






% Clear the workspace
close all; clearvars; clc;

rc = RunConfig.getInstance();
rc.inputValidation = false;
rc.enableLogging = false;

% Logging
logFileName = 'logs/pca_init.txt';
if ~exist('logs', 'dir')
    mkdir('logs');
end

diary(logFileName); % start logging

% Generate dataset
[X, D] = Datasets.generateSyntheticBPCAData();

% PPCA
[W_PPCA, sigmaSq] = PPCA(X, D - 1); % PPCA expects X in [N x D] format


%% BPCA - no initialization setup
stabilityRun = 3;
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
    disp(['The best model converged in ', num2str(convIt), ' iterations.']);

    % arrW{s} = bestW;

    %% Visualize
    Visualization.plotHintonDiagrams({W_PPCA, bestW}, {'PPCA', 'BPCA'}, '', ['W_ppca_bpca_', num2str(s)], mfilename);
end

elapsedTime = toc;
fprintf('\n\n\nElapsed time: %.4f [s]\n', elapsedTime);
fprintf('Average number of iterations: %d\n', round(convItAvg / stabilityRun));

diary off; 
