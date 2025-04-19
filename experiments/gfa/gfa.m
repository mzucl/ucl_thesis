% Clear the workspace
close all; clearvars; clc;

% Logging
logFileName = 'logs/gfa_2G.txt';
if ~exist('logs', 'dir')
    mkdir('logs');
end

diary(logFileName); % start logging

% Model settings
settings = ModelSettings.getInstance();
% settings.VALIDATE = false;
% settings.DEBUG = false;



%% Generate data and train the model
data = Datasets.generateSyntheticGFAData(2);

X1 = data.X_train{1}; % [D1 x N]
X2 = data.X_train{2}; % [D2 x N]

K = 10;

stabilityRun = 2;
modelSelectionIter = 2;
convItAvg = 0;

tic;

for s = 1:stabilityRun 
    maxElbo = -Inf;
    bestW = NaN;
    convIt = NaN;

    for i = 1:modelSelectionIter
        gfaModel =  GFA({X1, X2}, K);
        [elboVals, it] = gfaModel.fit(10);
    
        if elboVals(end) > maxElbo
            maxElbo = elboVals(end);
            bestModel = gfaModel;
            convIt = it;
        end
    end

    convItAvg = convItAvg + convIt;
    disp(['The best model converged in ', num2str(convIt), 'iterations.\n']);
end

elapsedTime = toc;
fprintf('\n\n\nElapsed time: %.4f [s]\n', elapsedTime);
fprintf('Average number of iterations: %.4f\n', convItAvg / stabilityRun);

diary off; 

%% Visualize true and inferred latent factors
Visualization.plotLatentFactors(data.Z, 'True latent factors', '', mfilename);
Visualization.plotLatentFactors(bestModel.Z.E, 'Inferred latent factors', '', mfilename);


%% Visualize loadings and alpha
Visualization.plotFactorLoadingsAndAlpha(data.W, data.D, data.alpha, 'bottom', '', 2.5, 'True $\mathbf{W}^\top$ and \boldmath{$\alpha$}$^\top$', '', mfilename);
Visualization.plotFactorLoadingsAndAlpha(bestModel.W, bestModel.D, bestModel.alpha, 'bottom', '', 2.5, 'Inferred $\mathbf{W}^\top$ and \boldmath{$\alpha$}$^\top$', '', mfilename);