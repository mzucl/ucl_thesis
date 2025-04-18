% Clear the workspace
close all; clearvars; clc;

%% Generate data and run the experiment
data = Datasets.generateSyntheticGFAData(2);

bestOverallModel = Experiment('sgfa', 10, {data.X_train{1}, data.X_train{2}}, mfilename).run();

%% [NOTE] The visualization functionality is deliberately kept outside of 
% the `Experiment` class to allow for greater flexibility and control over 
% aspects such as figure title, name, and other presentation settings.
%% Visualize true and inferred latent factors
Visualization.plotLatentFactors(data.Z, 'True latent factors', '', mfilename);
Visualization.plotLatentFactors(bestOverallModel.Z.E, 'Inferred latent factors', '', mfilename);

%% Visualize loadings and alpha
Visualization.plotFactorLoadingsAndAlpha(data.W, data.D, data.alpha, 'bottom', '', 2.5, 'True $\mathbf{W}^\top$ and \boldmath{$\alpha$}$^\top$');
Visualization.plotFactorLoadingsAndAlpha(bestOverallModel.W, bestOverallModel.D, bestOverallModel.alpha, 'bottom', '', 2.5, 'Inferred $\mathbf{W}^\top$ and \boldmath{$\alpha$}$^\top$');