% Clear the workspace
close all; clear all; clc;

% Logging
logFileName = 'logs/sgfa_mnist.txt';
if ~exist('logs', 'dir')
    mkdir('logs');
end

% Figures folder
% figsSubfolder = 'sgfa_mnist';

% diary(logFileName); % start logging

% Model settings
settings = ModelSettings.getInstance();
settings.VALIDATE = false;
settings.DEBUG = false;

folderName = 'mnist18'; % mnist38

%% Import data and train the model
[X, y] = Datasets.getMNISTData(folderName, false, 500);
shuffle_indices = randperm(size(X, 1));
            
% Shuffle rows of both matrices using the same permutation
X = X(shuffle_indices, :);
y = y(shuffle_indices, :);


% Number of folds
k = 5;


stabilityRuns = 10;
K = 400;

%% Storing results
Ks = cell(k, 1);
specificitys = cell(k, 1);
balanced_accuracys = cell(k, 1);
precisions = cell(k, 1);
recalls = cell(k, 1);
f1_scores = cell(k, 1);
npvs = cell(k, 1);
iters = cell(k, 1);
elbos = cell(k, 1);



% Create the k-fold partition
cv = cvpartition(y, 'KFold', k);


%% 
for i = 1:k
    disp(['cv', num2str(i)]);
    % [X_tr, y_tr, X_te, y_te] = Datasets.trainTestSplit(X, y, true);


    % Get the training and test indices for the i-th fold
    trainIdx = training(cv, i);
    testIdx = test(cv, i);
    
    % Split the data into training and test sets
    X_tr = X(trainIdx, :);
    y_tr = y(trainIdx);
    X_te = X(testIdx, :);
    y_te = y(testIdx);


    bestModel = NaN;
    bestElbo = -inf;
    bestIter = 0; % conv. iter of the best model

    for s = 1:stabilityRuns
        model = SGFA({X_tr', y_tr'}, K, 10000, 1e-5); % Views are expected in DxN
        [elboVals, iter] = model.fit(10);
        if elboVals(end) > bestElbo
            bestModel = model;
            bestElbo = elboVals(end);
            bestIter = iter;
        end
    end

    %% Make predictions with the 'bestModel'
    Z = bestModel.Z.E;
    K_eff = size(Z, 1);
    W1 = bestModel.views(1).W.E;
    W2 = bestModel.views(2).W.E;
    mu1 = bestModel.views(1).mu.E;
    mu2 = bestModel.views(2).mu.E;
    T1 = bestModel.views(1).tau.E * eye(bestModel.D(1));
    
    sigma_Z = Utility.matrixInverse(eye(K_eff) + W1' * T1 * W1);

   
    % Extract the best threshold
    train_best_threshold = 0;
    
    
    % Test data

    MU_Z = sigma_Z * (W1' * T1 * (X_te' - mu1));
    pred_te = W2 * MU_Z + mu2;
    
    pred_te = pred_te >= train_best_threshold;



    %% Results
    % Compute the confusion matrix
    y_te = y_te >= train_best_threshold;
    confMat = confusionmat(y_te, pred_te');
    
    % Extract confusion matrix elements
    TN = confMat(1,1);  % True Negatives
    FP = confMat(1,2);  % False Positives
    FN = confMat(2,1);  % False Negatives
    TP = confMat(2,2);  % True Positives
    
    % 1. Balanced Accuracy
    sensitivity = TP / (TP + FN);     % True Positive Rate (Recall)
    specificity = TN / (TN + FP);     % True Negative Rate
    balanced_accuracy = (sensitivity + specificity) / 2;
    
    % 2. Precision
    precision = TP / (TP + FP);
    
    % 3. Recall
    recall = sensitivity;  % Already computed as sensitivity (TP / (TP + FN))
    
    % 4. F1 Score
    f1_score = 2 * (precision * recall) / (precision + recall);
    
    % 5. Negative Predictive Value (NPV)
    npv = TN / (TN + FN);
    
    % Display results
    disp(['Balanced Accuracy: ', num2str(balanced_accuracy)]);
    disp(['Precision: ', num2str(precision)]);
    disp(['Recall: ', num2str(recall)]);
    disp(['F1 Score: ', num2str(f1_score)]);
    disp(['Negative Predictive Value: ', num2str(npv)]);


    % Store results
    Ks{i} = K_eff;
    specificitys{i} = specificity;
    balanced_accuracys{i} = balanced_accuracy;
    precisions{i} = precision;
    recalls{i} = recall;
    f1_scores{i} = f1_score;
    npvs{i} = npv;
    iters{i} = bestIter;
    elbos{i} = bestElbo;
end

%%
clc
numericArray = cell2mat(iters);

% Compute the mean and standard deviation
meanValue = mean(numericArray(:));  % Compute the mean of all elements
stdValue = std(numericArray(:));    % Compute the standard deviation of all elements

% Display results
disp(['Mean: ', num2str(meanValue)]);
disp(['Standard Deviation: ', num2str(stdValue)]);