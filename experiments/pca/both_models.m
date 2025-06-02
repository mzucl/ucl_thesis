% Clear the workspace
close all; clearvars; clc;

rc = RunConfig.getInstance();
rc.inputValidation = true;
rc.enableLogging = true;

% Logging
logFileName = 'logs/pca_no_init.txt';
if ~exist('logs', 'dir')
    mkdir('logs');
end

diary(logFileName); % start logging


% Generate dataset
[X, D] = Datasets.generateSyntheticBPCAData();

% Save X to CSV for use with datastore
T = array2table(X);
writetable(T, 'Xdata.csv');

T = readtable('Xdata.csv');
numChunks = 3;
chunkSize = height(T) / numChunks;

if isfolder("dataChunks")
    rmdir("dataChunks", 's');  % 's' stands for recursive delete (subfolders & files)
end


mkdir('dataChunks');
for i = 1:numChunks
    idx = (1:chunkSize) + (i-1)*chunkSize;
    chunk = T(idx, :);
    writetable(chunk, sprintf('dataChunks/Xdata_part%d.csv', i));
end

% PPCA
[W_PPCA, sigmaSq] = PPCA(X, D - 1); % PPCA expects X in [N x D] format

% obj = BPCA(X');
obj = BPCA_mr('fgg');
        
[elboVals, it] = obj.fit();

%% Visualize
Visualization.plotHintonDiagrams({W_PPCA, obj.W.E}, {'PPCA', 'BPCA'}, '', ['W_ppca_bpca_', '12345'], 'lala');