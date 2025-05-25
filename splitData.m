[X, D] = Datasets.generateSyntheticBPCAData();

% Save X to CSV for use with datastore
T = array2table(X);
writetable(T, 'Xdata.csv');

T = readtable('Xdata.csv');
numChunks = 3;
chunkSize = height(T) / numChunks;

mkdir('dataChunks');
for i = 1:numChunks
    idx = (1:chunkSize) + (i-1)*chunkSize;
    chunk = T(idx, :);
    writetable(chunk, sprintf('dataChunks/Xdata_part%d.csv', i));
end