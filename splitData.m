T = readtable('Xdata.csv');
numChunks = 3;
chunkSize = width(T) / numChunks;

for i = 1:numChunks
    idx = (1:chunkSize) + (i-1)*chunkSize;
    chunk = T(:, idx);
    writetable(chunk, sprintf('dataChunks/Xdata_part%d.csv', i));
end