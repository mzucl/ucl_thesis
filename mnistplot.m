% Author: Mediha Zukic
% Contact: mediha.zukic.23@alumni.ucl.ac.uk
% Date: 2025-05-13

folderNames = ["mnist18", "mnist38"];

for i = 1:length(folderNames) 
    folderName = char(folderNames(i));

    X = readmatrix(['datasets/', folderName, '/pixels.tsv'], 'FileType', 'text'); % [N x D];
    img1 = reshape(X(1, :), [28, 28])';
    img2 = reshape(X(end, :), [28, 28])';
    
    hfig = figure;
    
    subplot(1, 2, 1);
    imshow(img1);
    
    subplot(1, 2, 2);
    imshow(img2);
    
    Visualization.formatFigure(hfig);
    Visualization.exportFigure(hfig, folderName, 'mnist');
end



