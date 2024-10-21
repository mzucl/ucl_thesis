% Folder path where .mat files are stored
folderPath = 'datasets/ehssan/demographic';

% Get list of all .mat files in the folder
matFiles = dir(fullfile(folderPath, '*.mat'));

% Create a new CSV file to store the output
outputCsvFile = 'demographic.csv';

% Open file for writing
fid = fopen(outputCsvFile, 'w');

% Loop through each .mat file
for i = 1:length(matFiles)
    disp(i);
    % Get the file name
    fileName = matFiles(i).name;
    
    % Load the .mat file
    fileData = load(fullfile(folderPath, fileName));
    
    % Assuming there is only one variable in each .mat file
    % Get the field names of the loaded data
    vars = fieldnames(fileData);
    
    % Assuming the first variable is the matrix to flatten
    matrixData = fileData.(vars{1});

    matrixData = matrixData(1:2);
    
    % Flatten the matrix to a row vector
    % flattenedMatrix = matrixData(:)';  % Transpose to ensure it's a row
    
    fileNameParts = split(fileName, '.');

    % Combine the file name with the flattened matrix
    rowData = [fileNameParts(1), num2cell(matrixData)];
    
    % Convert rowData to string format for writing to CSV
    rowString = strjoin(string(rowData), ',');
    
    % Write the row to the CSV file
    fprintf(fid, '%s\n', rowString);
end

% Close the CSV file
fclose(fid);

disp('CSV file has been created successfully.');
