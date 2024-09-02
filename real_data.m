% Import CSV file
data = readtable('control_eeg_questionnaire.csv');

colsToDrop = {'ppid', 'isControl', 'gender', 'age'};

data(:, colsToDrop) = [];

dataArray = table2array(data);

colNumEEG = 155;

X1 = dataArray(:, 1:colNumEEG);
X2 = dataArray(:, colNumEEG+1:end);
