% Author: Mediha Zukic
% Contact: mediha.zukic.23@alumni.ucl.ac.uk
% Date: 2025-05-13

% NaN == NaN -> 0
% NaN ~= NaN -> 1

% Save value to workspace
% assignin('base', 'objT', obj);

% Set the background color of the plot
set(gca, 'Color', [0.8 0.8 0.8]); % Light gray background color

% Turn off the axis
axis off;

% Optionally, turn off the figure background color
set(gcf, 'Color', [1 1 1]); % White background color for the figure


% Use testCase.verifyEqual with Tolerances for Float Types
% testCase.verifyEqual(result1, result2, 'AbsTol', 1e-10)