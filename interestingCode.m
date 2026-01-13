% Copyright (c) 2025 Mediha Zukic
% Licensed under the MIT License.


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

    methods (Static)
        function obj = loadobj(s)
            % LOADOBJ Loads a saved Gamma object from a struct
            obj = Gamma(s.a, s.b);
            if isfield(s, 'priorClass') && strcmp(s.priorClass, 'Gamma')
                % Recursively load prior if it exists
                obj.prior = Gamma.loadobj(s.prior);  
            end
        end
    end



    methods
        %% Save object
        function s = saveobj(obj)
            % SAVEOBJ Converts Gamma object to a struct for saving
            s.a = obj.a;
            s.b = obj.b;
            if isa(obj.prior, 'Gamma')
                s.prior = saveobj(obj.prior);  % recursively save prior
                s.priorClass = 'Gamma';
            end
        end


                function loss = BCE(yTrue, yPred)
            % yTrue: vector of true labels (0 or 1)
            % yPred: vector of predicted probabilities (between 0 and 1)
            
            % Ensure predictions are in a valid range [epsilon, 1 - epsilon] to avoid log(0)
            epsilon = 1e-15; 
            yPred = max(min(yPred, 1 - epsilon), epsilon);
        
            loss = -mean(yTrue .* log(yPred) + (1 - yTrue) .* log(1 - yPred));
        end