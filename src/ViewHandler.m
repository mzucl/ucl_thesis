% Class to handle the dataset, in GFA that it is termed view
% TODO (medium): Extend this to enable passing in the file path instead of
% passing in the data matrix

% [NOTE] Not sure if it is a good idea to merge generation of the view data
% (by e.g. calling 'generateSyntheticData' and 'handling'; probably better
% to separate the two, but something to think about.

classdef ViewHandler
    properties
        data
    end
    
    properties (Dependent)
        N       % Number of samples
        D       % Dimensionality
    end

    methods
        function obj = ViewHandler(data)
            % TODO (medium): We assume here that the observations (N of them) are in
            % the columns of 'data', thus data is DxN matrix. Implement something to get rid of this
            % assumption.
            if nargin < 1
                error(['Error in class ' class(obj) ': Too few arguments passed.']);
            else
                obj.data = data;
            end
        end
        
        function col = getObservation(obj, idx, tr)
            % Default value for 'tr' is False
            if nargin < 3
                tr = 0;
            end
            if idx > 0 && idx <= obj.N
                observation = obj.data(:, idx);
                col = Utility.ternary(tr == 0, observation, observation');
            else
                error(['Error in class ' class(obj) ': Index out of bounds.']);
            end
        end

        function normSq = getObservationNormSq(obj, idx)
            % TODO: Add try-catch blocks
            col = obj.getObservation(idx);
            normSq = norm(col) ^ 2;
        end
      
        function el = getObservationEntry(obj, idx, d)
            % 'd' is the dimension we are interested in
            col = obj.getObservation(idx);
            if d > 0 && d <= length(col)
                el = col(d);
            else
                error(['Error in class ' class(obj) ': Index out of bounds.']);
            end
        end


        %% Getters
        function value = get.N(obj)
            value = size(obj.data, 2);
        end

        function value = get.D(obj)
            value = size(obj.data, 1);
        end
    end
end
