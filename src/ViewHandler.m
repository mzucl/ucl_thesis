% Format of the view: [D x N], observations are columns of the matrix
classdef ViewHandler
    properties
        data
    end
    
    properties (Dependent)
        N       % Number of samples
        D       % Dimensionality
    end

    methods(Access = private)
        % Helper method that throws an error if index (of the observation) is not valid
        function validateIndex(obj, idx)
            if idx < 1 || idx > obj.N 
                error(['Error in ' class(obj) ': Index out of range.']); 
            end
        end
    end

    methods
        function obj = ViewHandler(data)
            % TODO (high): We assume here that the observations (N of them) are in
            % the columns of 'data', thus data is DxN matrix. Implement something to get rid of this
            % assumption. I could pass in number of samples as well and
            % check which dimension of 'data' corresponds to it, but this
            % will make problems when we import views from external files (maybe).
            if nargin < 1
                error(['Error in class ' class(obj) ': Too few arguments passed.']);
            else
                obj.data = data;
            end
        end
        
        function observation = getObservation(obj, idx, tr)
            if nargin < 2
                error(['Error in class ' class(obj) ': Too few arguments passed.']);
            end
            % Default value for 'tr' is false
            if nargin < 3
                tr = false;
            end

            validateIndex(obj, idx);

            observation = obj.data(:, idx);
            observation = Utility.ternary(tr, observation', observation);
        end

        function normSq = getObservationNormSq(obj, idx)
            observation = obj.getObservation(idx);
            normSq = norm(observation) ^ 2;
        end
      
        function el = getObservationEntry(obj, idx, d)
            % 'd' is the dimension we are interested in
            observation = obj.getObservation(idx);
            if d > 0 && d <= length(observation)
                el = observation(d);
            else
                error(['Error in class ' class(obj) ': Index out of bounds.']);
            end
        end


        %% Getters
        function value = get.D(obj)
            value = size(obj.data, 1);
        end

        function value = get.N(obj)
            value = size(obj.data, 2);
        end
    end
end
