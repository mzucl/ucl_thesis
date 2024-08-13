classdef ViewHandler
    properties
        X
    end
    
    properties (Dependent)
        N       % Number of samples
        D       % Dimensionality/number of features
        TrXtX   % Tr(X^TX)
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
        function obj = ViewHandler(data, featuresInCols)
            % 'featuresInCols' (default: true) tells if the features are
            % stored in columns or rows of the 'data' matrix.
            
            % 'X' that is the part of the ViewHandler is in the [D x N]
            % format, that is observations are stored in the columns of 'X'
            if nargin < 1
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
            
            % [NOTE] I choose default to be 'true' because usually tabular
            % data has features in the columns
            elseif nargin < 2
                featuresInCols = true;
            end

            obj.X = Utility.ternary(featuresInCols, data', data);
        end
        
        

        %% Retreive methods
        function observation = getObservation(obj, idx, tr)
            if nargin < 2
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
            end

            % Default value for 'tr' is false
            if nargin < 3
                tr = false;
            end

            validateIndex(obj, idx);

            observation = obj.X(:, idx);
            observation = Utility.ternary(tr, observation', observation);
        end

        function normSq = getObservationNormSq(obj, idx)
            observation = obj.getObservation(idx);
            normSq = observation' * observation;
        end
      
        function el = getObservationEntry(obj, idx, d)
            % 'd' is the dimension we are interested in
            observation = obj.getObservation(idx);
            if d > 0 && d <= length(observation)
                el = observation(d);
            else
                error(['##### ERROR IN THE CLASS ' class(obj) ': Index out of bounds.']);
            end
        end



        %% Getters
        function value = get.D(obj)
            value = size(obj.X, 1);
        end

        function value = get.N(obj)
            value = size(obj.X, 2);
        end

        function value = get.TrXtX(obj)
            value = trace(obj.X' * obj.X);
        end
    end
end
