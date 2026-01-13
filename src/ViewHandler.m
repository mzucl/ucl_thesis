classdef ViewHandler < handle
    properties
        X

        %% Constant Dependent properties
        % Defined like this because they never change after initialization
        N           % Number of samples
        D           % Dimensionality/number of features
        Tr_XtX      % Tr(X^TX)
        XXt         % XX^T
    end



    
    methods(Access = private)
        % Helper method that throws an error if index (of the observation) is not valid
        function validateIndex(obj, idx)
            if idx < 1 || idx > obj.N 
                error(['Error in ' class(obj) ': Index out of range.']); 
            end
        end

        function validateRowIndex(obj, idx)
            if idx < 1 || idx > obj.D 
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

            obj.X = LogicUtils.ternary(featuresInCols, data', data);

            % Set dependent properties
            obj.D = size(obj.X, 1);
            obj.N = size(obj.X, 2);
            obj.Tr_XtX = dot(obj.X(:), obj.X(:));
            obj.XXt = obj.X * obj.X';
        end
        


        

        %% Retreive methods
        function observation = getObservation(obj, idx, tr)
            if RunConfig.getInstance().inputValidation && nargin < 2
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
            end

            % Default value for 'tr' is false
            if nargin < 3
                tr = false;
            end

            validateIndex(obj, idx);

            observation = obj.X(:, idx);
            observation = LogicUtils.ternary(tr, observation', observation);
        end

        % Returns the row at index idx; if 'tr' is true then the row is
        % returned as a column vector
        function observation = getRow(obj, idx, tr)
            if RunConfig.getInstance().inputValidation && nargin < 2
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
            end

            % Default value for 'tr' is false
            if nargin < 3
                tr = false;
            end

            validateRowIndex(obj, idx);

            observation = obj.X(idx, :);
            observation = LogicUtils.ternary(tr, observation', observation);
        end

        function normSq = getObservationNormSq(obj, idx)
            observation = obj.getObservation(idx);
            normSq = observation' * observation;
        end
      
        function el = getObservationEntry(obj, idx, d)
            % 'd' is the dimension we are interested in
            observation = obj.getObservation(idx);
            if RunConfig.getInstance().inputValidation
                if d < 1 || d > length(observation)
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Index out of bounds.']);
                end
            end
            el = observation(d);
        end
    end
end
