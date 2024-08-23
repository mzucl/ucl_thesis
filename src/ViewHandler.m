classdef ViewHandler < handle
    % TODO: I just added < handle -> check if that makes any issues for GFA
    % model!!!

    % THIS IS OVERKILL FOR THIS CLASS!
    properties
        X

        N           % Number of samples
        D           % Dimensionality/number of features
        Tr_XtX      % Tr(X^TX)
        XXt         % XX^T
    end
    
    % properties (Dependent)
    %     N           % Number of samples
    %     D           % Dimensionality/number of features
    %     Tr_XtX      % Tr(X^TX)
    %     XXt         % XX^T
    % end

    % Private properties for caching
    properties (Access = private)
        Tr_XtX_Cached
        XXt_Cached
    end

    events
        XChanged
    end

    properties (Access = private)
        controller
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
    
        function setFlags(obj)
            obj.controller.setDirtyFlag('Tr_XtX');
            obj.controller.setDirtyFlag('XXt');
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

            obj.D = size(obj.X, 1);
    
      
            obj.N = size(obj.X, 2);
        
            obj.Tr_XtX = dot(obj.X(:), obj.X(:));
            obj.XXt = obj.X * obj.X';
    
        


            % obj.controller = Controller();
            % addlistener(obj, 'XChanged', @(src, evt) obj.setFlags());
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

        % Returns the row at index idx; if 'tr' is true then the row is
        % returned as a column vector
        function observation = getRow(obj, idx, tr)
            if nargin < 2
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
            end

            % Default value for 'tr' is false
            if nargin < 3
                tr = false;
            end

            validateRowIndex(obj, idx);

            observation = obj.X(idx, :);
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


        %% Setters
        % function set.X(obj, value)
        %     obj.X = value;
        %     notify(obj, 'XChanged');
        % end


        %% Getters
        % function value = get.D(obj)
        %     value = size(obj.X, 1);
        % end
        % 
        % function value = get.N(obj)
        %     value = size(obj.X, 2);
        % end
        % 
        % function value = get.Tr_XtX(obj)
        %     if obj.controller.isDirty('Tr_XtX')
        %         obj.Tr_XtX_Cached = dot(obj.X(:), obj.X(:));
        %         obj.controller.clearDirtyFlag('Tr_XtX');
        %     end
        %     value = obj.Tr_XtX_Cached;
        % end
        % 
        function value = get.XXt(obj)
            if obj.controller.isDirty('XXt')
                obj.XXt_Cached = obj.X * obj.X';
                obj.controller.clearDirtyFlag('XXt');
            end
            value = obj.XXt_Cached;
        end
    end
end
