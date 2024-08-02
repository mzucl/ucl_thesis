classdef GammaDistributionContainer < handle
    properties
        distributions
    end
    
    properties (Dependent)
        Size            % Number of distributions in the container
        Expectation     % Under the assumption of independence of each component
                        % This is a cell array where each entry is an
                        % expectation of one component; 
                        % obj.Expectation{i}
        Value           
    end

    methods(Access = private)
        % Helper method that throws an error if index is not valid
        function validateIndex(obj, idx)
            if idx < 1 || idx > obj.Size 
                error(['Error in ' class(obj) ': Index out of range.']); 
            end
        end
    end

    methods
        function obj = GammaDistributionContainer(dists, b, numDistributions)
            % Constructor with 1 parameter: array of distributions
            if nargin == 1 && Utility.areAllInstancesOf(dists, 'GammaDistribution')
                obj.distributions = dists;
            
            % Constructor with 2 parameters: different values for parameters of each
            % distribution; 'dists' is now array of 'a' parameters
            % a: []
            % b: []
            elseif nargin == 2 && Utility.isArray(dists) && Utility.isArray(b)
                if length(dists) ~= length(b)
                    error(['Error in ' class(obj) ': Arrays of values' ...
                        ' for a and b should be of the same size.']);
                end
                a = dists;
                numInstances = length(a);
                obj.distributions = repmat(GammaDistribution(), numInstances, 1); % Preallocate
 
                for i = 1:numInstances
                    obj.distributions(i) = GammaDistribution(a(i), b(i));
                end
            
            % Constructor with 2 parameters: different values for 'b' parameters of each
            % distribution; 'dists' is now the value of 'a' parameter that
            % is the same for all distibutions
            % a: scalar
            % b: []
            elseif nargin == 2 && isscalar(dists) && Utility.isArray(b)
                a = dists;
                numInstances = length(b);
                obj.distributions = repmat(GammaDistribution(), numInstances, 1); % Preallocate
 
                for i = 1:numInstances
                    obj.distributions(i) = GammaDistribution(a, b(i));
                end

            % Constructor with 2 parameters: different values for 'a' parameter for each
            % distribution; 'dists' is now array of 'a' parameters and 'b'
            % is the same for all distributions
            % a: []
            % b: scalar
            elseif nargin == 2 && Utility.isArray(dists) && isscalar(b)
                a = dists;
                numInstances = length(a);
                obj.distributions = repmat(GammaDistribution(), numInstances, 1); % Preallocate
 
                for i = 1:numInstances
                    obj.distributions(i) = GammaDistribution(a(i), b);
                end
            
            % Constructor with 3 parameters: all distributions have the
            % same parameter values for 'a' and 'b'
            % a: scalar
            % b: scalar 
            elseif nargin == 3 && isnumeric(dists) && isnumeric(b) && isnumeric(numDistributions)
                a = dists;

                obj.distributions = repmat(GammaDistribution(), numDistributions, 1); % Preallocate
 
                for i = 1:numDistributions
                    obj.distributions(i) = GammaDistribution(a, b);
                end

            else
                error(['Error in ' class(obj) ': Invalid input arguments.']);
            end
        end

        %% Single distibution methods
        function dist = getDistribution(obj, idx)
            % Returns the distribution at index 'idx'
            obj.validateIndex(idx);

            dist = obj.distributions(idx);
        end

        function obj = updateDistribution(obj, idx, dist)
            % Updates the distribution at index 'idx'
            obj.validateIndex(idx);

            obj.distributions(idx) = dist;
        end

        function obj = updateDistributionParams(obj, idx, a, b, inc)
            if nargin < 4
                error(['Error in class ' class(obj) ': Too few arguments passed.']);
            end
            obj.validateIndex(idx);

            % Default for 'inc' is true
            if nargin == 4
                inc = true;
            end

            obj.distributions(idx).updateParameters(a, b, inc);
        end

        %% All distibutions methods
        function obj = updateAllDistributionsParams(obj, a, b, inc)
            if nargin < 3
                error(['Error in class ' class(obj) ': Too few arguments passed.']);
            end
            % Default for 'inc' is true
            if nargin == 3
                inc = true;
            end

            % Check 'a'
            if isscalar(a)
                a = a * ones(1, obj.Size);
            elseif Utility.isArray(a)
                if length(a) ~= obj.Size
                    error(['Error in ' class(obj) ': Dimensions do not match.']);
                end
            end

            % Check 'b'
            if isscalar(b)
                b = b * ones(1, obj.Size);
            elseif Utility.isArray(b)
                if length(b) ~= obj.Size
                    error(['Error in ' class(obj) ': Dimensions do not match.']);
                end
            end

            for i=1:obj.Size
                obj.distributions(i).updateParameters(a(i), b(i), inc);
            end
        end



        %% Getters
        function value = get.Size(obj)
            value = length(obj.distributions);
        end

        function value = get.Expectation(obj)
            value = cell(1, obj.Size);
            for i = 1:obj.Size
                value{i} = obj.distributions(i).Expectation;
            end
        end

        function value = get.Value(obj)
            value = zeros(obj.Size, 1);
            for i = 1:obj.Size
                value(i) = obj.distributions(i).Value;
            end
        end
    end
end