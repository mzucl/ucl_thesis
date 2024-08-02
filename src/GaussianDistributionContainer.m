classdef GaussianDistributionContainer
    properties
        distributions
    end
    
    properties (Dependent)
        Size                % Number of distributions in the container
        Expectation         % Under the assumption of independence of each component; 
                            % This is a cell array where each entry is an
                            % expectation of one component;
        ExpectationXXt      % Similar to above
        ExpectationC        % 'C' stands for container
        ExpectationCt 
        ExpectationCtC
        ExpectationXSqNorm  % This is a cell array where each entry is E[|component|^2]
    end

    methods
        function obj = GaussianDistributionContainer(dim, numDistributions, prec)
            % Constructor with 2 scalar parameters: 'dim',
            % 'numDistributions'
            %   All distributions are standard normal multivariate
            %   Gaussians with the same dimensionality
            
            % TODO: Add more initialization options, similar to
            % GammaDistribution;

            % [NOTE] All distibutions in the container have the same 'dim'.

            % 'prec' has to be an array of values - one per distribution ->
            % assert that this is true

            if nargin < 2
                error(['Error in class ' class(obj) ': Too few arguments passed.']);
            elseif nargin == 2
                obj.distributions = repmat(GaussianDistribution(), numDistributions, 1); % Preallocate
                for i = 1:numDistributions
                    obj.distributions(i) = GaussianDistribution(0, 1, dim);
                end
            elseif nargin == 3
                obj.distributions = repmat(GaussianDistribution(), numDistributions, 1); % Preallocate
                for i = 1:numDistributions
                    obj.distributions(i) = GaussianDistribution(0, 1./prec ...
                        (i), dim);
                end
            end
        end
       


        %% Methods
        function dist = getDistribution(obj, idx)
            % Returns the distribution at index 'idx'
            if idx > obj.Size
                error(['Error in ' class(obj) ': Index out of range.']); 
                % TODO: Add a separate function for index validation
            end
            dist = obj.distributions(idx);
        end

        function obj = updateDistribution(obj, idx, dist)
            if idx > obj.Size
                error(['Error in ' class(obj) ': Index out of range.']);
            end
            % Updates the distribution at index 'idx'
            obj.distributions(idx) = dist;
        end

        function obj = updateDistributionParams(obj, idx, mu, cov)
            if nargin < 3
                error(['Error in class ' class(obj) ': Too few arguments passed.']);
            end
            if idx > obj.Size
                error(['Error in ' class(obj) ': Index out of range.']);
            end

            obj.distributions(idx).updateParameters(mu, cov);
        end

        function obj = updateDistributionMu(obj, idx, mu)
            if nargin < 2
                error(['Error in class ' class(obj) ': Too few arguments passed.']);
            end
            if idx > obj.Size
                error(['Error in ' class(obj) ': Index out of range.']);
            end

            obj.distributions(idx).updateMu(mu);
        end

        % [NOTE] The type of these update methods depends on the current
        % need (e.g. in the update equations all covariances are set to the
        % same value). In future more general methods can be added (maybe
        % needed in the update equations of more complex models like GFA).

        % [NOTE] We shouldn't update parameters ('mu' and 'cov')
        % separately, so this should be refactored (as well as the code in
        % the GaussianDistribution class.
        function obj = updateAllDistributionsCovariance(obj, cov)
            % Update covariance of all distributions
            if nargin < 1
                error(['Error in class ' class(obj) ': Too few arguments passed.']);
            end
            
            % TODO: Check if 'cov' is a valid covariance matrix
            for i=1:obj.Size
                obj.distributions(i).updateCovariance(cov);
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

        function value = get.ExpectationC(obj)
            value = cell2mat(obj.Expectation);
        end

        function value = get.ExpectationCt(obj)
            value = obj.ExpectationC';
        end

        function value = get.ExpectationXXt(obj)
            value = cell(1, obj.Size);
            for i = 1:obj.Size
                value{i} = obj.distributions(i).ExpectationXXt;
            end
        end
        
        function value = get.ExpectationCtC(obj)
            value = zeros(obj.Size, obj.Size);
            for i = 1:obj.Size
                for j = 1:obj.Size
                    value(i, j) = obj.distributions(i).ExpectationXt * ...
                        obj.distributions(j).Expectation;
                end
            end
        end

        function value = get.ExpectationXSqNorm(obj)
            value = cell(1, obj.Size);
            for i = 1:obj.Size
                value{i} = obj.distributions(i).ExpectationXtX;
            end
        end
    end
end