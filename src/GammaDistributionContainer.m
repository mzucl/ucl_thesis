%% [NOTE]
% Components inside the container are always independent, thus always
% appear as factor inside the product, e.g. p(Z) is product of p(zn) for
% each zn (zn is a latent variable corresponding to the observation xn).
%%
classdef GammaDistributionContainer < handle
    properties
        distributions
    end
    
    properties (Dependent)
        Size            % Number of distributions in the container
        Expectation     % Under the assumption of independence of each component
                        % This is a cell array where each entry is an
                        % expectation of one component; 
                        % e.g. access ith component expectation using obj.Expectation{i}
        ExpectationC    % Array of expectations of all components
        H               % This is a cell array where each entry is an
                        % entropy of one component; 
        HC              % This is the entropy of the collection, under the assumptions given in the [NOTE] above
        
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
        function obj = GammaDistributionContainer(a, b, priors, numDistributions)
            % Left 'numDistributions' as a last parameter because it can be
            % infered from other parameters in some cases. 'priors' is the
            % same for all distributions for now (or NaN can be passed in), but code can be extended
            % to accomodate for different priors if needed (thus the name in plural).
            % (Utility.areAllInstancesOf(dists, 'GammaDistribution')) - if
            % the 'priors' is an array, all elements have to be instances
            % of 'GammaDistribution' class.

            if nargin < 2
                error(['Error in class ' class(obj) ': Too few arguments passed.']);

            elseif nargin == 2 || nargin == 3
                hasPriors = nargin == 3 && ~(isnumeric(priors) && isnan(priors)); % priors is not 'NaN'
                if nargin == 2 && ~hasPriors % This is necessary because when nargin == 2 then 'priors' is not defined
                    priors = NaN;
                end
                % Different values for parameters of each distribution;
                % a: []
                % b: []
                if Utility.isArray(a) && Utility.isArray(b)
                    if length(a) ~= length(b)
                        error(['Error in ' class(obj) ': Arrays of values' ...
                            ' for a and b should be of the same size.']);
                    end
                    numDistributions = length(a);
                    obj.distributions = repmat(GammaDistribution(), numDistributions, 1); % Preallocate
     
                    for i = 1:numDistributions
                        obj.distributions(i) = Utility.ternary(hasPriors, ...
                            GammaDistribution(a(i), b(i), priors), GammaDistribution(a(i), b(i)));
                    end
                end


                % Different values of 'b' parameter for each distribution;
                % 'a' is the same for all distributions.
                % a: scalar
                % b: []
                if isscalar(a) && Utility.isArray(b)
                    numDistributions = length(b);
                    obj.distributions = repmat(GammaDistribution(), numDistributions, 1); % Preallocate
     
                    for i = 1:numDistributions
                        obj.distributions(i) = Utility.ternary(hasPriors, ...
                            GammaDistribution(a, b(i), priors), GammaDistribution(a, b(i)));
                    end
                end

                % Different values of 'a' parameter for each
                % distribution; 'b' is the same for all distributions.
                % a: []
                % b: scalar
                if Utility.isArray(a) && isscalar(b)
                    numDistributions = length(a);
                    obj.distributions = repmat(GammaDistribution(), numDistributions, 1); % Preallocate
     
                    for i = 1:numDistributions
                        obj.distributions(i) = Utility.ternary(hasPriors, ...
                            GammaDistribution(a(i), b, priors), GammaDistribution(a(i), b));
                    end
                end
                % Container with only 1 distribution
                % a: scalar
                % b: scalar
                if isscalar(a) && isscalar(b)
                    numDistributions = 1;
                    obj.distributions = repmat(GammaDistribution(), numDistributions, 1); % Preallocate

                    for i = 1:numDistributions
                        obj.distributions(i) = Utility.ternary(hasPriors, GammaDistribution(a, b, priors), ...
                            GammaDistribution(a, b));
                    end
                end
            
            % Constructor with 4 parameters: all distributions have the
            % same parameter values for 'a', 'b' and 'prior' parameter. 'a'
            % and 'b' are scalar in this case because the number of
            % distributions is defined, so it could result in conflicting
            % value if I allow them (one or both) to be arrays. This can be
            % added later as an extension if needed.
            % a: scalar
            % b: scalar 
            elseif nargin == 4 && ...
                isnumeric(a) && ...
                isnumeric(b) && ...
                Utility.isNaNOrInstanceOf(priors, 'GammaDistribution') && ...
                isnumeric(numDistributions)

                obj.distributions = repmat(GammaDistribution(), numDistributions, 1); % Preallocate
                
                hasPriors = ~(isnumeric(priors) && isnan(priors)); % 'priors' is set
                for i = 1:numDistributions
                    obj.distributions(i) = Utility.ternary(hasPriors, GammaDistribution(a, b, priors), ...
                        GammaDistribution(a, b));
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

        function value = get.H(obj)
            value = cell(1, obj.Size);
            for i = 1:obj.Size
                value{i} = obj.distributions(i).H;
            end
        end

        function value = get.HC(obj)
            value = 0;
            for i = 1:obj.Size
                value = value + obj.distributions(i).H;
            end
        end

        function value = get.ExpectationC(obj)
            value = cell2mat(obj.Expectation);
        end

        function value = get.Value(obj)
            value = zeros(obj.Size, 1);
            for i = 1:obj.Size
                value(i) = obj.distributions(i).Value;
            end
        end
    end
end