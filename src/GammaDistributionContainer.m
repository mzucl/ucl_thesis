%% [NOTE]
% Components inside the container are always independent, thus they always
% appear as factors inside the product, e.g. p(Z) is product of p(zn) for
% each zn (zn is a latent variable corresponding to the observation xn).
%%
classdef GammaDistributionContainer < handle
    properties
        distributions
    end
    
    properties (Dependent)
        Size                % Number of distributions in the container
        Expectation         % This is a cell array where each entry is an
                            % expectation of one component (under the assumption of independence); 
                            % e.g. access ith component expectation using obj.Expectation{i}
        ExpectationC        % Array of expectations of all components
        H                   % This is a cell array where each entry is an
                            % entropy of one component; 
        HC                  % This is the entropy of the collection (sum of entries in H)
        ExpectationLnP      % Cell array where each entry is ExpectationLnP
        ExpectationLnPC     % The sum of all entries in ExpectationLnP
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
            % infered from other parameters in some cases. 
            %
            % 'priors' is either single GammaDistribution in which case all
            % distributions have the same prior, or an array of
            % GammaDistributions - one for each of the components. It can
            % also be set to NaN in which case all distributions have NaN
            % priors.
            %
            % Similar to the implementation of GammaDistribution class, we
            % can pass only one parameter to the constructor. In that case
            % 'a' is the array of priors used to initialize components.
            % -------------------------------------------------------------------------------------------
            
            switch nargin
                % Default constructor, container will have a single default
                % GammaDistribution
                % ------------------------------------------------------------------
                case 0
                    numDistributions = 1;

                    obj.distributions = repmat(GammaDistribution(), numDistributions, 1); % Preallocate
                    obj.distributions(1) = GammaDistribution();


                % If it is a scalar then we create a container that has
                % that many default GammaDistribution instances, but if
                % it is an array of priors (or a single
                % GammaDistribution) then it is used to initialize the
                % parameters of the components.
                % ------------------------------------------------------------------
                case 1
                    % 'a' is a number
                    if Utility.isSingleNumber(a)
                        numDistributions = a;
                        obj.distributions = repmat(GammaDistribution(), numDistributions, 1); % Preallocate

                        for i=1:numDistributions
                            obj.distributions(i) = GammaDistribution();
                        end
                    % 'a' is a list of priors
                    elseif Utility.areAllInstancesOf(a, 'GammaDistribution')
                        numDistributions = length(a);
                        obj.distributions = repmat(GammaDistribution(), numDistributions, 1); % Preallocate

                        for i=1:numDistributions
                            obj.distributions(i) = GammaDistribution(a(i));
                        end
                    else
                        error(['Error in class ' class(obj) ': Invalid argument passed in.']);
                    end

                case {2, 3}
                    % Prepare first 'a' and 'b' - they hold the parameter
                    % values and are set for both types of constructor

                    % OPTION 1
                    % a: []
                    % b: []
                    if Utility.isArray(a) && Utility.isArray(b)
                        if length(a) ~= length(b)
                            error(['Error in ' class(obj) ': Arrays of values' ...
                                ' for a and b should be of the same size.']);
                        end

                        numDistributions = length(a);
                    end
    
    
                    % OPTION 2
                    % a: scalar
                    % b: []
                    if isscalar(a) && Utility.isArray(b)
                        numDistributions = length(b);

                        a = repmat(a, 1, numDistributions); % repmat(a, 1, N);
                    end
    
                    % OPTION 3
                    % a: []
                    % b: scalar
                    if Utility.isArray(a) && isscalar(b)
                        numDistributions = length(a);

                        b = repmat(b, 1, numDistributions); % repmat(a, 1, N);
                    end
                    % Container with only 1 distribution
                    % a: scalar
                    % b: scalar
                    if isscalar(a) && isscalar(b)
                        numDistributions = 1;
                    end

                    % Priors
                    hasPriors = nargin == 3 && ~Utility.isNaN(priors);

                    if hasPriors && length(priors) > 1 && length(priors) ~= length(a)
                        error(['Error in ' class(obj) ': Arrays of values' ...
                            ' for a and priors should be of the same size.']);
                    end

                    obj.distributions = repmat(GammaDistribution(), numDistributions, 1); % Preallocate

                    if nargin == 2 % This is necessary because when nargin == 2 then 'priors' is not defined
                        priors = NaN;
                    end

                    for i = 1:numDistributions
                        if ~hasPriors
                            obj.distributions(i) = GammaDistribution(a(i), b(i));
                        else
                            % "isscalar(priors)" is equivalent to "length(priors) == 1", even though priors is an object!
                            obj.distributions(i) = Utility.ternaryOpt(isscalar(priors), ...
                            @() GammaDistribution(a(i), b(i), priors), @() GammaDistribution(a(i), b(i), priors(i)));
    
                        end
                    end


                % Constructor with 4 parameters: all distributions have the
                % same parameter values for 'a', 'b' and 'prior' parameter. 'a'
                % and 'b' are scalar in this case because the number of
                % distributions is defined, so it could result in conflicting
                % value if I allow them (one or both) to be arrays. This can be
                % added later as an extension if needed. The potential
                % conflict between the number of priors in 'priors' and
                % 'numDistributions' results in an error.
                % a: scalar
                % b: scalar
                % ------------------------------------------------------------------
                case 4
                    hasPriors = ~Utility.isNaN(priors); % 'priors' is set

                    validParameters = Utility.isSingleNumber(a) && Utility.isSingleNumber(b) && ...
                        Utility.isSingleNumber(numDistributions) && ...
                        (~hasPriors || Utility.areAllInstancesOf(priors, 'GammaDistribution') && ...
                        (length(priors) == numDistributions) || isscalar(priors));
     
                    if ~validParameters
                        error(['Error in ' class(obj) ': Invalid input arguments.'])
                    end
                    
                    obj.distributions = repmat(GammaDistribution(), numDistributions, 1); % Preallocate
                    
                    for i = 1:numDistributions
                        if ~hasPriors
                            obj.distributions(i) = GammaDistribution(a, b);
                        else
                            obj.distributions(i) = Utility.ternaryOpt(isscalar(priors), @() GammaDistribution(a, b, priors), ...
                            @() GammaDistribution(a, b, priors(i)));
                        end
                    end

                otherwise
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

            % Default for 'inc' is false
            if nargin == 4
                inc = false;
            end

            obj.distributions(idx).updateParameters(a, b, inc);
        end

        %% All distibutions methods
        function obj = updateAllDistributionsParams(obj, a, b, inc)
            if nargin < 3
                error(['Error in class ' class(obj) ': Too few arguments passed.']);
            end
            % Default for 'inc' is false
            if nargin == 3
                inc = false;
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

        function value = get.ExpectationLnP(obj)
            value = cell(1, obj.Size);
                for i = 1:obj.Size
                    value{i} = obj.distributions(i).ExpectationLnP;
                end
        end

        function value = get.ExpectationLnPC(obj)
            value = 0;
            for i = 1:obj.Size
                value = value + obj.distributions(i).ExpectationLnP;
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