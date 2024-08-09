%% [NOTE]
% Components inside the container are always independent and they
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
                            % expectation of one component. To access ith
                            % component expectation use obj.Expectation{i}

        ExpectationC        % Array (column vector) of expectations of all components

        H                   % This is a cell array where each entry is an
                            % entropy of one component; 

        HC                  % This is the entropy of the collection (sum of entries in H)

        ExpectationLn       % Cell array where each entry is ExpectationLn

        ExpectationLnC      % Sum of all entries in ExpectationLn

        ExpectationLnP      % Cell array where each entry is ExpectationLnP

        ExpectationLnPC     % The sum of all entries in ExpectationLnP

        Value           
    end

    properties(Access = private)
        expCInit
    end

    methods(Access = private)
        function validateIndex(obj, idx)
            if idx < 1 || idx > obj.Size 
                error(['##### ERROR IN THE CLASS ' class(obj) ': Index out of range.']); 
            end
        end
    end

    %% Options for the constructor GammaDistributionContainer
    % SINGLE PARAMETER
    % ()                        -> single default GammaDistribution object
    % (4)                       -> 4 default GammaDistribution objects
    % (GammaDistribution)       -> single GammaDistribution object with prior
    % ([GammaDistribution])     -> length([]) GammaDistribution objects with
    %                               priors given by []

    % 2 PARAMETERS
    % (a, b)                    -> single GammaDistribution(a, b) object
    % ([a], b)                  -> length(a) GammaDistribution(a(i), b) objects
    % (a, [b])                  -> length(b) GammaDistribution(a, b(i)) objects
    % ([a], [b])                -> length(a) == length(b) GammaDistribution(a(i), b(i)) objects

    % 3 PARAMETERS
    % Same as for the previous constructor, but with the addition of the
    % priors, which can be set to:
    %       NaN                 -> all components have NaN prior
    %       GammaDistribution   -> all components have same prior
    %       [GammaDistribution] -> all components have their own prior
    
    % 4 PARAMETERS
    % (a, b, priors, numDistributions)  -> all components (# of them in numDistributions)
    %                                       have the same 'a' and 'b', and
    %                                       priors defined by 'priors'
    %                                       which can be a single object or
    %                                       NaN, or an array of objects in
    %                                       which case all distributions
    %                                       have different priors

    methods
        %% Constructor
        function obj = GammaDistributionContainer(a, b, priors, numDistributions)
            % Left 'numDistributions' as a last parameter because it can be
            % infered from other parameters in some cases. 
            %
            % 'priors' is either single GammaDistribution/NaN in which case all
            % distributions have the same prior, or an array of
            % GammaDistributions - one for each of the components.
            %
            % Similar to the implementation of GammaDistribution class, we
            % can pass a single parameter to the constructor. In that case
            % 'a' is the array of priors used to initialize components and their priors.
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
                    % 'a' is a list of GammaDistribution objects or a
                    % single instance of it
                    elseif Utility.areAllInstancesOf(a, 'GammaDistribution')
                        numDistributions = length(a);
                        obj.distributions = repmat(GammaDistribution(), numDistributions, 1); % Preallocate

                        for i=1:numDistributions
                            obj.distributions(i) = GammaDistribution(a(i));
                        end
                    else
                        error(['##### ERROR IN THE CLASS ' class(obj) ': Invalid argument passed in.']);
                    end

                case {2, 3}
                    % Prepare first 'a' and 'b' - they hold the parameter
                    % values and are set for nargin = 2 and nargin = 3. We
                    % infer the number of components from 'a' and 'b'
                    % values

                    % OPTION 1
                    % a: []
                    % b: []
                    if Utility.isArray(a) && Utility.isArray(b)
                        if length(a) ~= length(b)
                            error(['##### ERROR IN THE CLASS ' class(obj) ': Arrays of values' ...
                                ' for a and b should be of the same size.']);
                        end
                        numDistributions = length(a);
    
                    % OPTION 2
                    % a: scalar
                    % b: []
                    elseif isscalar(a) && Utility.isArray(b)
                        numDistributions = length(b);
                        a = repmat(a, 1, numDistributions);
    
                    % OPTION 3
                    % a: []
                    % b: scalar
                    elseif Utility.isArray(a) && isscalar(b)
                        numDistributions = length(a);
                        b = repmat(b, 1, numDistributions);

                    % OPTION 4
                    % a: scalar
                    % b: scalar
                    elseif Utility.isSingleNumber(a) && Utility.isSingleNumber(b)
                        numDistributions = 1;
                    else 
                        error(['##### ERROR IN THE CLASS ' class(obj) ': Invalid argument passed in.']);
                    end

                    % Priors
                    hasPriors = nargin == 3 && ~Utility.isNaN(priors);

                    if hasPriors && length(priors) > 1 && length(priors) ~= length(a)
                        error(['##### ERROR IN THE CLASS ' class(obj) ': Arrays of values' ...
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
                            obj.distributions(i) = Utility.ternaryOpt(isscalar(priors), ...
                            @() GammaDistribution(a(i), b(i), priors), @() GammaDistribution(a(i), b(i), priors(i)));
    
                        end
                    end


                % Constructor with 4 parameters: all distributions have the
                % same parameter values for 'a', 'b' and 'prior'. 'a'
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
                    hasPriors = ~Utility.isNaN(priors); % 'priors' is not NaN

                    % Validate parameters: separated into multiple lines
                    % for readability
                    validParameters = Utility.isSingleNumber(a) && Utility.isSingleNumber(b) && ...
                        Utility.isSingleNumber(numDistributions);

                    validPriors = Utility.areAllInstancesOf(priors, 'GammaDistribution') && ...
                        (isscalar(priors) || length(priors) == numDistributions);

                    validParameters = validParameters && (~hasPriors || validPriors);
     
                    if ~validParameters
                        error(['##### ERROR IN THE CLASS ' class(obj) ': Invalid input arguments.'])
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
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Invalid input arguments.']);
            end
            % Set initial expectation to the real expectation
            obj.setExpCInit(obj.ExpectationC);
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
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
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
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
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
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Dimensions do not match.']);
                end
            end

            % Check 'b'
            if isscalar(b)
                b = b * ones(1, obj.Size);
            elseif Utility.isArray(b)
                if length(b) ~= obj.Size
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Dimensions do not match.']);
                end
            end

            for i=1:obj.Size
                obj.distributions(i).updateParameters(a(i), b(i), inc);
            end
        end



        %% Setters
        function obj = setExpCInit(obj, value)
            if ~all(value > 0)
                error(['##### ERROR IN THE CLASS ' class(obj) ': Expectation is a strictly positive number.']);
            end
            if length(value) ~= obj.Size
                error(['##### ERROR IN THE CLASS ' class(obj) ': Number of elements in the expectation must be equal to the number of ' ...
                    'components. ']);
            end
            obj.expCInit = value;
        end



        %% Getters
        % 'expInit' is a private property -> needs a getter for the access
        function value = getExpCInit(obj)
            value = obj.expCInit;
        end

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
            value = cell2mat(obj.Expectation)';
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

        function value = get.ExpectationLn(obj)
            value = cell(1, obj.Size);
            for i = 1:obj.Size
                value{i} = obj.distributions(i).ExpectationLn;
            end
        end

        function value = get.ExpectationLnC(obj)
            value = 0;
            for i = 1:obj.Size
                value = value + obj.distributions(i).ExpectationLn;
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