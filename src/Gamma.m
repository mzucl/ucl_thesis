% The Gamma distribution is the distribution of the sum of 
% independent exponential random variables. It is characterized by:
% 
% - The shape parameter `a`, representing the number of 
%   exponential random variables being summed.
% - The rate parameter `b` (inverse scale), which governs
%   the rate of exponential decay.
% 
% Note:
% - Both `a` and `b` must always be strictly positive.
%
% Example usage:
%   g = Gamma(2, 3);       % a = 2, b = 3
%   E = g.E;               % expected value
classdef Gamma < handle
    properties (Access = private)
        a        % Shape parameter
        b        % Rate (inverse scale) parameter
        prior    % Prior distribution (Gamma instance)
    end

    properties (Access = private)
        initialExp              % Explicit expectation initialization

        % Stores computed dependent properties to avoid recomputation
        cache = struct(...
            'E', NaN, ...       % Cached expected value
            'H', NaN, ...       % Cached entropy
            'E_LnX', NaN, ...   % Cached E[ln(x)] w.r.t q(x)
            'E_LnP', NaN)       % Cached E[ln(p(x))] w.r.t q(x)

        % Flags to indicate whether cached values are valid
        cacheFlags = false(1, 4)
    end

    % Short property names are intentional to mirror mathematical notation.
    properties (Dependent)
        E           % Expected value
        H           % Entropy
        E_LnX       % E[ln(x)] w.r.t q(x)
        E_LnP       % E[ln(p(x))] w.r.t q(x)
    end

    methods (Access = private)
        function clearCache(obj)
            % CLEARCACHE Invalidates all dependent property cache flags
            obj.cacheFlags = false(1, 4);
        end

        function obj = validateAndSetProperty(obj, value, propName)
            % Validate number of arguments
            CustomError.validateNumberOfParameters(nargin, 3, 3);
    
            % Input validation
            if RunConfig.getInstance().validateInput
                if ~Utility.isSingleNumber(value)
                    CustomError.raiseError(CustomError.ERR_TYPE_INPUT_VALIDATION, ...
                        CustomError.ERR_INVALID_PARAMETER_NUMERICAL);
                end
                if value <= 0
                    CustomError.raiseError(CustomError.ERR_TYPE_INPUT_VALIDATION, ...
                        CustomError.ERR_INVALID_PARAMETER_POSITIVE);
                end
            end
    
            % Set the property
            obj.(propName) = value;
    
            % Clear cache
            obj.clearCache();
        end
    end

    methods
        %% Deep copy
        function newObj = copy(obj)
            % COPY Creates a deep copy of the Gamma object
            newObj = Gamma();
            newObj.a = obj.a;
            newObj.b = obj.b;
            if ~Utility.isNaN(obj.prior)
                newObj.prior = Gamma();
                newObj.prior.a = obj.prior.a;
                newObj.prior.b = obj.prior.b;
            end
        end

        %% Equality operators
        function isEqual = eq(obj1, obj2)
            % EQ Checks equality between two Gamma objects
            if ~Utility.isNaNOrInstanceOf(obj1, 'Gamma') || ~Utility.isNaNOrInstanceOf(obj2, 'Gamma')
                isEqual = false;
                return;
            end
            if Utility.isNaN(obj1) && Utility.isNaN(obj2)
                isEqual = true;
                return;
            elseif xor(Utility.isNaN(obj1), Utility.isNaN(obj2))
                isEqual = false;
                return;
            end
            isEqual = obj1.a == obj2.a && obj1.b == obj2.b;
            if ~isEqual || (Utility.isNaN(obj1.prior) && Utility.isNaN(obj2.prior))
                return;
            end
            isEqual = obj1.prior == obj2.prior;
        end

        function isNotEqual = ne(obj1, obj2)
            % NE Checks inequality between two Gamma objects
            isNotEqual = ~eq(obj1, obj2);
        end

        % Default constructor: initializes with default values
        function obj = Gamma()
            obj.a = Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_A');
            obj.b = Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_B');
            obj.prior = NaN;
        end

        % function obj = Gamma(a, b, prior)
        %     % GAMMA Constructs a Gamma object
        %     %
        %     % Usage:
        %     %   obj = Gamma()               % default a and b
        %     %   obj = Gamma(a)              % sets a=b=a
        %     %   obj = Gamma(a, b)           % sets a and b
        %     %   obj = Gamma(a, b, prior)    % sets a, b, and prior
        %     %   obj = Gamma(priorGammaObj)  % copy constructor with prior
        % 
        %     obj.a = Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_A');
        %     obj.b = Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_B');
        %     obj.prior = NaN;
        % 
        %     switch nargin
        %         case 1
        %             if Utility.areAllInstancesOf(a, 'Gamma')
        %                 obj = a.copy();
        %                 obj.prior = a.copy();
        %             elseif Utility.isSingleNumber(a)
        %                 obj.updateParameters(a, a);
        %             else
        %                 error(['##### ERROR IN CLASS ' class(obj) ': Invalid arguments.']);
        %             end
        %         case {2,3}
        %             obj.updateParameters(a, b);
        %             if nargin == 3
        %                 if RunConfig.getInstance().validateInput && ...
        %                         ~Utility.isNaNOrInstanceOf(prior, 'Gamma')
        %                     error(['##### ERROR IN CLASS ' class(obj) ': Invalid prior parameter.']);
        %                 end
        %                 obj.prior = prior.copy();
        %             end
        %     end
        % 
        %     obj.setInitialExp(obj.E);
        % end

        %% Update methods
        function obj = updateParameters(obj, a, b)
            % UPDATEPARAMETERS Updates both shape (a) and rate (b) parameters
            if RunConfig.getInstance().validateInput
                if nargin < 3
                    error(['##### ERROR IN CLASS ' class(obj) ': Too few arguments.']);
                end
                if ~Utility.isSingleNumber(a) || ~Utility.isSingleNumber(b)
                    error(['##### ERROR IN CLASS ' class(obj) ': Parameters must be numeric.']);
                end
                if (a <= 0 || b <= 0)
                    error(['##### ERROR IN CLASS ' class(obj) ': Parameters must be positive.']);
                end
            end
            obj.a = a;
            obj.b = b;
            obj.clearCache();
        end

        function obj = updateA(obj, a)
            % UPDATEA Updates the shape parameter `a`
            obj = obj.validateAndSetProperty(a, 'a');
        end
        
        function obj = updateB(obj, b)
            % UPDATEB Updates the rate parameter `b`
            obj = obj.validateAndSetProperty(b, 'b');
        end
        
        function obj = setInitialExp(obj, initialExp)
            % SETINITIALEXP Sets the initial expectation value
            obj = obj.validateAndSetProperty(initialExp, 'initialExp');
        end


        %% Getter for initial expectation
        function value = getInitialExp(obj)
            % GETINITIALEXP Returns the initial expectation value
            value = obj.initialExp;
        end

        %% Dependent properties
        function value = get.E(obj)
            % GET.E Returns the expected value
            if ~obj.cacheFlags(1)
                obj.cache.E = obj.a / obj.b;
                obj.cacheFlags(1) = true;
            end
            value = obj.cache.E;
        end

        function value = get.H(obj)
            % GET.H Returns the entropy of the Gamma distribution
            if ~obj.cacheFlags(2)
                obj.cache.H = gammaln(obj.a) - (obj.a - 1) * psi(obj.a) - log(obj.b) + obj.a;
                obj.cacheFlags(2) = true;
            end
            value = obj.cache.H;
        end

        function value = get.E_LnX(obj)
            % GET.E_LNX Returns E[ln(x)] w.r.t q(x)
            if ~obj.cacheFlags(3)
                obj.cache.E_LnX = psi(obj.a) - log(obj.b);
                obj.cacheFlags(3) = true;
            end
            value = obj.cache.E_LnX;
        end

        function value = get.E_LnP(obj)
            % GET.E_LNP Returns E[ln(p(x))] w.r.t q(x)
            if ~obj.cacheFlags(4)
                if isa(obj.prior, 'Gamma')
                    obj.cache.E_LnP = -gammaln(obj.prior.a) + obj.prior.a * log(obj.prior.b) + ...
                        (obj.prior.a - 1) * (psi(obj.a) - log(obj.b)) - obj.prior.b * obj.a / obj.b;
                else
                    obj.cache.E_LnP = NaN;
                end
                obj.cacheFlags(4) = true;
            end
            value = obj.cache.E_LnP;
        end
    end

    methods (Static)
        function obj = fromParameters(a, b, prior)
            % FROMPARAMETERS Creates Gamma object with given a, b, and optional prior
            CustomError.validateNumberOfParameters(nargin, 1, 3);

            if nargin == 1
                b = a;
            end

            obj = Gamma();  % default initialization
            
            
            
            
            obj.updateParameters(a, b);

            if nargin == 3
                if ~isa(prior, 'Gamma') && ~isnan(prior)
                    error('Gamma:InvalidPrior', 'Prior must be a Gamma object or NaN.');
                end
                obj.prior = prior.copy();
            end

            obj.setExpInit(obj.E);
        end

        function obj = fromPrior(prior)
            % FROMPRIOR Creates a Gamma object as a copy of a prior
            if ~isa(prior, 'Gamma')
                error('Gamma:InvalidInput', 'Input must be a Gamma object.');
            end
            obj = prior.copy();
            obj.prior = prior.copy();
        end
    end

      % function obj = Gamma(a, b, prior)
        %     % GAMMA Constructs a Gamma object
        %     %
        %     % Usage:
        %     %   obj = Gamma()               % default a and b
        %     %   obj = Gamma(a)              % sets a=b=a
        %     %   obj = Gamma(a, b)           % sets a and b
        %     %   obj = Gamma(a, b, prior)    % sets a, b, and prior
        %     %   obj = Gamma(priorGammaObj)  % copy constructor with prior
        % 
        %     obj.a = Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_A');
        %     obj.b = Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_B');
        %     obj.prior = NaN;
        % 
        %     switch nargin
        %         case 1
        %             if Utility.areAllInstancesOf(a, 'Gamma')
        %                 obj = a.copy();
        %                 obj.prior = a.copy();
        %             elseif Utility.isSingleNumber(a)
        %                 obj.updateParameters(a, a);
        %             else
        %                 error(['##### ERROR IN CLASS ' class(obj) ': Invalid arguments.']);
        %             end
        %         case {2,3}
        %             obj.updateParameters(a, b);
        %             if nargin == 3
        %                 if RunConfig.getInstance().validateInput && ...
        %                         ~Utility.isNaNOrInstanceOf(prior, 'Gamma')
        %                     error(['##### ERROR IN CLASS ' class(obj) ': Invalid prior parameter.']);
        %                 end
        %                 obj.prior = prior.copy();
        %             end
        %     end
        % 
        %     obj.setInitialExp(obj.E);
        % end



    
    % NOTE:
    % These getters are primarily intended for unit testing and debugging.
    % Core algorithms should rely on dependent properties and public methods.
    methods (Access = {?matlab.unittest.TestCase})
        function value = getA(obj)
            % GETA Returns the shape parameter `a`
            value = obj.a;
        end

        function value = getB(obj)
            % GETB Returns the rate parameter `b`
            value = obj.b;
        end

        function value = getPrior(obj)
            % GETPRIOR Returns the prior
            value = obj.prior;
        end
    end
end