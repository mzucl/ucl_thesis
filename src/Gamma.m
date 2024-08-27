% The Gamma distribution is defined as the distribution of a sum of 
% independent exponential random variables, where the shape parameter 'a'
% corresponds to the number of such variables and the rate parameter 'b' 
% (or the inverse scale parameter) controls the rate of the exponential decay.
%   -> a, b must always be strictly positive 
%
%
%%
classdef Gamma < handle
    properties
        a  
        b
        prior
    end
    
    properties (Access = private)
        expInit % Explicit expectation initialization
        cache = struct(...
            'E', NaN, ...
            'H', NaN, ...
            'E_LnX', NaN, ...
            'E_LnP', NaN);

        cacheFlags = false(1, 4);
    end



    properties (Dependent)
        E
        H          % Entropy
        E_LnX      % E[ln(x)] wrt to q(x)
        E_LnP      % E[ln(p(x))] wrt to the q(x)

        Var
        Val        % Sample from the distribution
    end
    


    methods(Access = private)
        function clearCache(obj)
            obj.cacheFlags = false(1, 4);
        end
    end



    methods
        %% Deep copy and operators overloading
        function newObj = copy(obj)
            newObj = Gamma();
            
            newObj.a = obj.a;
            newObj.b = obj.b;
            
            % Copy the prior (manually)
            if ~Utility.isNaN(obj.prior)
                newObj.prior = Gamma();
                newObj.prior.a = obj.prior.a;
                newObj.prior.b = obj.prior.b;
            end
        end

        function isEqual = eq(obj1, obj2)
            if ~Utility.isNaNOrInstanceOf(obj1, 'Gamma') || ~Utility.isNaNOrInstanceOf(obj2, 'Gamma')
                isEqual = false;
                return;
            end

            % Both are NaN
            if Utility.isNaN(obj1) && Utility.isNaN(obj2)
                isEqual = true;
                return;

            % Only one is NaN
            elseif xor(Utility.isNaN(obj1), Utility.isNaN(obj2))
                isEqual = false;
                return;
            end
                
            % Both are set - compare them!    
            isEqual = obj1.a == obj2.a && obj1.b == obj2.b;

            % If parameters a and b are not equal, they are different objects
            % regardless of the prior!
            %   Second part of the condition is added because obj1.prior ==
            %   obj2.prior won't call 'isEqual' if both are NaN, and we are
            %   relaying on that call for comparing priors.
            if ~isEqual || Utility.isNaN(obj1.prior) && Utility.isNaN(obj2.prior)
                return;
            end

            isEqual = obj1.prior == obj2.prior;
        end

        function isNotEqual = ne(obj1, obj2)
            isNotEqual = ~eq(obj1, obj2);
        end


        
        %% Options for the constructor Gamma
        % ZERO PARAMETERS
        % -> default values for 'a' and 'b'; prior = NaN;
        %
        % 1 PARAMETER: a
        % OPTION 1: 'a' is an instance of Gamma
        %       -> set 'obj' and 'prior' to that value
        % OPTION 2: 'a' is a scalar
        %       -> set 'a' and 'b' to that value
        % 
        % 2 PARAMETERS: a, b
        % -> Gamma(a, b)
        %
        % 3 PARAMETERS: a, b, prior
        % -> same as previous constructor, but the prior is set
        %
        %
        %%
        function obj = Gamma(a, b, prior)
            % Default param values
            obj.a = Constants.DEFAULT_GAMMA_A;
            obj.b = Constants.DEFAULT_GAMMA_B;
            obj.prior = NaN;

            switch nargin
                case 1
                    % If the ONLY parameter is of a type Gamma that is
                    % the prior and we initialize the 'obj' and its prior
                    % using that value
                    if Utility.areAllInstancesOf(a, 'Gamma')
                        obj = a.copy();
                        obj.prior = a.copy();

                    % The one parameter passed in is the value for 'a', set
                    % both 'a' and 'b' to that value
                    elseif Utility.isSingleNumber(a)
                        obj.updateParameters(a, a);
                    else
                        error(['##### ERROR IN THE CLASS ' class(obj) ': Invalid arguments passed.']);
                    end
            
                case {2, 3} % a, b
                    obj.updateParameters(a, b);

                    if nargin == 3 % prior
                        if Constants.VALIDATE && ~Utility.isNaNOrInstanceOf(prior, 'Gamma')
                            error(['##### ERROR IN THE CLASS ' class(obj) ': Invalid prior parameter.']);
                        end
                        obj.prior = prior.copy();
                    end
            end

            % Set initial expectation to the actual expectation
            obj.setExpInit(obj.E);
        end



        %% Update methods
        function obj = updateParameters(obj, a, b)
            if Constants.VALIDATE
                if nargin < 3
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
                end
                if ~Utility.isSingleNumber(a) || ~Utility.isSingleNumber(b)
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Parameters must be numerical values.']);
                end
                if (a <= 0 || b <= 0)
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Parameters are strictly positive number.']);
                end
            end

            obj.a = a;
            obj.b = b;

            % Clear cache
            obj.clearCache();
        end
        
        function obj = updateA(obj, a)
            if Constants.VALIDATE
                if nargin < 2
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
                end
                if ~Utility.isSingleNumber(a)
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Parameter must be a numerical value.']);
                end
                if (a <= 0)
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Parameter a is a strictly positive number.']);
                end
            end

            obj.a = a;
            
            % Clear cache
            obj.clearCache();
        end

        function obj = updateB(obj, b)
            if Constants.VALIDATE
                if nargin < 2
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
                end
                if ~Utility.isSingleNumber(b)
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Parameter must be a numerical value.']);
                end
                if (b <= 0)
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Parameter b is a strictly positive number.']);
                end
            end

            obj.b = b;

            % Clear cache
            obj.clearCache();
        end


        
        %% Setters
        function obj = setExpInit(obj, value)
            if Constants.VALIDATE && value <= 0
                error(['##### ERROR IN THE CLASS ' class(obj) ': Expectation is a strictly positive number.']);
            end
            obj.expInit = value;
        end



        %% Getters
        function value = getExpInit(obj)
            value = obj.expInit;
        end

        

        %% Dependent properties
        function value = get.E(obj)
            if ~obj.cacheFlags(1)
                obj.cache.E = obj.a / obj.b;
                obj.cacheFlags(1) = true;
            end
            value = obj.cache.E;
        end
        
        function value = get.H(obj)
            % From MATLAB help
            % Y = gammaln(X) computes the natural logarithm of the gamma function for each element of X.
            % log(X) is the natural logarithm of the elements of X.
            % psi(X) evaluates the psi function (also know as the digamma function) for each element of X.
            if ~obj.cacheFlags(2)
                obj.cache.H = gammaln(obj.a) - (obj.a - 1) * psi(obj.a) - log(obj.b) + obj.a;
                obj.cacheFlags(2) = true;
            end
            value = obj.cache.H;
        end

        function value = get.E_LnX(obj)
            if ~obj.cacheFlags(3)
                obj.cache.E_LnX = psi(obj.a) - log(obj.b);
                obj.cacheFlags(3) = true;
            end
            value = obj.cache.E_LnX;
        end

        function value = get.E_LnP(obj)
            if ~obj.cacheFlags(4)
                if isa(obj.prior, 'Gamma') % The value is set only when prior is defined
                    obj.cache.E_LnP = -gammaln(obj.prior.a) + obj.prior.a * log(obj.prior.b) + ...
                        (obj.prior.a - 1) * (psi(obj.a) - log(obj.b)) - obj.prior.b * obj.a / obj.b;
                else
                    obj.cache.E_LnP = NaN;
                end
                obj.cacheFlags(4) = true;
            end
            value = obj.cache.E_LnP;
        end

        function value = get.Var(obj)
            value = obj.a / obj.b^2;
        end

        function value = get.Val(obj)
            value = gamrnd(obj.a, obj.b);
        end
    end
end