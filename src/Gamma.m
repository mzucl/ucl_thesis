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
    

    properties (Constant)
        validate = Constants.VALIDATE;
    end


    properties (Access = private)
        expInit % Explicit expectation initialization
        cache = struct(...
            'E', NaN, ...
            'Var', NaN, ...
            'H', NaN, ...
            'E_LnX', NaN, ...
            'E_LnP', NaN);
    end



    properties (Dependent)
        E
        Var
        H          % Entropy
        E_LnX      % E[ln(x)] wrt to q(x)
        E_LnP      % E[ln(p(x))] wrt to the q(x)

        Val        % Sample from the distribution
    end
    


    methods(Access = private)
        function clearCache(obj)
            fields = fieldnames(obj.cache);
            
            for i = 1:length(fields)
                obj.cache.(fields{i}) = NaN;
            end
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



        %% Constructors
        % TODO (high): Use update methods to set the parameters, because in
        % those methods parameter values are validated unlike here
        function obj = Gamma(a, b, prior)
            % Optional parameters: a, b, prior
            % -------------------------------------------------------------
            obj.prior = NaN; % It is only set when the number of parameters passed in is 3

            switch nargin
                case 0
                    obj.a = Constants.DEFAULT_GAMMA_A;
                    obj.b = Constants.DEFAULT_GAMMA_B;

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
                        obj.a = a;
                        obj.b = a;
                    else
                        error(['##### ERROR IN THE CLASS ' class(obj) ': Invalid arguments passed.']);
                    end
            
                case {2, 3}
                    obj.a = a;
                    obj.b = b;
                    if nargin == 3 && ~Utility.isNaN(prior) % 'prior' is passed in
                        obj.prior = prior.copy();
                    end
                otherwise
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Too many arguments passed into the constructor.']);
            end
            % Set initial expectation to the real expectation
            obj.setExpInit(obj.E);
        end



        %% Update methods
        function obj = updateParameters(obj, a, b)
            if Gamma.validate
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
        end
        
        function obj = updateA(obj, a)
            if Gamma.validate
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
        end

        function obj = updateB(obj, b)
            if Gamma.validate
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
        end


        
        %% Setters
        function obj = setExpInit(obj, value)
            if Gamma.validate && value <= 0
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
            value = obj.a / obj.b;
        end
        
        function value = get.Var(obj)
            value = obj.a / obj.b^2;
        end

        function value = get.H(obj)
            % From MATLAB help
            % Y = gammaln(X) computes the natural logarithm of the gamma function for each element of X.
            % log(X) is the natural logarithm of the elements of X.
            % psi(X) evaluates the psi function (also know as the digamma function) for each element of X.
            value = gammaln(obj.a) - (obj.a - 1) * psi(obj.a) - log(obj.b) + obj.a;
        end

        function value = get.E_LnX(obj)
            value = psi(obj.a) - log(obj.b);
        end

        % TODO (low/medium): Using expressions for obj.E_LnX and obj.E
        % could speed this up.
        function value = get.E_LnP(obj)
            if Gamma.validate && ~isa(obj.prior, 'Gamma')
                error(['##### ERROR IN THE CLASS ' class(obj) ': Prior must be defined.']);
            end
            value = -gammaln(obj.prior.a) + obj.prior.a * log(obj.prior.b) + ...
                (obj.prior.a - 1) * obj.E_LnX - obj.prior.b * obj.E;
        end

        function value = get.Val(obj)
            value = gamrnd(obj.a, obj.b);
        end
    end
end