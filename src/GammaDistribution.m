% The Gamma distribution is defined as the distribution of a sum of 
% independent exponential random variables, where the shape parameter a
% corresponds to the number of such variables and the rate parameter b 
% (or the inverse scale parameter) controls the rate of the exponential decay.
%   -> a, b must always be strictly positive 
%%
classdef GammaDistribution < handle
    properties
        a  
        b
        prior
    end

    properties(Access = private)
        expInit
    end

    % Private properties for caching
    properties (Access = private)
        Tr_XtX_Cached
        XXt_Cached
    end
    
    % TODO (medium): Naming for 'ExpectationLn' and 'ExpectationLnP' don't follow
    % the 'convention' of adding X, e.g. 'ExpectationXtX' for E[x^Tx] in
    % GaussianDistribution. But given that I am not too happy with those
    % names either, let's leave it like this for now.
    properties (Dependent)
        E
        Var
        H          % Entropy
        E_LnX       % E[ln(x)]
        Val        % Sample from the distribution
        E_LnP      % E[ln(p(x))] wrt to the q(x)
    end
    
    methods
        %% Deep copy and operators overloading
        function newObj = copy(obj)
            newObj = GammaDistribution();
            
            newObj.a = obj.a;
            newObj.b = obj.b;
            
            % Copy the prior (manually)
            if ~Utility.isNaN(obj.prior)
                newObj.prior = GammaDistribution();
                newObj.prior.a = obj.prior.a;
                newObj.prior.b = obj.prior.b;
            end
        end

        function isEqual = eq(obj1, obj2)
            if ~Utility.isNaNOrInstanceOf(obj1, 'GammaDistribution') || ~Utility.isNaNOrInstanceOf(obj2, 'GammaDistribution')
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
        function obj = GammaDistribution(a, b, prior)
            % Optional parameters: a, b, prior
            % -------------------------------------------------------------
            obj.prior = NaN; % It is only set when the number of parameters passed in is 3

            switch nargin
                case 0
                    obj.a = Constants.DEFAULT_GAMMA_A;
                    obj.b = Constants.DEFAULT_GAMMA_B;

                case 1
                    % If the ONLY parameter is of a type GammaDistribution that is
                    % the prior and we initialize the obj and its prior
                    % using that value
                    if Utility.areAllInstancesOf(a, 'GammaDistribution')
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
        % 'inc' (increment) is an optional parameter with default value of false
        function obj = updateParameters(obj, a, b, inc)
            switch nargin
                case {1, 2}
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
                case 3
                    inc = false;
            end
            if ~Utility.isSingleNumber(a) || ~Utility.isSingleNumber(b)
                error(['##### ERROR IN THE CLASS ' class(obj) ': Parameters must be numerical values.']);
            end

            % TODO (high): All these update methods implemented like this
            % will leave the object in half-updated state in the case of an
            % error (e.g. 'a' will be updated and b will throw an error).
            obj.updateA(a, inc);
            obj.updateB(b, inc);
        end
        
        function obj = updateA(obj, a, inc)
            switch nargin
                case 1
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
                case 2
                    inc = false;
            end

            if ~Utility.isSingleNumber(a)
                error(['##### ERROR IN THE CLASS ' class(obj) ': Parameter must be a numerical value.']);
            end
            
            if (inc == true && obj.a + a <= 0 || inc == false && a <= 0)
                error(['##### ERROR IN THE CLASS ' class(obj) ': Parameter a is a strictly positive number.']);
            end

            obj.a = Utility.ternary(inc, obj.a + a, a);
        end

        function obj = updateB(obj, b, inc)
            switch nargin
                case 1
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
                case 2
                    inc = false;
            end

            if ~Utility.isSingleNumber(b)
                error(['##### ERROR IN THE CLASS ' class(obj) ': Parameter must be a numerical value.']);
            end

            if (inc == true && obj.b + b <= 0 || inc == false && b <= 0)
                error(['##### ERROR IN THE CLASS ' class(obj) ': Parameter b is a strictly positive number.']);
            end

            obj.b = Utility.ternary(inc, obj.b + b, b);
        end


        
        %% Setters
        function obj = setExpInit(obj, value)
            if value <= 0
                error(['##### ERROR IN THE CLASS ' class(obj) ': Expectation is a strictly positive number.']);
            end
            obj.expInit = value;
        end



        %% Getters
        % 'expInit' is a private property -> needs a getter for the access
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

        % E[ln \tau]; tau is the name of the random variable (as in the Appendix), not the noise
        % from the models;
        function value = get.E_LnX(obj)
            value = psi(obj.a) - log(obj.b);
        end

        % E[ln p(\tau)]
        function value = get.E_LnP(obj)
            if ~isa(obj.prior, 'GammaDistribution')
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