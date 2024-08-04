classdef GammaDistribution < handle
    properties
        a  
        b
        prior
    end
    
    properties (Dependent)
        Expectation
        Variance
        H                   % Entropy
        ExpectationLn       % E[ln(x)]; used often
        Value               % Sample from the distribution
        ExpectationLnP      % E[ln(p(x))] wrt to the q(x)
    end
    
    methods
        function obj = GammaDistribution(a, b, prior)
            % Optional parameters: a, b, prior
            obj.prior = NaN;
            switch nargin
                case 0
                    % Default values
                    obj.a = Constants.DEFAULT_GAMMA_A;
                    obj.b = Constants.DEFAULT_GAMMA_B;
                case 1
                    obj.a = a;
                    obj.b = a;
                case 2
                    obj.a = a;
                    obj.b = b;
                case 3
                    obj.a = a;
                    obj.b = b;
                    obj.prior = prior;
                otherwise
                    error(['Error in class ' class(obj) ': Too many arguments passed into the constructor.']);
            end
        end

        %% Update methods: we return obj so we can chain these
        function obj = updateParameters(obj, a, b, inc)
            % Optional parameters: inc
            switch nargin
                case {1, 2}
                    error(['Error in class ' class(obj) ': Too few arguments passed.']);
                case 3
                    % Default for 'inc' is true
                    inc = true;
            end
            obj.a = Utility.ternary(inc, obj.a + a, a);
            obj.b = Utility.ternary(inc, obj.b + b, b);  
        end
        
        function obj = updateA(obj, a, inc)
            % Optional parameters: inc
            switch nargin
                case 1
                    error(['Error in class ' class(obj) ': Too few arguments passed.']);
                case 2
                    % Default for 'inc' is true
                    inc = true;
            end

            obj.a = Utility.ternary(inc, obj.a + a, a);
        end

        function obj = updateB(obj, b, inc)
            % Optional parameters: inc
            switch nargin
                case 1
                    error(['Error in class ' class(obj) ': Too few arguments passed.']);
                case 2
                    % Default for 'inc' is true
                    inc = true;
            end

            obj.b = Utility.ternary(inc, obj.b + b, b);
        end

        % Getters
        function value = get.Expectation(obj)
            value = obj.a / obj.b;
        end
        
        function value = get.Variance(obj)
            value = obj.a / obj.b^2;
        end

        function value = get.H(obj)
            % From MATLAB help
            % Y = gammaln(X) computes the natural logarithm of the gamma function for each element of X.
            % log(X) is the natural logarithm of the elements of X.
            % psi(X) evaluates the psi function (also know as the digamma function) for each element of X.

            value = gammaln(obj.a) - (obj.a - 1) * psi(obj.a) - log(obj.b) + obj.a;
        end

        function value = get.ExpectationLn(obj)
            value = psi(obj.a) - log(obj.b);
        end

        function value = get.ExpectationLnP(obj)
            if ~isa(obj.prior, 'GammaDistribution')
                error(['Error in class ' class(obj) ': Prior must be defined.']);
            end
            value = -gammaln(obj.prior.a) + obj.prior.a * log(obj.prior.b) + ...
                (obj.prior.a - 1) * obj.ExpectationLn - obj.prior.b * obj.Expectation;
        end

        function value = get.Value(obj)
            value = gamrnd(obj.a, obj.b);
        end
    end
end