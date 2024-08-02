classdef GammaDistribution < handle
    properties
        a  
        b
    end
    
    properties (Dependent)
        Expectation
        Variance
        Value
    end
    
    methods
        function obj = GammaDistribution(a, b)
            % Optional parameters: a, b
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

        function value = get.Value(obj)
            value = gamrnd(obj.a, obj.b);
        end
    end
end