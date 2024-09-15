classdef (Abstract) Bound
    properties
        a
    end
    
    % These methods are defined as 'Static' because they are general utility 
    % functions hat are not tied to any specific instance of the 'Bound' object.
    methods (Static)
        function value = sigma(x)
            value = exp(x)./(exp(x) + 1);
        end

        function value = lse(x)
            value = log(exp(x) + 1);
        end
    end

    methods
        function obj = Bound(a)
            % Default constructor
            if nargin == 0
                obj.a = 0;
            else
                obj.a = a;
            end
        end

        function value = c(obj)
            value = Bound.lse(obj.a);
        end
    
        function value = g(obj)
            % TODO (medium): The second term is same as c, maybe this can 
            % be used for improving performance?
            value = exp(obj.a - Bound.lse(obj.a));
        end
    end

    methods(Abstract)
        h(obj)
    end
end
