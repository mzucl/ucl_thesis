classdef (Abstract) Bound < handle
    properties
        xi
    end
    
    % These methods are defined as 'Static' because they are general utility 
    % functions hat are not tied to any specific instance of the 'Bound'
    % object, but still are related to the bounds.
    methods (Static)
        function value = sigma(x)
            value = exp(x)./(exp(x) + 1);
        end

        function value = lse(x)
            value = log(exp(x) + 1);
        end
    end



    methods
        function obj = Bound(xi)
            % Default constructor
            if nargin == 0
                obj.xi = 0;
            else
                obj.xi = xi;
            end
        end



        %% Update methods
        function obj = updateXi(newXi)
            obj.xi = newXi;
        end



        %% Getters
        function value = c(obj)
            value = Bound.lse(obj.xi);
        end
    
        function value = g(obj)
            % TODO (medium): The second term is same as c, maybe this can 
            % be used for improving performance?
            value = exp(obj.xi - Bound.lse(obj.xi));
        end

        function value = t(obj)
            value = obj.h() .* obj.xi - obj.g();
        end

        function value = getElboConst(obj)
            value = sum(sum(-obj.c() + obj.xi .* obj.g() - 1/2 * obj.xi.^2 .* obj.h()));
        end
    end

    methods(Abstract)
        h(obj)
    end
end
