classdef JaakkolaBound < Bound
    methods(Access = private)
        function value = lambda(obj)
            value = (1./(1 + exp(-obj.xi)) - 0.5) ./ (2 * obj.xi);
        end

        function result = g(obj)
            result =  0.5 + 2 * obj.lambda() .* obj.xi;
        end

        function result = h(obj)
            result = 2 * obj.lambda();
        end
    end
end