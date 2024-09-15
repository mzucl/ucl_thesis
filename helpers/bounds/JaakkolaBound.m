classdef JaakkolaBound < Bound
    methods
        function value = lambda(obj)
            value = (1./(1 + exp(-obj.a)) - 0.5) ./ (2 * obj.a);
        end

        function result = g(obj)
            result =  0.5 + 2 * obj.lambda() .* obj.a;
        end

        function result = h(obj)
            result = 2 * obj.lambda();
        end
    end
end