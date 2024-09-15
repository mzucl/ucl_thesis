classdef LaplaceApproximation < Bound
    methods
        function result = h(obj)
            result = 1./(2 * (cosh(obj.a) + 1));
        end
    end
end