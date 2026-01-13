classdef TypeValidation
    % TYPEVALIDATION Helper class for validating object types and equality
    %
    % Provides safe utilities for:
    %   - Type checks on arrays and scalars
    %   - Handling NaN values alongside class instances
    %   - Equality checks with NaN-aware semantics

    methods (Static)
        %% Check whether all elements are instances of a given class
        function res = areAllInstancesOf(arr, className)
            arguments
                arr
                className (1,1) string
            end
            % Returns true for:
            %   - a single instance of the class
            %   - an array where all elements are instances of the class
            res = all(arrayfun(@(x) isa(x, className), arr));
        end

        %% Check whether input is NaN or instance(s) of a given class
        function res = isNaNOrInstanceOf(obj, className)
            arguments
                obj
                className (1,1) string
            end
            % True if:
            %   - obj is numeric NaN
            %   - obj is a single instance of className
            %   - obj is an array of instances of className
            res = (isnumeric(obj) && isnan(obj)) || ...
                  TypeValidation.areAllInstancesOf(obj, className);
        end

        %% NaN-safe equality comparison
        function res = areEquivalent(obj1, obj2)
            arguments
                obj1
                obj2
            end
            % Equality semantics:
            %   - NaN == NaN  -> true
            %   - NaN == x    -> false
            %   - otherwise  -> use overloaded '==' operator
            if NumericValidation.isNaN(obj1) && NumericValidation.isNaN(obj2)
                res = true;
            elseif xor(NumericValidation.isNaN(obj1), NumericValidation.isNaN(obj2))
                res = false;
            else
                res = (obj1 == obj2);
            end
        end
    end
end