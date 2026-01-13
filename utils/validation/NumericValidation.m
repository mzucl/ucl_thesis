classdef NumericValidation
    % NUMERICVALIDATION Helper class for numeric input validation
    %
    % Provides safe numeric checks that avoid errors for non-numeric or
    % object inputs.

    methods (Static)
        %% Safe NaN check for numeric inputs
        function res = isNaN(obj)
            arguments
                obj
            end
            % NOTE:
            %   Built-in isnan throws an error for non-numeric inputs or
            %   instances of user-defined classes. This method safely
            %   returns false in those cases.
            res = isnumeric(obj) && isnan(obj);
        end

        %% Check for finite numeric scalar
        function res = isFiniteNumericScalar(x)
            arguments
                x
            end
            % Returns true only if x is:
            %   - numeric
            %   - scalar
            %   - finite (not NaN or Inf)
            res = isnumeric(x) && isscalar(x) && isfinite(x);
        end
    end
end
