classdef LogicUtils
    % LOGICUTILS Utility class for logical and conditional helpers
    %
    % Provides static methods for:
    %   - Conditional ternary operations
    %   - Lazy evaluation ternary operations
    %
    % Example usage:
    %   result = LogicUtils.ternary(x > 0, 1, -1);
    %   result = LogicUtils.ternaryOpt(x > 0, @() expensiveComputation1(), @() expensiveComputation2());

    methods(Static)
        %% Simple ternary operator
        function result = ternary(cond, valTrue, valFalse)
            % TERNARY Returns valTrue if cond is true, valFalse otherwise
            %
            % Usage:
            %   result = LogicUtils.ternary(condition, valueIfTrue, valueIfFalse);
            arguments
                cond (1,1) logical
                valTrue
                valFalse
            end

            if cond
                result = valTrue;
            else
                result = valFalse;
            end
        end

        %% Lazy evaluation ternary
        function result = ternaryOpt(cond, valTrueFunc, valFalseFunc)
            % TERNARYOPT Returns valTrueFunc() if cond is true, valFalseFunc() otherwise
            %
            % Both valTrueFunc and valFalseFunc are function handles that are
            % ONLY evaluated if selected.
            %
            % Use when computing valTrue or valFalse is expensive or cannot/should not 
            % be evaluated immediately.
            %
            % WARNING: This function is SLOWER than 'ternary'. Do NOT use
            % inside loops!
            %
            % Usage:
            %   result = LogicUtils.ternaryOpt(condition, @() computeTrue(), @() computeFalse());
            %   result = LogicUtils.ternaryOpt(isscalar(priors), @() GaussianDistribution(priors), @() GaussianDistribution(priors(i)));
            arguments
                cond (1,1) logical
                valTrueFunc (1,1) function_handle
                valFalseFunc (1,1) function_handle
            end

            if cond
                result = valTrueFunc();
            else
                result = valFalseFunc();
            end
        end
    end
end