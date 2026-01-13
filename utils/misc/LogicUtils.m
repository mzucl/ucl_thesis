classdef LogicUtils

    methods (Static)
        function result = ternary(cond, valTrue, valFalse)
            if cond
                result = valTrue;
            else
                result = valFalse;
            end
        end
        
        % Optimized version that doesn't evaluate the unnecessary value,
        % either 'valTrue' or 'valFalse' depending on the 'cond'
        % Also, slower compared to 'ternary' for simple stuff, so it should
        % be used just in case both true and false statements can't/shouldn't be
        % evaluated or outside loops.
        %
        % EXAMPLE: Utility.ternaryOpt(isscalar(priors), @() GaussianDistribution(priors), @() GaussianDistribution(priors(i)));
        function result = ternaryOpt(cond, valTrueFunc, valFalseFunc)
            if cond
                result = valTrueFunc();
            else
                result = valFalseFunc();
            end
        end
    end
end