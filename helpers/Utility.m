% TODO (medium): Write tests for methods in this file
% [NOTE] isnumeric and isscalar are different methods
%   isnumeric checks if the variable is of a numeric type (even array).
%   isscalar checks if the variable is a single element (1x1 array), regardless of its type.
% 
%   isnumeric(NaN) -> true
%   isnumeric([any element is not numeric]) -> fals [single boolean value]
%   
%   isscalar(NaN) -> true
%   
%   isnan(obj) -> ERROR
%   isnan([]) -> returns elements wise isnan check!
%%


classdef Utility
    methods (Static)




        % ISSINGLENUMBER Check if input is a finite scalar numeric value.
        %
        %   RES = ISSINGLENUMBER(X) returns true if X is a 1x1 numeric value
        %   (double, single, or integer type) that is finite (not NaN or Inf).
        %   Otherwise, it returns false.
        %
        function res = isSingleNumber(x)
            res = isscalar(x) && isnumeric(x) && ~isnan(x);
        end



        % Built-in 'isnan' results in error for instances of a class
        function res = isNaN(obj)
            res = isnumeric(obj) && isnan(obj);
        end

        % This will return true if the obj is NaN, a single instance of the
        % class, or an array of instances of a class
        function res = isNaNOrInstanceOf(obj, className)
            res = isnumeric(obj) && isnan(obj) || Utility.areAllInstancesOf(obj, className);
        end

        % Returns true even if arr is just a single instance of the class
        function res = areAllInstancesOf(arr, className)
            res = all(arrayfun(@(x) isa(x, className), arr));
        end

        % Compare obj1 and obj2 that can be NaN or instances of a class
        % (that overloads '==')
        function res = areEqual(obj1, obj2)
            if Utility.isNaN(obj1) && Utility.isNaN(obj2) % both are 'NaN'
                res = true;
            elseif xor(Utility.isNaN(obj1), Utility.isNaN(obj2)) % one is 'NaN', the other is not
                res = false;
            else 
                res = obj1 == obj2;
            end
        end
    
        function result = importCSV(fileName)
            result = readmatrix(fileName);
        end

        function loss = BCE(yTrue, yPred)
            % yTrue: vector of true labels (0 or 1)
            % yPred: vector of predicted probabilities (between 0 and 1)
            
            % Ensure predictions are in a valid range [epsilon, 1 - epsilon] to avoid log(0)
            epsilon = 1e-15; 
            yPred = max(min(yPred, 1 - epsilon), epsilon);
        
            loss = -mean(yTrue .* log(yPred) + (1 - yTrue) .* log(1 - yPred));
        end

        
    end
end


%% 
% True only for NaN
% isnan(NaN)    % true
% isnan(Inf)    % false
% isnan(5)      % false

% True only for real, finite numbers: Rejects both NaN and ±Inf
% isfinite(NaN)   % false
% isfinite(Inf)   % false
% isfinite(5)     % true


