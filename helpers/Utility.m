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
        function res = isSingleNumber(x)
            res = isscalar(x) && isnumeric(x) && ~isnan(x);
        end

        % ismatrix(3) -> true
        function res = isArray(x)
            res = ~isscalar(x) && ismatrix(x) && numel(size(x)) == 2 && (size(x, 1) == 1 || size(x, 2) == 1);
        end

        function res = isMatrix(x)
            res = ~isscalar(x) && ismatrix(x) && ~Utility.isArray(x);
        end

        % This will return true if the obj is NaN, a single instance of the
        % class, or an array of instances of a class
        function res = isNaNOrInstanceOf(obj, className)
            res = isnumeric(obj) && isnan(obj) || Utility.areAllInstancesOf(obj, className);
        end

        % Built-in 'isnan' results in error for instances of a class
        function res = isNaN(obj)
            res = isnumeric(obj) && isnan(obj);
        end

        % Compare obj1 and obj2 that can be NaN or instances of a
        % class
        function res = areEqual(obj1, obj2)
            if Utility.isNaN(obj1) && Utility.isNaN(obj2)
                res = true;
            elseif ~Utility.isNaN(obj1) && Utility.isNaN(obj2) || ...
                    Utility.isNaN(obj1) && ~Utility.isNaN(obj2)
                res = false;
            else 
                res = obj1 == obj2;
            end
        end

        % Returns true even if arr is just a single instance of the class
        function res = areAllInstancesOf(arr, className)
            res = all(arrayfun(@(x) isa(x, className), arr));
        end
        
        function result = ternary(cond, valTrue, valFalse)
            if cond
                result = valTrue;
            else
                result = valFalse;
            end
        end
        
        % Optimized version that doesn't evaluate the unnecessary value,
        % either 'valTrue' or 'valFalse' depending on the 'cond'
        function result = ternaryOpt(cond, valTrueFunc, valFalseFunc)
            if cond
                result = valTrueFunc();
            else
                result = valFalseFunc();
            end
        end



        %% Algebra
        function res = isSquareMatrix(matrix)
            [rows, cols] = size(matrix);
            res = (rows == cols);
        end

        function res = isSymmetricMatrix(matrix)
            res = Utility.isSquareMatrix(matrix) && isequal(matrix, matrix.');
        end

        function res = isPositiveDefinite(matrix)
            % Try to perform Cholesky decomposition
            try
                chol(matrix);
                res = true;
            catch e
                if strcmp(e.identifier, 'MATLAB:posdef')
                    res = false;
                else
                    rethrow(e);
                end
            end
        end

        function res = isSemiPositiveDefinite(matrix)
            if ~Utility.isSquareMatrix(matrix)
                error(['##### ERROR IN THE CLASS' class(obj) ': Input must be a square matrix.']);
            end
  
            eigenvalues = eig(matrix);
            
            % Check if all eigenvalues are non-negative
            res = all(eigenvalues >= -1e-10);
        end

        function res = isValidCovarianceMatrix(matrix)
            % res = Utility.isSymmetricMatrix(matrix) && Utility.isSemiPositiveDefinite(matrix);
            res = Utility.isSemiPositiveDefinite(matrix);
        end

        function [isDiagonal, diagElementsOrValue] = checkAndExtractDiagonal(A)
            isDiagonal = isequal(A, diag(diag(A)));
            
            if isDiagonal
                diagElements = diag(A);

                diagElementsOrValue = Utility.ternary(all(diagElements == diagElements(1)), ...
                    diagElements(1), diagElements);
            else
                diagElementsOrValue = [];
            end
        end

        function invA = choleskyInverse(A)
            if ~Utility.isSquareMatrix(A)
                error(['##### ERROR IN THE CLASS' class(obj) ': Matrix must be square for inversion.']);
            end
            
            % Perform Cholesky decomposition
            try
                L = chol(A, 'lower');
            catch
                error(['##### ERROR IN THE CLASS' class(obj) ': Matrix is not positive definite.']);
            end
            
            invL = inv(L);
            invA = invL' * invL;
        end

        function invA = matrixInverse(A)
            % Compute the inverse of matrix A using LU decomposition
            if ~Utility.isSquareMatrix(A)
                error(['##### ERROR IN THE CLASS' class(obj) ': Matrix must be square for inversion.']);
            end
            
            % Compute the inverse using LU decomposition
            [L, U, P] = lu(A);
            invA = U \ (L \ P);
        end

        function A = generateRandomSPDMatrix(n)
            R = randn(n);

            A = R' * R;
        end

        function A = generateRandomIntMatrix(m, n)
            if nargin < 2 % Square matrix
                n = m;
            end
            
            minValue = 1;
            maxValue = 10;
            
            A = randi([minValue, maxValue], m, n);
        end



        %% Visualization and debugging
        function isIncreasing = isMonotonicIncreasing(arr)
            % Check if the array is monotonically increasing
            isIncreasing = all(diff(arr) > 0);
        end

        function plotStructVariables(resArr)
            numIterations = length(resArr);
            
            % Check the first struct to get the field names
            firstRes = resArr{1};
            fieldNames = fieldnames(firstRes);
            
            numFields = length(fieldNames);
            
            % Determine the number of rows needed for subplots with 2 columns
            numRows = ceil(numFields / 2);
            
            figure;
            
            for i = 1:numFields
                subplot(numRows, 2, i);
                
                data = zeros(1, numIterations);
                
                % Collect data across all iterations
                for j = 1:numIterations
                    data(j) = resArr{j}.(fieldNames{i});
                end
                
                % Plot the data
                plot(1:numIterations, data, 'LineWidth', 1.5);
                
                % Set the title and labels
                title(['Variable: ', fieldNames{i}]);
                xlabel('Iteration');
                ylabel(fieldNames{i});
                grid on;
            end
            
            % Adjust layout for better visibility
            sgtitle('Evolution of ELBO variables over iteration');
        end

        % X is in DxN format
        function [X, Z, W] = generateToyDataset(N, D, K, noiseVariance) % noiseVariance - Variance of the Gaussian noise
            Z = randn(K, N);
            W = randn(D, K); 
        
            % Generate the data matrix X without noise
            X = W * Z;
        
            % Add Gaussian noise
            noise = sqrt(noiseVariance) * randn(D, N);
            X = X + noise;
        end
    
        function hintonDiagram(matrix, t)
            maxWeight = max(abs(matrix(:)));

            for i = 1:size(matrix, 1)
                for j = 1:size(matrix, 2)
                    % Determine the size of the square
                    weight = matrix(i, j);
                    height = sqrt(abs(weight) / maxWeight);
                    width = height;
        
                    % White for positive, black for negative
                    color = [1 1 1] * (weight >= 0);
        
                    rectangle('Position', [j - width / 2, i - height / 2, width, height], 'FaceColor', color, 'EdgeColor', 'none');
                end
            end
            set(gca, 'YDir', 'reverse', 'XAxisLocation', 'top');
            title(t);
        end
    end
end
