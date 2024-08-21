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
        
        % Check if the array is monotonically increasing
        function isIncreasing = isMonotonicIncreasing(arr)
            isIncreasing = all(diff(arr) > 0);
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

        function res = isRotationMatrix(R)
            % Check if the matrix is square
            [rows, cols] = size(R);
            if rows ~= cols
                res = false;
                return;
            end
            
            % Orthogonality check: R' * R should be close to the identity matrix
            shouldBeIdentity = R' * R;
            identityMatrix = eye(rows);
            orthoCheck = norm(shouldBeIdentity - identityMatrix, 'fro') < 1e-6;
            
            % Determinant check: det(R) should be close to 1
            detCheck = abs(det(R) - 1) < 1e-6;

            res = orthoCheck && detCheck;
        end

        function res = isSemiPositiveDefinite(matrix)
            if ~Utility.isSquareMatrix(matrix)
                error(['##### ERROR IN THE CLASS Utility' ': Input must be a square matrix.']);
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
                error('##### ERROR IN THE CLASS Utility : Matrix must be square for inversion.');
            end
            
            % Perform Cholesky decomposition
            try
                L = chol(A, 'lower');
            catch
                error(['##### ERROR IN THE CLASS Utility' ': Matrix is not positive definite.']);
            end
            
            invL = inv(L);
            invA = invL' * invL;
        end

        function invA = matrixInverse(A)
            % Compute the inverse of matrix A using LU decomposition
            if ~Utility.isSquareMatrix(A)
                error(['##### ERROR IN THE CLASS Utility' ': Matrix must be square for inversion.']);
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

        function R = generateRandomRotationMatrix(n)
            [Q, ~] = qr(randn(n));  % QR decomposition of a random matrix
            
            % Ensure the determinant is 1 (not -1)
            if det(Q) < 0
                Q(:,1) = -Q(:,1);  % Flip the sign of the first column
            end
            
            R = Q;
        end
    end
end
