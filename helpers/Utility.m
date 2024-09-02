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
    properties (Constant)
        SETTINGS = ModelSettings.getInstance();
    end





    methods (Static)
        %% General
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





        %% Check type
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

        

        

        %% Tests with matrices
        function res = isSquareMatrix(matrix)
            [rows, cols] = size(matrix);
            res = (rows == cols);
        end

        function res = isSymmetricMatrix(matrix)
            res = Utility.isSquareMatrix(matrix) && isequal(matrix, matrix.'); % TODO (high): consider using norm(..., 'fro') for this check
        end

        function res = isSingular(matrix)
            res = ~Utility.isSquareMatrix(matrix) || Utility.ternary(det(matrix) == 0, true, false);
        end
 
        function res = isRotationMatrix(matrix)
            if ~Utility.isSquareMatrix(matrix)
                res = false;
                return;
            end
            
            % Orthogonality check: R' * R should be close to the identity matrix
            shouldBeIdentity = matrix' * matrix;
            identityMatrix = eye(rows);

            orthoCheck = norm(shouldBeIdentity - identityMatrix, 'fro') < Constants.TOL;
            
            % Determinant check: det(R) should be close to 1
            detCheck = abs(det(matrix) - 1) < Constants.TOL;

            res = orthoCheck && detCheck;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % [NOTE] In general matrix doesn't need to be symmetric to be PD or
        % PSD, but we use these tests for the covariance matrices are need
        % to be symmetric!
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Cholesky Decomposition: A = LL^T
        % 1. Square
        % 2. Symmetric
        % 3. Positive Definite
        %
        % Matrix cannot be positive definite if it is not symmetric!
        function res = isPositiveDefinite(matrix)
            if ~Utility.isSymmetricMatrix(matrix)
                res = false;
                return;
            end
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

        function res = isPositiveSemiDefinite(matrix)
            if ~Utility.isSymmetricMatrix(matrix)
                res = false;
                return;
            end
            
            eigenvalues = eig(matrix);
            
            % Check if all eigenvalues are non-negative
            res = all(eigenvalues > 0);
        end

        function res = isValidCovarianceMatrix(matrix)
            res = Utility.isPositiveSemiDefinite(matrix);
        end





        %% Matrix inverse
        % For these two methods, the VALIDATE flag is checked within the methods
        % themselves rather than at the point of invocation. This is done for 
        % convenience, especially since these methods are often called within loops.
        function invMatrix = choleskyInverse(matrix)
            if Utility.SETTINGS.VALIDATE 
                if Utility.isSingular(matrix)
                    error(['##### ERROR IN THE CLASS ' mfilename('class') ': Matrix must be non-singular for choleskyInverse.']);
                elseif ~Utility.isSymmetricMatrix(matrix)
                    error(['##### ERROR IN THE CLASS ' mfilename('class') ': Matrix must be symmetric for choleskyInverse.']);
            
                end
            end
            
            % Perform Cholesky decomposition
            try
                L = chol(matrix, 'lower');
            catch
                if strcmp(e.identifier, 'MATLAB:posdef')
                    error(['##### ERROR IN THE CLASS ' mfilename('class') ': Matrix is not positive definite.']);
                else
                    rethrow(e);
                end
            end
            
            invL = inv(L);
            invMatrix = invL' * invL;
        end

        % LU Decomposition: A = LU
        % 1. Square
        % 2. Non-Singular
        % 3. Not Necessarily Symmetric or Positive Definite (unlike
        % Cholesky decomposition)
        function invMatrix = matrixInverse(matrix)
            if Utility.SETTINGS.VALIDATE && Utility.isSingular(matrix)
                    error(['##### ERROR IN THE CLASS ' mfilename('class') ': Matrix must be non-singular for matrixInverse.']);
            end

            [L, U, P] = lu(matrix);
            invMatrix = U \ (L \ P);
        end





        %% Generate different matrices
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

        function matrix = generateRandomBinaryMatrix(m, n)
            matrix = randi([0, 1], m, n);
        end
    
    



        %% Misc
        function logDetA = logDetUsingCholesky(A)
            % Perform Cholesky decomposition
            try
                L = chol(A, 'lower');
            catch
                error(['##### ERROR IN THE CLASS ' mfilename('class') ': Matrix is not positive definite.']);
            end

            logDetA = 2 * sum(log(diag(L)));
        end

        % Check if matrix is diagonal 'isDiagonal' and if is returns the diagonal elements
        % (element in case matrix is a scalar multiple of the identity matrix)
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
    end
end
