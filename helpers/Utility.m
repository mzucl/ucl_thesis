% TODO (medium): Write tests for methods in this file
classdef Utility
    methods (Static)
        function res = isArray(x)
            res = ~isscalar(x) && ismatrix(x) && numel(size(x)) == 2 && (size(x, 1) == 1 || size(x, 2) == 1);
        end

        function res = isMatrix(x)
            res = ~isscalar(x) && ismatrix(x) && ~Utility.isArray(x);
        end

        function res = areAllInstancesOf(arr, className)
            res = all(arrayfun(@(x) isa(x, className), arr));
        end

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
                error(['Error in class ' class(obj) ': Input must be a square matrix.']);
            end
  
            eigenvalues = eig(matrix);
            
            % Check if all eigenvalues are non-negative
            res = all(eigenvalues >= -1e-10);
        end

        function res = isValidCovarianceMatrix(matrix)
            res = Utility.isSymmetricMatrix(matrix) && Utility.isSemiPositiveDefinite(matrix);
        end

        function invA = matrixInverse(A)
            % Compute the inverse of matrix A using LU decomposition
            if ~Utility.isSquareMatrix(A)
                error(['Error in class ' class(obj) ': Matrix must be square for inversion.']);
            end
            
            % Compute the inverse using LU decomposition
            [L, U, P] = lu(A);
            invA = U \ (L \ P);
        end

        function result = ternary(cond, valTrue, valFalse)
            if cond
                result = valTrue;
            else
                result = valFalse;
            end
        end

        function A = generateRandomSPDMatrix(n)
            R = randn(n);

            A = R' * R;
        end
    end
end
