classdef MatrixValidation
    % MATRIXVALIDATION Helper class for validating matrix properties
    % 
    % Provides static methods to check:
    %   - Numeric vectors
    %   - 2D numeric matrices
    %   - Square, symmetric, singular, rotation matrices
    %   - Positive definite and covariance matrices
    %
    % Example:
    %   MatrixValidation.isNumericVector([1, 2, 3])
    %   MatrixValidation.isRotationMatrix(eye(3))
    
    methods(Static)
        %% Check if input is a numeric vector (1D)
        function res = isNumericVector(x)
            arguments
                x
            end
            % True only for 1D numeric vectors (row or column)
            res = isvector(x) && isnumeric(x) && ~isscalar(x);
        end

        %% Check if array is monotonically increasing
        function isIncreasing = isMonotonicIncreasing(arr)
            arguments
                arr {mustBeNumeric}
            end
            isIncreasing = all(diff(arr) > 0);
        end

        %% Check if input is a numeric 2D matrix
        function res = isNumeric2DMatrix(x)
            arguments
                x
            end
            % True for numeric matrices, excluding 1D vectors
            res = ismatrix(x) && ~MatrixValidation.isNumericVector(x);
        end

        %% Check if a matrix is square
        function res = isSquareMatrix(matrix)
            arguments
                matrix
            end
            [rows, cols] = size(matrix);
            res = (rows == cols);
        end

        %% Check if a matrix is symmetric
        function res = isSymmetricMatrix(matrix)
            arguments
                matrix
            end
            tol = ConfigUtils.getValue('General', 'TOL');
            res = MatrixValidation.isSquareMatrix(matrix) && ...
                  norm(matrix - matrix', 'fro') < tol;
        end

        %% Check if a matrix is singular
        function res = isSingularMatrix(matrix)
            arguments
                matrix
            end
            if ~MatrixValidation.isSquareMatrix(matrix)
                res = true;
                return;
            end
            % Condition number check
            condThreshold = ConfigUtils.getValue('General', 'COND_THRESHOLD');
            res = cond(matrix) > condThreshold;
        end

        %% Check if a matrix is a rotation matrix
        function res = isRotationMatrix(matrix)
            arguments
                matrix
            end
            if ~MatrixValidation.isSquareMatrix(matrix)
                res = false;
                return;
            end
            
            tol = ConfigUtils.getValue('General', 'TOL');
            
            % Orthogonality check: R' * R should be close to the identity matrix
            shouldBeIdentity = matrix' * matrix;
            identityMatrix = eye(size(matrix,1));
            orthoCheck = norm(shouldBeIdentity - identityMatrix, 'fro') < tol;
            
            % Determinant check: det(R) should be close to 1
            detCheck = abs(det(matrix) - 1) < tol;
            
            res = orthoCheck && detCheck;
        end

        %% Check if a matrix is positive definite
        % NOTE: In general, a matrix does not need to be symmetric to be positive 
        % definite (PD) or positive semi-definite (PSD). However, for covariance 
        % matrices, symmetry is required, so these checks enforce symmetry.
        function res = isPositiveDefinite(matrix)
            arguments
                matrix
            end
            if ~MatrixValidation.isSymmetricMatrix(matrix)
                res = false;
                return;
            end
            % Try Cholesky decomposition
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

        %% Check if a matrix is a valid covariance matrix
        function res = isCovarianceMatrix(matrix)
            arguments
                matrix
            end
            % Must be symmetric
            if ~MatrixValidation.isSymmetricMatrix(matrix)
                res = false;
                return;
            end
            % Check positive semi-definiteness
            eigenvalues = eig(matrix);
            res = all(eigenvalues >= 0);
        end
    end
end