classdef LinearAlgebra
    % LINEARALGEBRA Helper class for advanced linear algebra operations
    %
    % Provides static methods for:
    %   - Matrix inversion (Cholesky and LU)
    %   - Log-determinant computation
    %   - Diagonal validation/extraction
    %   - 3D tensor computations
    %
    % Notes:
    %   - Some methods validate inputs internally for convenience.
    %   - Useful for numerical stability and tensor manipulations.

    methods(Static)
        %% Cholesky-based matrix inverse
        % VALIDATE flag checked internally for convenience
        function invMatrix = inverseCholesky(matrix)
            arguments
                matrix
            end
            % [NOTE] Validation is performed here instead of outside loops
            if RunConfig.getInstance().validateInput
                if MatrixValidation.isSingularMatrix(matrix)
                    CustomError.raiseError(CustomError.ERR_TYPE_ARG_VALIDATION, ...
                        CustomError.ERR_NON_SINGULAR_ARG);
                elseif ~MatrixValidation.isSymmetricMatrix(matrix)
                    CustomError.raiseError(CustomError.ERR_TYPE_ARG_VALIDATION, ...
                        CustomError.ERR_NON_SYMMETRIC_ARG);
                end
            end
            
            % Perform Cholesky decomposition
            try
                L = chol(matrix, 'lower');
            catch e
                if strcmp(e.identifier, 'MATLAB:posdef')
                    CustomError.raiseError(CustomError.ERR_TYPE_ARG_VALIDATION, ...
                        CustomError.ERR_NON_PD_ARG);
                else
                    rethrow(e);
                end
            end
            
            invL = inv(L);
            invMatrix = invL' * invL;
        end

        %% LU-based matrix inverse
        % [NOTE] Square, non-singular matrices. No symmetry or PD required (unlike inverseCholesky).
        function invMatrix = inverseLU(matrix)
            arguments
                matrix
            end
            if RunConfig.getInstance().validateInput && MatrixValidation.isSingularMatrix(matrix)
                CustomError.raiseError(CustomError.ERR_TYPE_ARG_VALIDATION, ...
                        CustomError.ERR_NON_SINGULAR_ARG);
            end

            [L, U, P] = lu(matrix);
            invMatrix = U \ (L \ P);
        end

        %% Log-determinant using Cholesky (or SVD fallback)
        function logDetA = logDetCholesky(A)
            arguments
                A
            end
            try
                % Requires positive definiteness
                L = chol(A, 'lower');
                logDetA = 2 * sum(log(diag(L)));
            catch
                % Fallback for positive semi-definite matrices
                warning('Matrix is not positive definite. Falling back to SVD.');
                s = svd(A);
                tol = ConfigUtils.getValue('General', 'TOL'); % avoid log(0)
                s = s(s > tol);
                logDetA = sum(log(s));
            end
        end

        %% Check if matrix is diagonal and extract diagonal
        function [isDiagonal, diagElementsOrValue] = validateAndExtractDiagonal(A)
            arguments
                A
            end
            % Check if A is diagonal
            isDiagonal = isequal(A, diag(diag(A)));
            
            if isDiagonal
                diagElements = diag(A);
                diagElementsOrValue = LogicUtils.ternary(all(diagElements == diagElements(1)), ...
                    diagElements(1), diagElements);
            else
                diagElementsOrValue = [];
            end
        end

        %% 3D tensor: outer product of columns
        function result = outerProduct3D(MU)
            arguments
                MU
            end
            [n, m] = size(MU);
            MU_reshaped = reshape(MU, [n, 1, m]);
            result = bsxfun(@times, MU_reshaped, permute(MU_reshaped, [2, 1, 3]));
        end

        %% Flatten 3D tensor to 2D matrix
        function A = flatten3DTo2D(M)
            arguments
                M
            end
            [K, ~, D] = size(M);
            A = reshape(M, [K*K, D]);
        end

        %% Convert matrix columns to 3D diagonal tensor
        function result = colsTo3D(M)
            arguments
                M
            end
            [rows, cols] = size(M);
            I = repmat(eye(rows), [1, 1, cols]);
            M_3D = reshape(M, [rows, 1, cols]);
            result = bsxfun(@times, I, M_3D);
        end
    end
end
