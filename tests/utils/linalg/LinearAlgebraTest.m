classdef LinearAlgebraTest < matlab.unittest.TestCase
    % CONFIGUTILSTEST Unit tests for the ConfigUtils helper class
    %
    % These tests verify correct retrieval of configuration values and
    % descriptions from the project-level `config.txt` file.
    %
    % Example usage:
    %   runtests('ConfigUtilsTest')
    %
    % WARNING:
    %   These tests depend on the contents of `config.txt`.
    %   If configuration values or descriptions are modified, the tests
    %   may fail and should be updated accordingly.

    properties (TestParameter)
        n = {5, 10, 15};    % matrix size
    end

    methods(Test)
        function test_inverseCholesky(testCase, n)
            A = RandomMatrices.spdMatrix(n);
            
            tol = ConfigUtils.getValue('General', 'TOL');
            testCase.verifyLessThan(norm(inv(A) - LinearAlgebra.inverseCholesky(A), 'fro'), tol);
        end

        function test_inverseLU(testCase, n)
            A = RandomMatrices.intMatrix(n);
            
            tol = ConfigUtils.getValue('General', 'TOL');
            testCase.verifyLessThan(norm(inv(A) - LinearAlgebra.inverseLU(A), 'fro'), tol);
        end

        function test_logDetCholesky(testCase, n)
            A = RandomMatrices.spdMatrix(n);
            
            tol = ConfigUtils.getValue('General', 'TOL');
            testCase.verifyEqual(log(det(A)), LinearAlgebra.logDetCholesky(A), 'AbsTol', tol);
        end

        function test_logDetCholesky_FallbackSVD(testCase)
            testCase.assumeFail('TODO: Add tests for SVD fallback in logDetCholesky.');
        end

        function test_validateAndExtractDiagonal(testCase, n)
            A = RandomMatrices.intMatrix(n) + 1; % non-diagonal matrix
            [isDiagonal, diagElementsOrValue] = LinearAlgebra.validateAndExtractDiagonal(A);

            testCase.verifyEqual(isDiagonal, false);
            testCase.verifyEqual(diagElementsOrValue, []);

            A = RandomMatrices.diagonalMatrix(n); % diagonal matrix
            [isDiagonal, diagElementsOrValue] = LinearAlgebra.validateAndExtractDiagonal(A);

            testCase.verifyEqual(isDiagonal, true);
            testCase.verifyEqual(diagElementsOrValue, diag(A));

            scalar = 5;
            A = scalar * eye(n); % scalar multiple of the identity matrix
            [isDiagonal, diagElementsOrValue] = LinearAlgebra.validateAndExtractDiagonal(A);

            testCase.verifyEqual(isDiagonal, true);
            testCase.verifyEqual(diagElementsOrValue, scalar);
        end

        function test_outerProduct3D(testCase)
            testCase.assumeFail('TODO: ...');
        end

        function test_flatten3DTo2D(testCase)
            testCase.assumeFail('TODO: ...');
        end

        function test_colsTo3D(testCase)
            testCase.assumeFail('TODO: ...');
        end
    end
end