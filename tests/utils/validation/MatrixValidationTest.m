classdef MatrixValidationTest < matlab.unittest.TestCase
    % MATRIXVALIDATIONTEST Unit tests for MatrixValidation helper methods
    %
    % Tests for:
    %   - NumericValidation.isNumericVector
    %   - (other methods can be added here in the future)
    %
    % Example usage:
    %   runtests('MatrixValidationTest')
    %   runtests('MatrixValidationTest', 'Tag', 'isNumericVector')
    methods(Test, TestTags = {'isNumericVector'})
        function testRowVector(testCase)
            x = [1, 2, 3];
            testCase.verifyTrue(MatrixValidation.isNumericVector(x));
        end

        function testColumnVector(testCase)
            x = [1; 2; 3];
            testCase.verifyTrue(MatrixValidation.isNumericVector(x));
        end

        function testScalar(testCase)
            x = 5;
            testCase.verifyFalse(MatrixValidation.isNumericVector(x));
        end

        function test2DMatrix(testCase)
            x = [1, 2; 3, 4];
            testCase.verifyFalse(MatrixValidation.isNumericVector(x));
        end

        function test3DTensor(testCase)
            x = rand(2,3,4);
            testCase.verifyFalse(MatrixValidation.isNumericVector(x));
        end

        function testCharRowVector(testCase)
            x = ['a','b','c'];
            testCase.verifyFalse(MatrixValidation.isNumericVector(x));
        end

        function testCellColumnVector(testCase)
            x = {1;2;3};
            testCase.verifyFalse(MatrixValidation.isNumericVector(x));
        end

        function testEmptyVector(testCase)
            x = [];
            testCase.verifyFalse(MatrixValidation.isNumericVector(x));
        end
    end
end