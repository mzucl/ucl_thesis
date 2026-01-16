classdef MatrixValidationTest < matlab.unittest.TestCase
    % MATRIXVALIDATIONTEST Unit tests for MatrixValidation helper methods
    %
    % Example usage:
    %   runtests('MatrixValidationTest')
    
    methods(Test)
        function test_isNumericVector(testCase)
            x = [1, 2, 3];             % row vector
            testCase.verifyTrue(MatrixValidation.isNumericVector(x));
                
            x = [1; 2; 3];             % column vector
            testCase.verifyTrue(MatrixValidation.isNumericVector(x));

            x = ['a', 'b', 'c'];       % char row vector
            testCase.verifyFalse(MatrixValidation.isNumericVector(x));

            x = 5;                     % scalar
            testCase.verifyFalse(MatrixValidation.isNumericVector(x));
        
            x = [1, 2; 3, 4];          % 2D matrix
            testCase.verifyFalse(MatrixValidation.isNumericVector(x));

            x = rand(2, 3, 4);         % 3D tensor
            testCase.verifyFalse(MatrixValidation.isNumericVector(x));

            x = {1; 2; 3};             % cell column vector
            testCase.verifyFalse(MatrixValidation.isNumericVector(x));

            x = {1, 2, 3};             % cell row vector
            testCase.verifyFalse(MatrixValidation.isNumericVector(x));

            x = [];                    % empty vector
            testCase.verifyFalse(MatrixValidation.isNumericVector(x));
        end
    end
end