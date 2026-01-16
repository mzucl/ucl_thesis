classdef ConfigUtilsTest < matlab.unittest.TestCase
    % MATRIXVALIDATIONTEST Unit tests for MatrixValidation helper methods
    %
    % Example usage:
    %   runtests('MatrixValidationTest')
    % WARNING: These test are testing access to the `config.txt` file. If the values/descriptions there are changed they will fail. 

    methods(Test)
        function test_getValue(testCase)
            val = ConfigUtils.getValue('Distribution', 'DEFAULT_GAMMA_A');
            testCase.verifyEqual(val, 10^-14);
        end

        function test_getDescription(testCase)
            val = ConfigUtils.getDescription('Distribution', 'DEFAULT_GAMMA_A');
            testCase.verifyEqual(val, 'Shape parameter for default Gamma prior');
        end
    end
end