classdef ConfigUtilsTest < matlab.unittest.TestCase
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