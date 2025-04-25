classdef RunConfigTest < matlab.unittest.TestCase
    % These settings are copied from the `RunConfig` class to maintain
    % their `private` access level within that class while still allowing 
    % their use here.
    properties (Access = private, Constant)
        INPUT_VALIDATION_DEFAULT = true;
        ENABLE_LOGGING_DEFAULT = true;
        ELBO_RECALC_INTERVAL_DEFAULT = 1;
        NUM_STABILITY_RUNS_DEFAULT = 2;
        NUM_MODEL_SELECTION_RUNS_DEFAULT = 2;
    end

    methods(Test)
        function testAll(testCase)
            rc1 = RunConfig.getInstance();
            testCase.verifyTrue(rc1.inputValidation == RunConfigTest.INPUT_VALIDATION_DEFAULT);

            % Update `rc1`
            rc1.inputValidation = false;
            rc1.numModelSelectionRuns = 10;

            % Get another instance
            rc2 = RunConfig.getInstance();

            % Verify that both instances are updated
            testCase.verifyFalse(rc1.inputValidation);
            testCase.verifyFalse(rc2.inputValidation);
            testCase.verifyEqual(rc1.numModelSelectionRuns, 10);
            testCase.verifyEqual(rc2.numModelSelectionRuns, 10);

            % Reset one instance to default values
            rc2.resetToDefaults();

            % Verify that both instances are reset
            testCase.verifyEqual(rc1.inputValidation, RunConfigTest.INPUT_VALIDATION_DEFAULT);
            testCase.verifyEqual(rc2.inputValidation, RunConfigTest.INPUT_VALIDATION_DEFAULT);
            testCase.verifyEqual(rc1.numModelSelectionRuns, RunConfigTest.NUM_MODEL_SELECTION_RUNS_DEFAULT);
            testCase.verifyEqual(rc2.numModelSelectionRuns, RunConfigTest.NUM_MODEL_SELECTION_RUNS_DEFAULT);
        end
    end
end