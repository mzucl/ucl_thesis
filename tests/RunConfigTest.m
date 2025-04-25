classdef RunConfigTest < matlab.unittest.TestCase
    % These settings are copied from the `SettingsManager` class to maintain
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
            ms1 = RunConfig.getInstance();
            testCase.verifyTrue(ms1.input_validation == RunConfigTest.INPUT_VALIDATION_DEFAULT);

            % Update `ms1`
            ms1.input_validation = false;
            ms1.num_model_selection_runs = 10;

            % Get another instance
            ms2 = RunConfig.getInstance();
            
            % Verify that both instances are updated
            testCase.verifyTrue(ms1.input_validation == false);
            testCase.verifyTrue(ms2.input_validation == false);
            testCase.verifyTrue(ms1.num_model_selection_runs == 10);
            testCase.verifyTrue(ms2.num_model_selection_runs == 10);

            % Reset one instance to default values
            ms2.resetToDefaults();

            % Verify that both instances are updated
            testCase.verifyTrue(ms1.input_validation == RunConfigTest.INPUT_VALIDATION_DEFAULT);
            testCase.verifyTrue(ms2.input_validation == RunConfigTest.INPUT_VALIDATION_DEFAULT);
            testCase.verifyTrue(ms1.num_model_selection_runs == RunConfigTest.NUM_MODEL_SELECTION_RUNS_DEFAULT);
            testCase.verifyTrue(ms2.num_model_selection_runs == RunConfigTest.NUM_MODEL_SELECTION_RUNS_DEFAULT);
        end
    end
end
  
