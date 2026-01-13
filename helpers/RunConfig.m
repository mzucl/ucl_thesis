classdef RunConfig < handle
    properties (Access = private, Constant)
        VALIDATE_INPUT_DEFAULT = true;
        VERBOSE_DEFAULT = true;
        LOG_ERRORS_DEFAULT = true;
        
        
        ELBO_RECALC_INTERVAL_DEFAULT = 1;
        NUM_STABILITY_RUNS_DEFAULT = 2;
        NUM_MODEL_SELECTION_RUNS_DEFAULT = 2;
    end

    properties (Access = public)
        % [Validation / Logs]
        validateInput   % If true, performs input validation (checks dimensions, PSD of covariance matrices, etc.)
        verbose         % If true, prints runtime progress and messages to the command window
        logErrors       % If true, logs errors to a file

        
        
        % [Experiment Setup]
        elboRecalcInterval     % Frequency (in iterations) at which the ELBO is recomputed during training
        numStabilityRuns       % Number of times the experiment is repeated to assess the stability of results
        numModelSelectionRuns  % Number of times the model is trained; the one with the highest ELBO is selected
    end

    methods (Access = private)
        function obj = RunConfig() % private constructor;
            % Initialize all properties to their default values
            obj.setDefaults();
        end

        % Helper method to set all properties to default values
        function setDefaults(obj)
            obj.validateInput = RunConfig.VALIDATE_INPUT_DEFAULT;
            obj.verbose = RunConfig.VERBOSE_DEFAULT;
            obj.logErrors = RunConfig.LOG_ERRORS_DEFAULT;
            

            obj.elboRecalcInterval = RunConfig.ELBO_RECALC_INTERVAL_DEFAULT;
            obj.numStabilityRuns = RunConfig.NUM_STABILITY_RUNS_DEFAULT;
            obj.numModelSelectionRuns = RunConfig.NUM_MODEL_SELECTION_RUNS_DEFAULT;
        end
    end

    methods (Static)
        function obj = getInstance()
            persistent instance
            if isempty(instance)
                instance = RunConfig();
            end
            obj = instance;
        end

        function resetToDefaults()
            obj = RunConfig.getInstance();
            % Reset all properties to their default values
            obj.setDefaults();
        end
    end
end