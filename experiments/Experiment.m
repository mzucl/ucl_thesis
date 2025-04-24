classdef Experiment
    properties
        modelName
        K
        views

        elboIterStep
        stabilityRun
        modelSelectionIter

        validateInputs
        debugMode
        logFileName
        
        result
    end

    methods
        % Constructor
        function obj = Experiment(modelName, K, views, ...
                logFileName, elboIterStep, stabilityRun, modelSelectionIter, validateInputs, debugMode)

            CustomError.validateNumberOfParameters(nargin, 3, 9);

            obj.modelName = modelName;
            obj.K = K;
            obj.views = views;

            cfg = ModelSettings.getInstance();

            % Default values for optional parameters
            obj.logFileName = '';
            obj.elboIterStep = cfg.ELBO_ITER_STEP;
            obj.stabilityRun = cfg.STABILITY_RUN;
            obj.modelSelectionIter = cfg.MODEL_SELECTION_ITER;
            obj.validateInputs = cfg.VALIDATE;
            obj.debugMode = cfg.DEBUG;

            % Override default values with provided arguments if available
            if nargin >= 4 && ~isempty(logFileName)
                obj.logFileName = logFileName;
            end
            if nargin >= 5 && ~isempty(elboIterStep)
                obj.elboIterStep = elboIterStep;
            end
            if nargin >= 6 && ~isempty(stabilityRun)
                obj.stabilityRun = stabilityRun;
            end
            if nargin >= 7 && ~isempty(modelSelectionIter)
                obj.modelSelectionIter = modelSelectionIter;
            end
            if nargin >= 8 && ~isempty(validateInputs)
                obj.validateInputs = validateInputs;
            end
            if nargin >= 9 && ~isempty(debugMode)
                obj.debugMode = debugMode;
            end
        end
        

        % Run the experiment
        function bestOverallModel = run(obj)
            % Logging
            if ~isempty(obj.logFileName)
                if ~exist('logs', 'dir')
                    mkdir('logs');
                end
                diary(['logs/', obj.logFileName, '.txt']); % Start logging
            end
            
            % TODO: THIS!!!
            % % Model settings
            % settings = ModelSettings.getInstance();
            % % settings.VALIDATE = false;
            % settings.DEBUG = false;

            % Average convergence iteration of the best models selected during the 
            % `modelSelectionIter` loop runs.
            convItAvg = 0;
            
            tic;
            % The best model out `stabilityRun` best models
            maxOverallElbo = -Inf;
            bestOverallModel = NaN;
            bestOverallConvIt = NaN;
            
            for s = 1:obj.stabilityRun 
                maxElbo = -Inf;
                bestModel = NaN;
                bestConvIt = NaN;
            
                for i = 1:obj.modelSelectionIter
                    switch lower(obj.modelName)
                        case 'sgfa'
                            model = SGFA(obj.views, obj.K);
                        case 'gfa'
                            model = GFA(obj.views, obj.K);
                        otherwise
                            CustomError.raiseError('InputCheck', CustomError.ERR_UNKNOWN_MODEL);
                    end
                    
                        
                    [elboVals, convIt] = model.fit(obj.elboIterStep);
                
                    if elboVals(end) > maxElbo
                        maxElbo = elboVals(end);
                        bestModel = model;
                        bestConvIt = convIt;
                    end
                end
            
                if maxElbo > maxOverallElbo
                    maxOverallElbo = maxElbo;
                    bestOverallModel = bestModel;
                    bestOverallConvIt = bestConvIt;
                end
            
                convItAvg = convItAvg + bestConvIt;
                fprintf('The best model of run %d/%d converged in %d iterations.\n', s, obj.stabilityRun, bestConvIt);
            end
            
            elapsedTime = toc;
            fprintf('\n\n=== Experiment Summary ===\n');
            fprintf('Elapsed time: %.4f seconds\n', elapsedTime);
            fprintf('Average number of convergence iterations: %d\n', round(convItAvg / obj.stabilityRun));
            fprintf('Best overall model converged in %d iterations\n', bestOverallConvIt);
            
            diary off;           
        end
    end
end
