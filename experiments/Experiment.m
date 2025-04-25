classdef Experiment
    properties
        modelName
        K
        views

        logFileName
        numModelSelectionRuns
        numStabilityRuns
        elboRecalcInterval

        maxIter
        tol
        doRotation
        
        result
    end



    methods
        function obj = Experiment(modelName, K, views, ...
                logFileName, numModelSelectionRuns, numStabilityRuns, elboRecalcInterval, enableLogging, inputValidation, ...
                maxIter, tol, doRotation)
            CustomError.validateNumberOfParameters(nargin, 3, 9);

            obj.modelName = modelName;
            obj.K = K;
            obj.views = views;

            % Default values for optional parameters
            rc = RunConfig.getInstance();
            obj.numModelSelectionRuns = rc.numModelSelectionRuns;
            obj.numStabilityRuns = rc.numStabilityRuns;
            obj.elboRecalcInterval = rc.elboRecalcInterval;

            obj.maxIter = Utility.getConfigValue('Optimization', 'DEFAULT_MAX_ITER');
            obj.tol = Utility.getConfigValue('Optimization', 'DEFAULT_TOL');
            obj.doRotation = false;

            obj.logFileName = '';

            % Override default values with provided arguments if available
            if nargin >= 4 && ~isempty(logFileName)
                obj.logFileName = logFileName;
            end
            if nargin >= 5 && ~isempty(numModelSelectionRuns)
                obj.numModelSelectionRuns = numModelSelectionRuns;
            end
            if nargin >= 6 && ~isempty(numStabilityRuns)
                obj.numStabilityRuns = numStabilityRuns;
            end
            if nargin >= 7 && ~isempty(elboRecalcInterval)
                obj.elboRecalcInterval = elboRecalcInterval;
            end
            if nargin >= 8 && ~isempty(enableLogging)
                % Set `rc.enableLogging`
                rc.enableLogging = enableLogging;
            end
            if nargin >= 9 && ~isempty(inputValidation)
                % Set `rc.inputValidation`
                rc.inputValidation = inputValidation;
            end
            if nargin >= 10 && ~isempty(maxIter)
                obj.maxIter = maxIter;
            end
            if nargin >= 11 && ~isempty(tol)
                obj.tol = tol;
            end
            if nargin >= 12 && ~isempty(doRotation)
                obj.doRotation = doRotation;
            end
        end
        


        function bestOverallModel = run(obj)
            if ~isempty(obj.logFileName)
                if ~exist('logs', 'dir')
                    mkdir('logs');
                end
                diary(['logs/', obj.logFileName, '.txt']); % Start logging
            end
            
            % Average convergence iteration of the best models selected during the 
            % `numModelSelectionRuns` loop.
            convItAvg = 0;
            
            tic;
            % The best model out `numStabilityRuns` best models
            maxOverallElbo = -Inf;
            bestOverallModel = NaN;
            bestOverallConvIt = NaN;
            
            for s = 1:obj.numStabilityRuns 
                maxElbo = -Inf;
                bestModel = NaN;
                bestConvIt = NaN;
            
                for i = 1:obj.numModelSelectionRuns
                    switch lower(obj.modelName)
                        case 'sgfa'
                            model = SGFA(obj.views, obj.K, obj.maxIter, obj.tol, obj.doRotation);
                        case 'gfa'
                            model = GFA(obj.views, obj.K);
                        otherwise
                            CustomError.raiseError('InputCheck', CustomError.ERR_UNKNOWN_MODEL);
                    end
                    
                        
                    [elboVals, convIt] = model.fit(obj.elboRecalcInterval);
                
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
                fprintf('The best model of run %d/%d converged in %d iterations.\n', s, obj.numStabilityRuns, bestConvIt);
            end
            
            elapsedTime = toc;
            fprintf('\n\n### Experiment Summary ###\n');
            fprintf('Elapsed time: %.4f seconds\n', elapsedTime);
            fprintf('Average number of convergence iterations: %d\n', round(convItAvg / obj.numStabilityRuns));
            fprintf('Best overall model converged in %d iterations\n', bestOverallConvIt);
            
            diary off; % End logging          
        end
    end
end