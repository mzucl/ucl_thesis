classdef CustomError < MException
    properties (Constant)
        % Global flag to enable/disable logging
        ENABLE_LOGGING           = true;

        % FileIO
        % Validation
        ERR_NOT_ENOUGH_INPUT_ARG = "Not enough input arguments provided.";
        ERR_TOO_MANY_INPUT_ARG   = "Too many input arguments provided.";



        ERR_INVALID_ARGUMENTS = 'Invalid arguments provided. Please check your input values.';
        ERR_INVALID_PRIOR_PARAMETER = 'The prior parameter is invalid. It must be of type Gamma.';
        ERR_TOO_FEW_ARGUMENTS = 'Insufficient arguments provided. At least the required parameters must be specified.';
        ERR_INVALID_PARAMETERS = 'Invalid parameters. Both parameters must be numeric values.';
        ERR_PARAMETERS_POSITIVE = 'Parameters must be strictly positive values.';
        ERR_EXPECTATION_POSITIVE = 'The expectation must be a strictly positive value.';
        ERR_INVALID_PARAMETER_A = 'Parameter "a" must be a strictly positive value.';
        ERR_INVALID_PARAMETER_B = 'Parameter "b" must be a strictly positive value.';



        ERR_UNKNOWN_MODEL        = "Unknown model name";
    end
    
    methods
        function obj = CustomError(category, className, funcName, message)
            timestamp = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
            errID = sprintf('%s:%s:%s', category, className, funcName);
            fullMessage = sprintf('##### ERROR [%s] %s in %s > %s: %s', ...
                                  category, timestamp, className, funcName, message);

            obj = obj@MException(errID, fullMessage);
            
            % Check if logging is enabled
            if CustomError.ENABLE_LOGGING
                CustomError.logError(fullMessage);
            end
        end
    end

    methods (Static, Access = private)
        function logError(msg)
            fid = fopen('error_log.txt', 'a');
            if fid ~= -1
                fprintf(fid, '%s\n', msg);
                fclose(fid);
            else
                warning('Could not write to error_log.txt');
            end
        end
    end

    methods (Static, Access = public)
        function raiseError(category, message)
            % Extract call stack
            stack = dbstack(1);  % Skip this helper
            if isempty(stack)
                className = 'UnknownClass';
                funcName = 'UnknownFunction';
            else
                fullName = stack(1).name;
                splitName = split(fullName, '.');
        
                if numel(splitName) == 2
                    className = splitName{1};
                    funcName = splitName{2};
                else
                    className = 'main';
                    funcName = fullName;
                end
            end
        
            throw(CustomError(category, className, funcName, message));
        end

        function validateNumberOfParameters(actualNumArgs, minNumArgs, maxNumArgs)
            if actualNumArgs < minNumArgs
                CustomError.raiseError('InputCheck', CustomError.ERR_NOT_ENOUGH_INPUT_ARG);
            elseif actualNumArgs > maxNumArgs
                CustomError.raiseError('InputCheck', CustomError.ERR_TOO_MANY_INPUT_ARG);
            end
        end
    end
end
