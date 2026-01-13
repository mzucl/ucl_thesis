classdef CustomError < MException
    % CUSTOMERROR Enhanced exception class with logging and helper methods
    %
    %   Provides structured error messages, optional logging, and 
    %   input argument validation helpers.

    properties (Constant)
        % Error type
        ERR_TYPE_ARG_VALIDATION = 'ValidateInput';

        % Error messages
        ERR_NOT_ENOUGH_INPUT_ARG  = 'Not enough input arguments provided.';
        ERR_TOO_MANY_INPUT_ARG    = 'Too many input arguments provided.';
        ERR_INVALID_NUMERIC_ARG   = 'Input argument must be numeric.';
        ERR_INVALID_POSITIVE_ARG  = 'Input argument must be strictly positive.';
        ERR_UNKNOWN_MODEL         = 'Unknown model name';
    end

    methods
        function obj = CustomError(category, className, funcName, message)
            % Construct a CustomError with formatted message and optional logging
            %
            %   OBJ = CUSTOMERROR(CATEGORY, CLASSNAME, FUNCNAME, MESSAGE)

            timestamp = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
            errID = sprintf('%s:%s:%s', category, className, funcName);
            fullMessage = sprintf('##### ERROR [%s] %s %s > %s: %s', ...
                                  category, timestamp, className, funcName, message);

            obj = obj@MException(errID, fullMessage);

            if RunConfig.getInstance().logErrors
                CustomError.logError(fullMessage);
            end
        end
    end

    methods (Static, Access = private)
        function logError(msg)
            % Append error message to log file
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
            % Raise a CustomError using the calling function's stack
            %
            %   raiseError(CATEGORY, MESSAGE)

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
            % Validate number of input arguments
            %
            %   validateNumberOfParameters(ACTUALNUMARGS, MINNUMARGS, MAXNUMARGS)

            if actualNumArgs < minNumArgs
                CustomError.raiseError(ERR_TYPE_ARG_VALIDATION, CustomError.ERR_NOT_ENOUGH_INPUT_ARG);
            elseif actualNumArgs > maxNumArgs
                CustomError.raiseError(ERR_TYPE_ARG_VALIDATION, CustomError.ERR_TOO_MANY_INPUT_ARG);
            end
        end
    end
end
