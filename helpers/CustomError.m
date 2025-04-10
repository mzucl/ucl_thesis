classdef CustomError < MException
    properties (Constant)
        % Global flag to enable/disable logging
        ENABLE_LOGGING           = true;

        ERR_NOT_ENOUGH_INPUT_ARG = "Not enough input arguments provided.";
        ERR_TOO_MANY_INPUT_ARG   = "Too many input arguments provided.";
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
    end
end
