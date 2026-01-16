classdef ConfigUtils
    % CONFIGUTILS Utility class for reading configuration constants
    %
    % Provides static methods to:
    %   - Load constants from a configuration file
    %   - Retrieve specific values or descriptions
    %
    % Notes:
    %   - Assumes config file (`config.txt`) is located in the project root
    %   - Caches loaded constants/descriptions in persistent variables.

    methods (Static, Access = private)
        %% Load constants and descriptions from file
        function [constants, descriptions] = loadConstants()
            % CustomError.validateNumberOfParameters(nargin, 0, 0);

            % Use current working directory as project root
            filename = ProjectPaths.configFile();
                
            lines = readlines(filename);
            constants = struct();
            descriptions = struct();
            currentSection = '';

            for i = 1:length(lines)
                line = strtrim(lines(i));

                % Skip empty lines and comments
                if line == "" || startsWith(line, "%")
                    continue;
                end

                % Section header [SECTION]
                if startsWith(line, "[") && endsWith(line, "]")
                    currentSection = extractBetween(line, "[", "]");
                    currentSection = currentSection{1};
                    constants.(currentSection) = struct();
                    descriptions.(currentSection) = struct();
                    continue;
                end

                % Remove trailing semicolon
                if endsWith(line, ";")
                    line = extractBefore(line, ";");
                end

                % Handle key = value # comment
                if contains(line, "=")
                    kv = split(line, "=");
                    key = strtrim(kv(1));
                    rest = strtrim(kv(2));

                    % Split value and comment (on '%')
                    valueStr = rest;
                    commentSplit = regexp(rest, '[%]', 'split');
                    if ~isempty(commentSplit)
                        valueStr = strtrim(commentSplit{1});
                    end

                    % Convert string to numeric, boolean, or keep as string
                    if ismember(valueStr, {'true', 'false'})
                        value = strcmp(valueStr, 'true');
                    elseif contains(valueStr, '^')
                        value = eval(valueStr); % e.g., '10^-14'
                    else
                        value = str2double(valueStr);
                        if isnan(value)
                            value = valueStr; % fallback to string
                        end
                    end

                    constants.(currentSection).(key) = value;

                    % Extract description if present
                    descMatch = regexp(rest, '[%](.*)$', 'tokens');
                    if ~isempty(descMatch)
                        descriptions.(currentSection).(key) = strtrim(descMatch{1}{1});
                    else
                        descriptions.(currentSection).(key) = '';
                    end
                end
            end
        end
    end

    methods (Static)
        %% Get config value by section and key
        function val = getValue(section, key)
            arguments
                section (1,1) string
                key (1,1) string
            end

            persistent configVals

            if isempty(configVals)
                [configVals, ~] = ConfigUtils.loadConstants();
            end

            val = configVals.(section).(key);
        end

        %% Get description of a config entry
        function desc = getDescription(section, key)
            arguments
                section (1,1) string
                key (1,1) string
            end

            persistent descriptions

            if isempty(descriptions)
                [~, descriptions] = ConfigUtils.loadConstants();
            end

            desc = descriptions.(section).(key);
        end
    end
end