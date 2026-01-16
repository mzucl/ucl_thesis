classdef ProjectPaths
    % PROJECTPATHS Centralized project path management
    % Notes:
    %   - Assumes config file (`config.txt`) is located in the project root

    methods (Static)
        function root = projectRoot()
            % Returns absolute path to project root directory
            
            persistent rootDir
            if isempty(rootDir)
                % Anchor: this file lives in <root>/utils folder
                thisFile = mfilename('fullpath');
                rootDir  = fileparts(fileparts(thisFile));
            end
            root = rootDir;
        end

        function path = configFile()
            path = fullfile(ProjectPaths.projectRoot(), 'config.txt');
        end
    end
end