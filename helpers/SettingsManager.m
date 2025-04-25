classdef SettingsManager
    properties (Access = private)
        % Default settings
        defaultSettings
        % Current settings (can be modified)
        currentSettings
    end
    
    methods (Access = private)
        % Private constructor to prevent direct instantiation
        function obj = SettingsManager()
            % Set default settings
            obj.defaultSettings = struct('DEBUG', false, 'VALIDATE', true);
            % Copy default settings to current settings
            obj.currentSettings = obj.defaultSettings;
        end
    end
    
    methods (Static)
        % Singleton pattern to get the settings manager instance
        function obj = getInstance()
            persistent instance
            if isempty(instance)
                disp('Creating new instance of SettingsManager');
                instance = SettingsManager();
            end
            obj = instance;
        end
    end
    
    methods
        % Reset the current settings to default values
        function resetToDefaults(obj)
            disp('Resetting to default settings');
            obj.currentSettings = obj.defaultSettings;
        end
        
        % Get the current settings
        function settings = getSettings(obj)
            disp('Getting current settings');
            settings = obj.currentSettings;
        end
        
        % Update a specific setting
        function setSetting(obj, field, value)
            if isfield(obj.currentSettings, field)
                disp(['Setting ' field ' updated to ' num2str(value)]);
                obj.currentSettings.(field) = value;
            else
                error('Invalid setting: %s', field);
            end
        end
    end
end
