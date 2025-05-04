classdef StandardScaler
    properties
        mean_ % Stores the mean of each feature
        std_  % Stores the standard deviation of each feature
    end
    
    methods
        function obj = StandardScaler()
            obj.mean_ = [];
            obj.std_ = [];
        end
        
        % Method to fit the scaler on training data
        function obj = fit(obj, data)
            % [NOTE] Data is stored in [DxN] format
            obj.mean_ = mean(data, 2);        % Mean of each feature (columns)
            obj.std_ = std(data, 1, 2);       % Std deviation of each feature (columns)
            obj.std_(obj.std_ == 0) = eps;    % Avoid division by zero
        end
        
        % Method to transform the data using the fitted scaler
        function scaledData = transform(obj, data)
            scaledData = (data - obj.mean_) ./ obj.std_;  % Standard scaling
        end
        
        function originalData = inverse_transform(obj, scaledData)
            originalData = scaledData .* obj.std_ + obj.mean_;  % Reverse scaling
        end
    end
end
