% Author: Mediha Zukic
% Contact: mediha.zukic.23@alumni.ucl.ac.uk
% Date: 2025-05-13

% Set Python environment if needed
% pyenv('Version', '/path/to/python')

% Import the module if it's not automatically recognized
vis_mod = py.importlib.import_module('visualization');

% Create dataset variable (this needs to be made compatible with Python if it’s a MATLAB struct)
% Assuming dataset is already a py object:
N_FOLDS = int32(10);  % or whatever value
n_reps = int32(10);

% Create the Python Experiment object
GFA_experiment = experiment_mod.Experiment(N_FOLDS, n_reps, dataset, "./Models", "./W", "./tau");

% Call Python methods
GFA_experiment.load_factors();
GFA_experiment.load_taus();


binaryVars = {};
categoricalVars = {};
continuousVars = {};

for v = 1:width(T)
    col = T.(v);
    uniqueVals = unique(col);

    if iscellstr(col) || isstring(col) || iscategorical(col)
        % String-like columns → categorical
        categoricalVars{end+1} = T.Properties.VariableNames{v};
    elseif isnumeric(col)
        if numel(uniqueVals) == 2
            binaryVars{end+1} = T.Properties.VariableNames{v};
        elseif numel(uniqueVals) < 0.05 * height(T) % heuristic threshold
            categoricalVars{end+1} = T.Properties.VariableNames{v};
        else
            continuousVars{end+1} = T.Properties.VariableNames{v};
        end
    end
end

T_binary = T(:, binaryVars);
T_categorical = T(:, categoricalVars);
T_continuous = T(:, continuousVars);
