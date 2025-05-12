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
