% Settings are defined on the level of the model. These are general
% settings that are relevant for all models

% TODO (high): Set 'Access = private' for the properties and 
% add setters and getters.
classdef ModelSettings < handle
    properties (Access = public)
        % Gamma distribution
        DEFAULT_GAMMA_A
        DEFAULT_GAMMA_B

        % Gaussian distribution
        DEFAULT_GAUSS_MU
        DEFAULT_GAUSS_DIM
        DEFAULT_GAUSS_PRECISION

        % Noise
        DEFAULT_NOISE_VAL

        % Optimization parameters
        DEFAULT_MAX_ITER
        DEFAULT_TOL

        % Latent factors
        LATENT_FACTORS_THRESHOLD

        % Diagonal loading/regulatization
        EPSILON
        
        VALIDATE
        DEBUG
    end

    methods (Access = private)
        function obj = ModelSettings() % private constructor
            % Gamma distribution
            obj.DEFAULT_GAMMA_A = 10^-14;
            obj.DEFAULT_GAMMA_B = 10^-14;
    
            % Gaussian distribution
            obj.DEFAULT_GAUSS_MU = 0;
            obj.DEFAULT_GAUSS_DIM = 1;
            obj.DEFAULT_GAUSS_PRECISION = 1; % spherical prior covariance
    
            % Noise
            obj.DEFAULT_NOISE_VAL = 1e3;
    
            % Optimization parameters
            obj.DEFAULT_MAX_ITER = 5000;
            obj.DEFAULT_TOL = 1e-6; % relative tolerance for convergance
    
            % Latent factors
            obj.LATENT_FACTORS_THRESHOLD = 1e-6; % determine if factor is relevant
    
            % Diagonal loading/regulatization
            obj.EPSILON = 1e-2;
    
            obj.VALIDATE = true; % Validate inputs
            obj.DEBUG = true;
        end
    end

    methods (Static)
        function obj = getInstance()
            persistent instance
            if isempty(instance)
                instance = ModelSettings();
            end
            obj = instance;
        end
    end
end







        