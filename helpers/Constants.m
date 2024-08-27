classdef Constants
    properties (Constant)
        % Colors
        BLUE = '#DAE8FC';
        PURPLE = '#E1D5E7';
        RED = '#F8CECC'; 

        % Gamma distribution
        DEFAULT_GAMMA_A = 10^-3; % 10^-14;
        DEFAULT_GAMMA_B = 10^-3; % 10^-14;

        % Gaussian distribution
        DEFAULT_GAUSS_MU = 0;
        DEFAULT_GAUSS_DIM = 1;
        DEFAULT_GAUSS_PRECISION = 1; % spherical prior covariance

        % Default optimization parameters
        DEFAULT_MAX_ITER = 100;
        DEFAULT_TOL = 1e-6; % relative tolerance

        DEFAULT_NOISE_VAL = 1e3;

        LATENT_FACTORS_THRESHOLD = 1e-6;

        VALIDATE = false; % Validate inputs in classes

        DEBUG = false;
    end
end