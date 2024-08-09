classdef Constants
    properties (Constant)
        % Gamma distribution
        DEFAULT_GAMMA_A = 10^-3; % 10^-14;
        DEFAULT_GAMMA_B = 10^-3; % 10^-14;

        % Gaussian distribution
        DEFAULT_GAUSS_MU = 0;
        DEFAULT_GAUSS_DIM = 1;
        DEFAULT_GAUSS_PRECISION = 10^-1;

        % Default optimization parameters
        DEFAULT_MAX_ITER = 1;
        DEFAULT_TOL = 1e-6;

        % Tau
        DEFAULT_NOISE_VAL = 1e2;
    end
end