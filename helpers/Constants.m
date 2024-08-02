classdef Constants
    properties (Constant)
        % Gamma distribution
        DEFAULT_GAMMA_A = 10^-3;
        DEFAULT_GAMMA_B = 10^-3;

        % Gaussian distribution
        DEFAULT_GAUSS_MU = 0;
        DEFAULT_GAUSS_COV = 1;
        DEFAULT_GAUSS_DIM = 1;
        DEFAULT_GAUSS_PRECISION = 10^-3;
    end
end