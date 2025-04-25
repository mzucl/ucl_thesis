classdef BaseView < handle
    properties
        X       % ViewHandler instance
        W       % GaussianContainer
        alpha   % GammaContainer
        mu      % Gaussian
        Z       % Latent variables
        K       % Number of components
        D       % Data dimensionality
        N       % Number of data points
    end

    methods
        function obj = BaseView(data, Z, K, featuresInCols)
            if nargin < 4
                featuresInCols = true;
            end


            % obj.X = ViewHandler(data, featuresInCols);
            obj.Z = Z;
            obj.K = K;

            obj.D = obj.X.D;
            obj.N = obj.X.N;

            obj.alpha = GammaContainer( ...
                "SD", ...
                obj.K.Val, ...
                Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_A'), ...
                Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_B') ...
                );

            obj.mu = Gaussian(obj.D, 0, eye(obj.D), 1e3);
        end

        function value = getExpectationLnPW(obj)
            value = -0.5 * (obj.W.E_SNC' * obj.alpha.E) + obj.D / 2 * (obj.alpha.E_LnX - obj.K.Val * log(2 * pi));
        end
    end
end
