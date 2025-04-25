classdef (Abstract) BaseView < handle & matlab.mixin.Heterogeneous
    properties
        K

        Z

        W               % [D x K] GaussianContainer
                        %       --- [size: D; each element corresponds to a row of the W matrix]
                        
                        % The prior over W is defined per column (each column
                        % has its own precision parameter), but the update
                        % equations operate on rows. Therefore, W is represented
                        % as a container of size D in a row-wise format.

        alpha           % [K x 1] GammaContainer         
                        %       --- [size: K]

        mu              % [D x 1] Gaussian
    end

    properties (Dependent, SetAccess = private)
        D               % Dimension of the view (#features)
        N               % Number of observations
        X               % ViewHandler instance for storing data
    end

    % Backing variables
    properties (Access = private) 
        D_
        N_
        X_
    end

    methods (Abstract)
        qWUpdate(obj, it) 
        qAlphaUpdate(obj)
        qMuUpdate(obj)
    end

    methods
        function obj = BaseView(data, Z, K, featuresInCols)
            obj.X_ = ViewHandler(data, featuresInCols);

            obj.D_ = obj.X.D;
            obj.N_ = obj.X.N;

            obj.Z = Z;
            obj.K = K;

            %                         type, size_, cols,   dim,        mu, cov, priorPrec
            obj.W = GaussianContainer("DS", obj.D, false, obj.K.Val, randn(obj.K.Val, obj.D));

            %  type, size_, a, b, prior
            obj.alpha = GammaContainer( ...
                "SD", ...
                obj.K.Val, ...
                Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_A'), ...
                Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_B'));
            
            %                  dim, mu,    cov,  priorPrec
            obj.mu = Gaussian(obj.D, 0, eye(obj.D), 10^3);
        end

        function value = getExpectationLnPW(obj)
            value = obj.W.E_SNC' * obj.alpha.E;
            value = -1/2 * value + obj.D/2 * (obj.alpha.E_LnX - obj.K.Val * log(2*pi));
        end

        %% Getters
        function value = get.D(obj)
            value = obj.D_;
        end

        function value = get.N(obj)
            value = obj.N_;
        end

        function value = get.X(obj)
            value = obj.X_;
        end
    end
end
