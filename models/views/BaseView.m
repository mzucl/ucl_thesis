classdef (Abstract) BaseView < handle & matlab.mixin.Heterogeneous
    properties
        Z
        K
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


            % if nargin < 4
            %     featuresInCols = true;
            % end
            % 
            % 
            % % obj.X = ViewHandler(data, featuresInCols);
            % obj.Z = Z;
            % obj.K = K;
            % 
            % obj.D = obj.X.D;
            % obj.N = obj.X.N;
            % 
            % obj.alpha = GammaContainer( ...
            %     "SD", ...
            %     obj.K.Val, ...
            %     Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_A'), ...
            %     Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_B') ...
            %     );
            % 
            % obj.mu = Gaussian(obj.D, 0, eye(obj.D), 1e3);
        end

        % function value = getExpectationLnPW(obj)
        %     value = -0.5 * (obj.W.E_SNC' * obj.alpha.E) + obj.D / 2 * (obj.alpha.E_LnX - obj.K.Val * log(2 * pi));
        % end




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
