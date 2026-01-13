classdef GFAGroup < handle
    properties         
        X               % ViewHandler instance to hold the data

        W               % [D x K] GaussianContainer      
                        %       --- [size: D; for each row in W matrix]

                        % Prior over W is defined per columns (each column
                        % has its own precision parameter, but update
                        % equations are defined by rows, so we are
                        % representing W as a size D container in a row
                        % format.

        alpha           % [K x 1] GammaContainer         
                        %       --- [size: K]

        T               % [D x 1] GammaContainer         
                        %       --- [size: D]

        Z
        K

        % CONSTANT (don't change after initialization) dependent properties
        D
        N
    end
   


    methods
        %% Constructors
        function obj = GFAGroup(data, Z, K, featuresInCols)
            if nargin < 3
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
            elseif nargin < 4
                featuresInCols = true;
            end

            obj.X = ViewHandler(data, featuresInCols);
            obj.Z = Z;
            obj.K = K;

            % Dependent properties
            obj.D = obj.X.D;
            obj.N = obj.X.N;

            %% Model setup and initialization
            %                          type, size_, a, b, prior
            obj.alpha = GammaContainer( ...
                "SD", ...
                obj.K.Val, ...
                Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_A'), ...
                Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_B') ...
                );

            %                      type, size_, a, b, prior
            obj.T = GammaContainer( ...
                "SD", ...
                obj.D, ...
                Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_A'), ...
                Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_B') ...
                );

            %                         type, size_, cols,   dim,   mu, cov, priorPrec
            obj.W = GaussianContainer("DD", obj.D, false, obj.K.Val, randn(obj.K.Val, obj.D));

            % Model initialization - second part
            % The first update equation is for W, so we need to initialize
            % everything that is used in those two equations and those
            % initilizations are given below.
            %   obj.T.expInit
            %   obj.alpha.ExpInit
            % ----------------------------------------------------------------
            obj.T.setExpInit(1000 * ones(obj.D, 1));        
            obj.alpha.setExpInit(repmat(1e-1, obj.K.Val, 1));

            % Performed only once
            obj.qConstantUpdates();
        end





        %% Update methods
        function obj = qConstantUpdates(obj)
            % alpha.a
            obj.alpha.updateAllDistributionsA(obj.alpha.prior.a + obj.D/2);
            % T.a
            obj.T.updateAllDistributionsA(obj.T.prior.a + obj.N/2);
        end


        function obj = qAlphaUpdate(obj)
            bNew = obj.alpha.prior.b + 1/2 * obj.W.E_SNC;
            obj.alpha.updateAllDistributionsB(bNew);
        end

        %% HUGE BUG HERE! Can't be separated like this under the expectation!!!

        function obj = qTauUpdate(obj)
            bNew = obj.T.prior.b + 1/2 * diag( ...
                obj.X.XXt - 2 * obj.W.E * obj.Z.E * obj.X.X' ...
                + obj.W.E * obj.Z.E_XXt * obj.W.E');

            obj.T.updateAllDistributionsB(bNew);
        end
    

        function obj = qWUpdate(obj, it)
            alphaExp = LogicUtils.ternary(it == 1, obj.alpha.getExpInit(true), obj.alpha.E_Diag);
            TExp = LogicUtils.ternary(it == 1, obj.T.getExpInit(), obj.T.E);

            % Update cov
            T_3D = reshape(TExp, 1, 1, []);

            ZZt_3D = repmat(obj.Z.E_XXt, 1, 1, obj.D);
            
            alpha_3D = repmat(alphaExp, 1, 1, obj.D);

            cov_inv = pagemtimes(ZZt_3D, T_3D) + alpha_3D;

            covNew = arrayfun(@(i) LinearAlgebra.inverseLU(cov_inv(:,:,i)), ...
                1:size(cov_inv, 3), 'UniformOutput', false);

            % Convert cell array to multidimensional matrix
            covNew = cat(3, covNew{:});

            V = reshape(obj.Z.E * obj.X.X' * diag(TExp), obj.K.Val, 1, obj.D); % Columns of the matrix will be in the third dimension
            muNew = squeeze(pagemtimes(covNew, V));
            
            obj.W.updateDistributionsParameters(muNew, covNew);
        end

        function value = getExpectationLnW(obj)
            % The sum in the eq implemented as a dot product
            value = obj.W.E_SNC' * obj.alpha.E;
            value = -1/2 * value + obj.D/2 * (obj.alpha.E_LnX - obj.K.Val * log(2*pi));
        end

        function value = getExpectationLnPX(obj)
            value = obj.N/2 * obj.T.E_LnX - obj.N * obj.D/2 * log(2 * pi) - 1/2 * ...
                obj.T.E' * diag( ...
                obj.X.XXt ...
                - 2 * obj.W.E * obj.Z.E * obj.X.X' ...
                + obj.W.E * obj.Z.E_XXt * obj.W.E');
        end
    end
end