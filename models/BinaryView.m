classdef BinaryView < handle
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

        mu              % [D x 1] Gaussian

        xi              % [D x N] Variational parameters

        bound           % Jaakkola or Bohning bound
    

        Z
        K

        % CONSTANT (don't change after initialization) dependent properties
        D
        N
    end
    

    properties(Access = private, Constant)
        SETTINGS = ModelSettings.getInstance();
    end

    properties(Access = private)
        expInit
        cache = struct(...
            'C', NaN, ...
            'G', NaN, ...
            'H', NaN);

        cacheFlags = false(1, 3); % I invalidate all entries in the cache at the same time
                                  % but I don't store valid values at the
                                  % same time, so I need a flag for each
                                  % entry.
                                  % Hardcoded for optimization purposes!
    end

    properties (Dependent)
        % Matrices that depend on \xi
        % TODO: Check if I need them??? I can access them directly through
        % the bound
        C
        G
        H
    end


    methods
        %% Constructors
        function obj = BinaryView(data, Z, K, featuresInCols)
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
            obj.alpha = GammaContainer("SD", obj.K.Val, SGFAGroup.SETTINGS.DEFAULT_GAMMA_A, SGFAGroup.SETTINGS.DEFAULT_GAMMA_B);
            
            %                  dim, mu,    cov,  priorPrec
            obj.mu = Gaussian(obj.D, 0, eye(obj.D), 10^3);

            %                         type, size_, cols,   dim,        mu, cov, priorPrec
            obj.W = GaussianContainer("DS", obj.D, false, obj.K.Val, randn(obj.K.Val, obj.D));

            % Model initialization - second part
            % The first update equation is for W, so we need to initialize
            % everything that is used in those two equations and those
            % initilizations are given below.
            %   obj.T.expInit
            %   obj.alpha.ExpInit
            % ----------------------------------------------------------------      
            obj.alpha.setExpInit(repmat(1e-1, obj.K.Val, 1));

            % Performed only once
            obj.qConstantUpdates();
        end





        %% Update methods
        function obj = qConstantUpdates(obj)
            % alpha.a
            obj.alpha.updateAllDistributionsA(obj.alpha.prior.a + obj.D/2);
        end

        function obj = qWUpdate(obj, it)
            alphaExp = Utility.ternary(it == 1, obj.alpha.getExpInit(true), obj.alpha.E_Diag);
            tauExp = Utility.ternary(it == 1, obj.tau.getExpInit(), obj.tau.E);

            covNew = Utility.matrixInverse(tauExp * obj.Z.E_XXt + alphaExp);
            muNew = covNew * tauExp * obj.Z.E * (obj.X.X - obj.mu.E)';
           
            obj.W.updateDistributionsParameters(muNew, covNew);
        end

        function obj = qAlphaUpdate(obj)
            bNew = obj.alpha.prior.b + 1/2 * obj.W.E_SNC;
            obj.alpha.updateAllDistributionsB(bNew);
        end

        function obj = qMuUpdate(obj)
            

            covNew = (1/(obj.mu.priorPrec + obj.N * obj.tau.E)) * eye(obj.D);
            muNew = obj.tau.E * (covNew * (obj.X.X * ones(obj.N, 1) - obj.W.E * obj.Z.E * ones(obj.N, 1)));
            
            obj.mu.updateParameters(muNew, covNew);
        end

        function value = getExpectationLnW(obj)
            value = obj.W.E_SNC' * obj.alpha.E;
            value = -1/2 * value + obj.D/2 * (obj.alpha.E_LnX - obj.K.Val * log(2*pi));
        end

        % TODO (medium): The value can be reused from the qTauUpdate method
        function value = getExpectationLnPX(obj)
            k = (obj.N * obj.D)/2;
            value = k * (obj.tau.E_LnX - log(2 * pi)) - obj.tau.E/2 * (obj.X.Tr_XtX - 2 * obj.mu.E_Xt * sum(obj.X.X, 2) + ...
                obj.N * obj.mu.E_XtX - 2 * sum(sum((obj.W.E * obj.Z.E) .* (obj.X.X - obj.mu.E))) + ...
                sum(sum(obj.W.E_XtX' .* obj.Z.E_XXt)));
        end



        %% Dependent properties
        function value = get.C(obj)
            if ~obj.cacheFlags(1)
                obj.cache.E = obj.a / obj.b;
                obj.cacheFlags(1) = true;
            end
            value = obj.cache.E;
        end

        function value = get.G(obj)
        end

        function value = get.H(obj)
        end
    end
end