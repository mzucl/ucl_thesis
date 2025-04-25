classdef SGFAGroup < handle
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

        tau             % [scalar] Gamma         

        Z
        K

        % CONSTANT (don't change after initialization) dependent properties
        D
        N
    end
    


    methods
        %% Constructors
        function obj = SGFAGroup(data, Z, K, featuresInCols)
            % obj@BaseView();
            % disp('SGFAGroup');
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
            %  type, size_, a, b, prior
            obj.alpha = GammaContainer( ...
                "SD", ...
                obj.K.Val, ...
                Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_A'), ...
                Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_B'));
            
            %                  dim, mu,    cov,  priorPrec
            obj.mu = Gaussian(obj.D, 0, eye(obj.D), 10^3);
            
            %
            tauPrior = Gamma( ...
                Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_A'), ...
                Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_B'));

            obj.tau = Gamma(tauPrior);

            %                         type, size_, cols,   dim,        mu, cov, priorPrec
            obj.W = GaussianContainer("DS", obj.D, false, obj.K.Val, randn(obj.K.Val, obj.D));

            % Model initialization - second part
            % The first update equation is for W, so we need to initialize
            % everything that is used in those two equations and those
            % initilizations are given below.
            %   obj.T.expInit
            %   obj.alpha.ExpInit
            % ----------------------------------------------------------------
            obj.tau.setExpInit(1000);        
            obj.alpha.setExpInit(repmat(1e-1, obj.K.Val, 1));

            % Performed only once
            obj.qConstantUpdates();
        end





        %% Update methods
        function obj = qConstantUpdates(obj)
            % alpha.a
            obj.alpha.updateAllDistributionsA(obj.alpha.prior.a + obj.D/2);
            % tau.a
            obj.tau.updateA(obj.tau.prior.a + (obj.N * obj.D)/2);
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

        function obj = qTauUpdate(obj)
            bNew = obj.tau.prior.b + 1/2 * ( ...
                obj.X.Tr_XtX - 2 * obj.mu.E_Xt * sum(obj.X.X, 2) + ...
                obj.N * obj.mu.E_XtX - 2 * sum(sum((obj.W.E * obj.Z.E) .* (obj.X.X - obj.mu.E))) + ...
                sum(sum(obj.W.E_XtX' .* obj.Z.E_XXt)));

            obj.tau.updateB(bNew);
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
    end
end