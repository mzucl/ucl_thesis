classdef SGFAGroup < BaseView
    properties
        tau             % [scalar] Gamma         
    end

  

    methods
        %% Constructors
        function obj = SGFAGroup(data, Z, K, featuresInCols)
            CustomError.validateNumberOfParameters(nargin, 3, 4);

            % Set default values
            if nargin < 4, featuresInCols = true; end

            obj@BaseView(data, Z, K, featuresInCols);

            %% Model setup and initialization
            tauPrior = Gamma( ...
                Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_A'), ...
                Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_B'));

            obj.tau = Gamma(tauPrior);

            % Model initialization – second stage
            % The first update is for W, so we initialize all components required for
            % the W and alpha update equations here. These include:
            %   obj.tau.expInit
            %   obj.alpha.expInit
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

        % TODO (medium): The value can be reused from the qTauUpdate method
        function value = getExpectationLnPX(obj)
            k = (obj.N * obj.D)/2;
            value = k * (obj.tau.E_LnX - log(2 * pi)) - obj.tau.E/2 * (obj.X.Tr_XtX - 2 * obj.mu.E_Xt * sum(obj.X.X, 2) + ...
                obj.N * obj.mu.E_XtX - 2 * sum(sum((obj.W.E * obj.Z.E) .* (obj.X.X - obj.mu.E))) + ...
                sum(sum(obj.W.E_XtX' .* obj.Z.E_XXt)));
        end
    end
end