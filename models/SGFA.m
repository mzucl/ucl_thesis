classdef SGFA < BaseModel
    methods
        function obj = SGFA(data, K, maxIter, tol, doRotation)
            CustomError.validateNumberOfParameters(nargin, 2, 5);
            
            % Set default values
            if nargin < 5, doRotation = false; end
            if nargin < 4, tol = Utility.getConfigValue('Optimization', 'DEFAULT_TOL'); end
            if nargin < 3, maxIter = Utility.getConfigValue('Optimization', 'DEFAULT_MAX_ITER'); end

            obj = obj@BaseModel(data, K, maxIter, tol, doRotation);

            %                         type, size_, cols, dim,     mu, cov, priorPrec
            initZMu = randn(obj.K.Val, 1);

            obj.Z = GaussianContainer("DS", obj.N, true, obj.K.Val, initZMu); % zeros(obj.K.Val, 1)


            for m = 1:obj.M
                obj.views(m) = SGFAGroup(data{m}, obj.Z, obj.K, false); % featuresInCols = false;
            end
        end



        %% Abstract methods
        function obj = qZUpdate(obj, it)
            covNew = zeros(obj.K.Val);
            muNew = zeros(obj.K.Val, obj.N);

            for m = 1:obj.M
                view = obj.views(m);

                tauExp = LogicUtils.ternary(it == 1, view.tau.getExpInit(), view.tau.E);

                covNew = covNew + tauExp * view.W.E_XtX;
                muNew = muNew + tauExp * view.W.E_Xt * (view.X.X - view.mu.E);
            end

            covNew = Utility.matrixInverse(eye(obj.K.Val) + covNew);
            muNew = covNew * muNew;

            obj.Z.updateDistributionsParameters(muNew, covNew);
        end


        function elbo = computeELBO(obj)
            elbo = 0;
            for m = 1:obj.M
                % p
                view = obj.views(m);
                elbo = elbo + view.getExpectationLnPX() + view.getExpectationLnPW() ... % p(.)
                    + view.alpha.E_LnP + view.mu.E_LnP + view.tau.E_LnP ... % p(.)
                    + view.W.H + view.alpha.H + view.mu.H + view.tau.H; % q(.)
            end

            elbo = elbo + obj.Z.H + obj.Z.E_LnP;
        end
    end   
end