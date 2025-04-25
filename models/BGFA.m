classdef BGFA < BaseModel
    properties
        Mc              % Number of continous views
        bound           % Bound used for binary views
    end



    methods
        function obj = BGFA(data, Mc, K, bound, maxIter, tol, doRotation)
            CustomError.validateNumberOfParameters(nargin, 3, 7);

            % Set default values
            if nargin < 7, doRotation = false; end
            if nargin < 6, tol = Utility.getConfigValue('Optimization', 'DEFAULT_TOL'); end
            if nargin < 5, maxIter = Utility.getConfigValue('Optimization', 'DEFAULT_MAX_ITER'); end
            if nargin < 4, bound = 'B'; end
            
            obj = obj@BaseModel(data, K, maxIter, tol, doRotation);

            % Validate `Mc`
            if Mc > obj.M
                CustomError.raiseError('Validation', 'The number of continuous views cannot exceed the total number of views.');
            end

            obj.Mc = Mc;
            obj.bound = bound;

            %% Model setup and initialization
            if bound == 'B'
                %                         type, size_, cols, dim,     mu, cov, priorPrec
                obj.Z = GaussianContainer("DS", obj.N, true, obj.K.Val, zeros(obj.K.Val, 1)); % STEP1
            elseif bound == 'J'
                %                         type, size_, cols,   dim,              mu,            -- cov, priorPrec
                obj.Z = GaussianContainer("DD", obj.N, true, obj.K.Val, randn(obj.K.Val, obj.N)); % STEP1
            end

            for m = 1:obj.Mc
                obj.views(m) = SGFAGroup(data{m}, obj.Z, obj.K, false); % featuresInCols = false;
            end

            for m = obj.Mc + 1:obj.M
                obj.views(m) = BinaryView(data{m}, obj.Z, obj.K, false, obj.bound); % featuresInCols = false;
            end
        end



        %% Abstract methods
        function obj = qZUpdate(obj)
            if obj.bound == 'B'
                covNew = zeros(obj.K.Val);
                muNew = zeros(obj.K.Val, obj.N);
    
                for m = 1:obj.Mc
                    view = obj.views(m);
                    covNew = covNew + view.tau.E * view.W.E_XtX;
                    muNew = muNew + view.tau.E * view.W.E_Xt * (view.X.X - view.mu.E);
                end

                for m = obj.Mc + 1:obj.M
                    view = obj.views(m);
                    covNew = covNew + 1/4 * view.W.E_XtX;
                    muNew = muNew + view.W.E_Xt * (view.X.X + view.bound.T - 1/4 * view.mu.E);
                end
    
                covNew = Utility.matrixInverse(eye(obj.K.Val) + covNew);
                muNew = covNew * muNew;
    
                obj.Z.updateDistributionsParameters(muNew, covNew);
            elseif obj.bound == 'J'
                % TODO
            end
        end

        function elbo = computeELBO(obj)
            elbo = 0;
            for m = 1:obj.M
                % p
                view = obj.views(m);
                elbo = elbo + view.getExpectationLnPX() + view.getExpectationLnPW() ... % p(.)
                    + view.alpha.E_LnP + view.mu.E_LnP + ... % p(.)
                    + view.W.H + view.alpha.H + view.mu.H; % q(.)
            end

            % 'tau' for continuous views only
            for m = 1:obj.Mc
                view = obj.views(m);
                elbo = elbo + view.tau.E_LnP + view.tau.H;
            end

            elbo = elbo + obj.Z.H + obj.Z.E_LnP;
        end



        %% Overridden Methods
        % Override `BaseModel` implementation to update `tau` only for continuous views.
        function obj = qTauUpdate(obj)
            for m = 1:obj.Mc % for continuous views only
                obj.views(m).qTauUpdate();
            end
        end

        function obj = qXiUpdate(obj)
            for m = obj.Mc + 1:obj.M % for binary views only
                obj.views(m).qXiUpdate();
            end
        end

        function stepUpdate(obj, it)
            obj.qZUpdate();
            obj.qWUpdate(it);
            obj.qMuUpdate();
            obj.qXiUpdate();
            % if it > 0
            %     obj.updateRotation();
            % end  
            obj.qAlphaUpdate();
            obj.qTauUpdate();
            obj.removeFactors(it);
        end
    end
end