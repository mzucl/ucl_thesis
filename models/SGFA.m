classdef SGFA < BaseModel
    methods
        % data = varargin{1}
        % K = varargin{2}
        function obj = SGFA(varargin)
            CustomError.validateNumberOfParameters(nargin, 2, 5);
            obj = obj@BaseModel(varargin{:});

            %                         type, size_, cols, dim,     mu, cov, priorPrec
            obj.Z = GaussianContainer("DS", obj.N, true, obj.K.Val, zeros(obj.K.Val, 1)); % STEP1

            data = varargin{1};
            for m = 1:obj.M
                obj.views(m) = SGFAGroup(data{m}, obj.Z, obj.K, false); % featuresInCols = false;
            end
        end



        %% Abstract methods
        function obj = qZUpdate(obj)
            covNew = zeros(obj.K.Val);
            muNew = zeros(obj.K.Val, obj.N);

            for m = 1:obj.M
                view = obj.views(m);
                covNew = covNew + view.tau.E * view.W.E_XtX;
                muNew = muNew + view.tau.E * view.W.E_Xt * (view.X.X - view.mu.E);
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
                elbo = elbo + view.getExpectationLnPX() + view.getExpectationLnW() ... % p(.)
                    + view.alpha.E_LnP + view.mu.E_LnP + view.tau.E_LnP ... % p(.)
                    + view.W.H + view.alpha.H + view.mu.H + view.tau.H; % q(.)
            end

            elbo = elbo + obj.Z.H + obj.Z.E_LnP;
        end
    end   
end