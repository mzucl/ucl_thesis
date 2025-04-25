classdef SGFA < BaseModel
    methods
        % data = varargin{1}
        % K = varargin{2}
        function obj = SGFA(varargin)
            CustomError.validateNumberOfParameters(nargin, 2, 5);
            obj = obj@BaseModel(varargin{:});

            %                         type, size_, cols, dim,     mu, cov, priorPrec
            obj.Z = GaussianContainer("DS", obj.N, true, obj.K.Val, zeros(obj.K.Val, 1)); % STEP1

            % Initialize views
            obj.views = SGFAGroup.empty(obj.M, 0);

            data = varargin{1};
            for i = 1:obj.M
                obj.views(i) = SGFAGroup(data{i}, obj.Z, obj.K, false); % featuresInCols = false;
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



        % Variables with '_' are expectations
        % X_tr and y_tr are used to set the threshold
        function [K_eff, predictions_te] = makePredictions(obj, X_tr, y_tr, X_te)
            Z_ = obj.Z.E;
            K_eff = size(Z_, 1);
            
            W1_ = obj.views(1).W.E;
            W2_ = obj.views(2).W.E;
            mu1_ = obj.views(1).mu.E;
            mu2_ = obj.views(2).mu.E;
            T1_ = obj.views(1).tau.E * eye(obj.D(1));
            
            sigma_Z = Utility.matrixInverse(eye(K_eff) + W1_' * T1_ * W1_);
    
            % Find the best threshold on the train data
            MU_Z = sigma_Z * (W1_' * T1_ * (X_tr' - mu1_));
            
            % [NOTE] Even though this is sigmoid, the value we get is not
            % probability, this is used just to clip it to the [0, 1] range
            predictions_tr = Bound.sigma(W2_ * MU_Z + mu2_);
            [fpr, tpr, thresholds, ~] = perfcurve(y_tr', predictions_tr, 1);
    
            % Calculate G-means
            gMeans = sqrt(tpr .* (1 - fpr));
            [~, idx] = max(gMeans);
            train_best_threshold = thresholds(idx);
            
            % Predictions on the test data
            MU_Z = sigma_Z * (W1_' * T1_ * (X_te' - mu1_));
            predictions_te = Bound.sigma(W2_ * MU_Z + mu2_); 
            predictions_te = predictions_te >= train_best_threshold;
            predictions_te = double(predictions_te');
        end
    end   
end