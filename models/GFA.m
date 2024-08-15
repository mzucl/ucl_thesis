% TODO (medium): Create a base class for all models with optimization params
% and some other stuff
classdef GFA < handle
    properties
        K               % Number of latent dimensions/principal components
    
        N               % Number of observations

        M               % Number of groups

        Z               % [K x N] GaussianDistributionContainer [size: N; for each latent variable zn]
        
        views           % An array of GFAGroup instances

        % Optimization parameters
        maxIter
        tol
        % They all share Z, also it shouldn't be a copy it should be a
        % reference!!!
    end
    
    % TODO (high): Implement
    % methods(Access = private)
    %     function [M, N] = validateSources(obj, idx)
    % 
    %     end
    % end

    methods
        % [NOTE] We need to deal with the form of the datasets here (e.g.
        % are they in table format, or we need to import
        % data from .csv files. Also, we should check if all datasets have the
        % same number of observations.
        %
        % [path, featureInCols: true/false, data] triplet can be used to
        % describe sources, where if path is None data is in 'data' and
        % viceversa

        % [NOTE] For now we have 2 datasets passed in as matrices
        %% Constructors
        function obj = GFA(data, K, maxIter, tol)
            if nargin < 2
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
            end

            % TODO (very high!): Implementent method below as soon as the
            % toy example with 2 views starts working
            % [obj.M, obj.N] = obj.validateSources(...)
            obj.M = length(data);
            obj.N = size(data{1}, 2); % Data passed in is DxN

            obj.K = DoubleWrapper(K);

            % Set parameters to default values that will be updated if value is
            % provided
            obj.maxIter = Constants.DEFAULT_MAX_ITER;
            obj.tol = Constants.DEFAULT_TOL;

            % Optional parameters: maxIter, tol
            switch nargin
                case 3
                    obj.maxIter = maxIter;
                case 4
                    obj.maxIter = maxIter;
                    obj.tol = tol;
            end

            %% Model setup and initialization
            % Z
            % ------------------------------------------------------ %
            % Initialize the model - set random values for the mean
            % This means we will run the update equation for W first and
            % that we should set some values for all the moments that are
            % in those update equations.
            initZMu = randn(obj.K.Val, 1);
            zPrior = GaussianDistribution(initZMu, eye(obj.K.Val));
            obj.Z = GaussianDistributionContainer(obj.N, zPrior, true);

            % Initialize views
            obj.views = GFAGroup.empty(obj.M, 0);

            for i = 1:obj.M
                obj.views(i) = GFAGroup(data{i}, obj.Z, obj.K, false); % featuresInCols = false;
            end
        end



        %% Update methods
        % obj.Z is GaussianDistributionContainer(cols = true)
        function obj = qZUpdate(obj)
            disp('qZUpdate');
            % Update covariance
            % All latent variables have the same covariance
            covNew = zeros(obj.K.Val);
            for m = 1:obj.M
                covNew = covNew + obj.views(m).W.E_Ct * ...
                    obj.views(m).T.E_Diag * obj.views(m).W.EC;
            end
            covNew = Utility.matrixInverse(eye(obj.K.Val) + covNew);
            obj.Z.updateAllDistributionsCovariance(covNew);

            % Update mu
            for n = 1:obj.N
                muNew = zeros(obj.K.Val, 1);
                for m = 1:obj.M
                    muNew = muNew + obj.views(m).W.E_Ct * ...
                        obj.views(m).T.E_Diag * obj.views(m).X.getObservation(n);
                end
                muNew = covNew * muNew;
                obj.Z.updateDistributionMu(n, muNew);
            end
        end

        function obj = qWUpdate(obj, it)
            disp('qWUpdate');
            for i = 1:obj.M
                obj.views(i).qWUpdate(it);
            end
        end

        function obj = qAlphaUpdate(obj)
            disp('qAlphaUpdate');
            for i = 1:obj.M
                obj.views(i).qAlphaUpdate();
            end
        end

        function obj = qTauUpdate(obj)
            disp('qTauUpdate');
            for i = 1:obj.M
                obj.views(i).qTauUpdate();
            end
        end



        %% fit() and ELBO
        function [elboVals, it, resArr] = fit(obj)
            elboVals = -Inf(1, obj.maxIter);
            resArr = cell(1, obj.maxIter);
        
            for it = 1:obj.maxIter
                obj.removeFactors();
                obj.qWUpdate(it);
                obj.qZUpdate();
                obj.qAlphaUpdate();
                obj.qTauUpdate();

                [currElbo, res] = obj.computeELBO();

                resArr{it} = res;

                % CHECK: ELBO has to increase from iteration to iteration
                if it ~= 1 && currElbo < elboVals(it - 1)
                    fprintf(2, 'ELBO decreased in iteration %d\n', it);
                end 

                elboVals(it) = currElbo;

                if it ~= 1
                    disp(['======= ELBO increased in iteration ', num2str(it), ' by: ', ...
                        num2str(currElbo - elboVals(it - 1))]);
                end

                % Check for convergence
                if it ~= 1 && abs(currElbo - elboVals(it - 1)) / abs(currElbo) < obj.tol
                    disp(['Convergence at iteration: ', num2str(it)]);
                    elboVals = elboVals(1:it); % cut the -Inf values at the end
                    resArr = resArr(1:it);
                    break;
                end
            end
        end

        function [elbo, res] = computeELBO(obj)
            % DEBUG
            res = {};

            % q
            qW = 0; qAlpha = 0; qTau = 0;
            % p
            pX = 0; pW = 0; pAlpha = 0; pTau = 0; 

            for i = 1:obj.M
                % p
                pX = pX + obj.views(i).getExpectationLnPX();
                pW = pW + obj.views(i).getExpectationLnW();
                pAlpha = pAlpha + obj.views(i).alpha.E_LnPC;
                pTau = pTau + obj.views(i).T.E_LnPC;
                % q
                qW = qW + obj.views(i).W.HC;
                qAlpha = qAlpha + obj.views(i).alpha.HC;
                qTau = qTau + obj.views(i).T.HC;
            end
            
            % Store to 'res'
            res.qZ = obj.Z.HC; 
            res.qW = qW; 
            res.qAlpha = qAlpha; 
            res.qTau = qTau;
            % p
            res.pX = pX;
            res.pZ = obj.Z.E_LnPC; 
            res.pW = pW;
            res.pAlpha = pAlpha; 
            res.pTau = pTau;

            % DEBUG

            elbo = res.pX + res.pZ + res.pW + res.pAlpha + res.pTau + ... % p(.)
                res.qZ + res.qW + res.qAlpha + res.qTau; % q(.)

            % DEBUG
            res.elbo = elbo;
            % DEBUG
        end


        function obj = removeFactors(obj, threshold)
            if nargin < 2
                threshold = Constants.LATENT_FACTORS_THRESHOLD;
            end
            % Calculate the average of the square of elements for each row of Z
            avgSquare = mean(obj.Z.EC.^2, 2);
        
            removeIdx = find(avgSquare < threshold);

            if isempty(removeIdx)
                return;
            end

            disp(['Removed ', num2str(length(removeIdx)), ' factors!']);
            % Update number of factors
            obj.K.Val = obj.K.Val - length(removeIdx);
        
            % Remove those rows from Z, corresponding columns from W, and elements from alpha
            obj.Z.removeDimensions(removeIdx);
            for m = 1:obj.M
                obj.views(m).alpha.removeDistributions(removeIdx);
                obj.views(m).W.removeDimensions(removeIdx);
            end
        end
    end
end