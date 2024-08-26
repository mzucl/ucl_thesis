% TODO (medium): Create a base class for all models with optimization params
% and some other stuff
classdef GFA < handle
    properties
        K               % Number of latent dimensions/principal components
    
        N               % Number of observations

        M               % Number of groups

        Z               % [K x N] GaussianContainer [size: N; for each latent variable zn]
        
        views           % An array of GFAGroup instances

        % Optimization parameters
        maxIter
        tol

        % CONSTANT (don't change after initialization) dependent properties
        D % Array containing dimensions for each view

        doRotation
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
        function obj = GFA(data, K, maxIter, tol, doRotation)
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
            obj.doRotation = false;

            % Optional parameters: maxIter, tol
            if nargin > 2
                obj.maxIter = maxIter;
                if nargin > 3
                    obj.tol = tol;
                    if nargin > 4
                        obj.doRotation = doRotation;
                    end
                end
            end

            %% Model setup and initialization
            % Z
            % ------------------------------------------------------ %
            % Initialize the model - set random values for the mean
            % This means we will run the update equation for W first and
            % that we should set some values for all the moments that are
            % in those update equations.
            initZMu = randn(obj.K.Val, 1);

            %                         type, size_, cols, dim,     mu, cov, priorPrec
            obj.Z = GaussianContainer("DS", obj.N, true, obj.K.Val, initZMu);

            % Initialize views
            obj.views = GFAGroup.empty(obj.M, 0);

            for i = 1:obj.M
                obj.views(i) = GFAGroup(data{i}, obj.Z, obj.K, false); % featuresInCols = false;
            end

            obj.D = [obj.views.D];
        end



        %% Update methods
        function obj = qZUpdate(obj)
            covNew = zeros(obj.K.Val);
            sum_WtTX = zeros(obj.K.Val, obj.N); % for mu update

            for m = 1:obj.M
                WtT = obj.views(m).W.E_Xt * obj.views(m).T.E_Diag;
                covNew = covNew + WtT * obj.views(m).W.E;
                sum_WtTX = sum_WtTX + WtT * obj.views(m).X.X;
            end

            covNew = Utility.matrixInverse(eye(obj.K.Val) + covNew);
            newMu = covNew * sum_WtTX;

            obj.Z.updateDistributionsParameters(newMu, covNew);
        end

        function obj = qWUpdate(obj, it)
            for i = 1:obj.M
                obj.views(i).qWUpdate(it);
            end
        end

        function obj = qAlphaUpdate(obj)
            for i = 1:obj.M
                obj.views(i).qAlphaUpdate();
            end
        end

        function obj = qTauUpdate(obj)
            for i = 1:obj.M
                obj.views(i).qTauUpdate();
            end
        end



        %% fit() and ELBO
        function [elboVals, it, resArr] = fit(obj)
            elboVals = -Inf(1, obj.maxIter);
            resArr = cell(1, obj.maxIter);
        
            for it = 1:obj.maxIter
                obj.removeFactors(it);
                obj.qWUpdate(it);
                obj.qZUpdate();
                % if it > 0
                %     obj.updateRotation();
                % end  
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
                pAlpha = pAlpha + obj.views(i).alpha.E_LnP;
                pTau = pTau + obj.views(i).T.E_LnP;
                % q
                qW = qW + obj.views(i).W.H;
                qAlpha = qAlpha + obj.views(i).alpha.H;
                qTau = qTau + obj.views(i).T.H;
            end
            
            % Store to 'res'
            res.qZ = obj.Z.H; 
            res.qW = qW; 
            res.qAlpha = qAlpha; 
            res.qTau = qTau;
            % p
            res.pX = pX;
            res.pZ = obj.Z.E_LnP; 
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



        %% Additional methods
        function obj = removeFactors(obj, it, threshold)
            if nargin < 3
                threshold = Constants.LATENT_FACTORS_THRESHOLD;
            end
            % Calculate the average of the square of elements for each row of Z
            avgSquare = mean(obj.Z.E.^2, 2);
        
            removeIdx = find(avgSquare < threshold);

            if isempty(removeIdx)
                return;
            end

            % if obj.DEBUG
            disp(['Removed ', num2str(length(removeIdx)), ' factors in iteration ', num2str(it)]);
            % end
            
            % Update number of factors
            obj.K.Val = obj.K.Val - length(removeIdx);
        
            % Remove those rows from Z, corresponding columns from W, and elements from alpha
            obj.Z.removeDimensions(removeIdx);
            for m = 1:obj.M
                obj.views(m).alpha.removeDimensions(removeIdx);
                obj.views(m).W.removeDimensions(removeIdx);
            end
        end
    
        %% Getters
        function value = get.D(obj)
            value = zeros(obj.M, 1);
            for m = 1:length(obj.views)
                value(m) = obj.views(m).D;
            end
        end
    end
end