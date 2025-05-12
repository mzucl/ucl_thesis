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

        W

        alpha

        doRotation
    end

    % TODO (high): Implement
    % methods(Access = private)
    %     function [M, N] = validateSources(obj, idx)
    % 
    %     end
    % end

    properties(Access = private, Constant)
        SETTINGS = ModelSettings.getInstance();
    end

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
            % Optional parameters: maxIter, tol
            if nargin < 2
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
            end

            % TODO (very high!): Implementent method below as soon as the
            % toy example with 2 views starts working
            % [obj.M, obj.N] = obj.validateSources(...)
            obj.M = length(data);
            obj.N = size(data{1}, 2); % Data passed in is DxN

            obj.K = DoubleWrapper(K);

            obj.maxIter = GFA.SETTINGS.DEFAULT_MAX_ITER;
            obj.tol = GFA.SETTINGS.DEFAULT_TOL;
            obj.doRotation = false;

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
            initZMu = randn(obj.K.Val, 1);

            %                         type, size_, cols, dim,     mu, cov, priorPrec
            obj.Z = GaussianContainer("DS", obj.N, true, obj.K.Val, zeros(obj.K.Val, 1)); % STEP1

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
                view = obj.views(m);
                WtT = view.W.E_Xt * view.T.E_Diag;
                covNew = covNew + view.W.getExpXtDX(view.T.E_Diag, true);
                sum_WtTX = sum_WtTX + WtT * view.X.X;
            end

            covNew = Utility.matrixInverse(eye(obj.K.Val) + covNew);
            muNew = covNew * sum_WtTX;

            obj.Z.updateDistributionsParameters(muNew, covNew);
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
        function [elboVals, it] = fit(obj, elboIterStep)
            if nargin < 2
                elboIterStep = 1;
            end

            elboVals = -Inf(1, obj.maxIter);
            % [NOTE] When elboIterStep ~= 1, indexing into elbo array is not done
            % using 'iter'; iter / elboIterStep + 1, but having independent
            % counter is cleaner; '+ 1' because we compute elbo in the
            % first iteration.
            elboIdx = 1;

            for it = 1:obj.maxIter

                obj.qZUpdate();
                obj.qWUpdate(it);
                % obj.qZUpdate();
                % if it > 0
                %     obj.updateRotation();
                % end  
                obj.qAlphaUpdate();
                obj.qTauUpdate();
                obj.removeFactors(it);

                if it ~= 1 && mod(it, elboIterStep) ~= 0
                    continue;
                end

                currElbo = obj.computeELBO();
                elboVals(elboIdx) = currElbo;

                if GFA.SETTINGS.DEBUG
                    if elboIdx ~= 1
                        disp(['======= ELBO increased by: ', num2str(currElbo - elboVals(elboIdx - 1))]);
                    end
                end

                % ELBO has to increase from iteration to iteration
                if elboIdx ~= 1 && currElbo < elboVals(elboIdx - 1)
                    fprintf(2, 'ELBO decreased in iteration %d by %f\n!!!', it, abs(currElbo - elboVals(elboIdx - 1)));
                end 

                % Check for convergence
                if elboIdx ~= 1 && abs(currElbo - elboVals(elboIdx - 1)) / abs(currElbo) < obj.tol
                    disp(['Convergence at iteration: ', num2str(it)]);
                    elboVals = elboVals(1:elboIdx); % cut the -Inf values at the end
                    break;
                end
                elboIdx = elboIdx + 1;

                if it == obj.maxIter
                    fprintf(2, 'Model did not converge in %d\n!!!', obj.maxIter);
                end
            end
        end

        function elbo = computeELBO(obj)
            elbo = 0;
            for i = 1:obj.M
                % p
                elbo = elbo + obj.views(i).getExpectationLnPX() + obj.views(i).getExpectationLnW() ...
                    + obj.views(i).alpha.E_LnP + obj.views(i).T.E_LnP + obj.views(i).W.H ...
                    + obj.views(i).alpha.H + obj.views(i).T.H;
            end

            elbo = elbo + obj.Z.H + obj.Z.E_LnP;
        end





        %% Additional methods
        function obj = removeFactors(obj, it, threshold)
            if nargin < 3
                threshold = GFA.SETTINGS.LATENT_FACTORS_THRESHOLD;
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

        % Big matrix W containing all views
        function value = get.W(obj)
            totalD = sum(obj.D);
            value = zeros(totalD, obj.K.Val);

            d = 0;
            for m = 1:obj.M
                Dm = obj.views(m).D;
                value(d + 1 : d + Dm, :) = obj.views(m).W.E;
                d = d + Dm;
            end
        end

        % Big matrix alpha, containing all alphas
        function value = get.alpha(obj)
            value = zeros(obj.K.Val, obj.M);
            for m = 1:obj.M
                value(:, m) = obj.views(m).alpha.E;
            end
        end
    end
end