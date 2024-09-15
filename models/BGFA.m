% TODO (medium): Create a base class for all models with optimization params
% and some other stuff
classdef BGFA < handle
    properties
        K               % Number of latent dimensions/principal components
    
        N               % Number of observations

        M               % Number of groups

        Mc              % Number of continous views

        Z               % [K x N] GaussianContainer [size: N; for each latent variable zn]
        
        views           % An array of GFAGroup instances

        bound           % Bound used for binary views

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
        function obj = BGFA(data, Mc, K, bound, maxIter, tol, doRotation) 
            % bound = 'J' or 'B'
            % Mc -> number of continous views

            % Optional parameters: maxIter, tol, doRotation
            if nargin < 4
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
            end

            % TODO (very high!): Implementent method below as soon as the
            % toy example with 2 views starts working
            % [obj.M, obj.N] = obj.validateSources(...)
            obj.Mc = Mc;
            obj.M = length(data);
            obj.N = size(data{1}, 2); % Data passed in is DxN

            obj.K = DoubleWrapper(K);

            obj.bound = bound;

            % Default values of optional parameters
            obj.maxIter = BGFA.SETTINGS.DEFAULT_MAX_ITER;
            obj.tol = BGFA.SETTINGS.DEFAULT_TOL;
            obj.doRotation = false;

            if nargin > 4
                obj.maxIter = maxIter;
                if nargin > 5
                    obj.tol = tol;
                    if nargin > 6
                        obj.doRotation = doRotation;
                    end
                end
            end

            %% Model setup and initialization
            if bound == 'B'
                %                         type, size_, cols, dim,     mu, cov, priorPrec
                obj.Z = GaussianContainer("DS", obj.N, true, obj.K.Val, zeros(obj.K.Val, 1)); % STEP1
            elseif bound == 'J'
                %                         type, size_, cols,   dim,   mu, cov, priorPrec
                obj.Z = GaussianContainer("DD", obj.N, true, obj.K.Val, randn(obj.K.Val, obj.D)); % STEP1
            end

            % Initialize views
            % TODO: If views had a base case this would be much easier
            % obj.views = SGFAGroup.empty(obj.M, 0);

            for i = 1:obj.Mc
                obj.views(i) = SGFAGroup(data{i}, obj.Z, obj.K, false); % featuresInCols = false;
            end

            for i = obj.Mc + 1:obj.M
                obj.views(i) = BinaryView(data{i}, obj.Z, obj.K, false, obj.bound); % featuresInCols = false;
            end

            obj.D = [obj.views.D];
        end



        %% Update methods
        function obj = qZUpdate(obj)
            if obj.bound == 'B'
                covNew = zeros(obj.K.Val);
                muNew = zeros(obj.K.Val, obj.N);
    
                for m = 1:obj.Mc
                    view = obj.views(m);
                    covNew = covNew + view.tau.E * view.W.E_XtX;
                    muNew = muNew + view.tau.E * view.W.E_Xt * (view.X.X - view.mu.E);
                end

                for m = obj.Mc:obj.M
                    view = obj.views(m);
                    covNew = covNew + 1/4 * view.W.E_XtX;
                    muNew = muNew + view.W.E_Xt * (view.X.X + view.bound.t() - 1/4 * view.mu.E);
                end
    
                covNew = Utility.matrixInverse(eye(obj.K.Val) + covNew);
                muNew = covNew * muNew;
    
                obj.Z.updateDistributionsParameters(muNew, covNew);
            elseif obj.bound == 'J'
            end
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

        function obj = qMuUpdate(obj)
            for i = 1:obj.M
                obj.views(i).qMuUpdate();
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
                obj.qMuUpdate();
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

                if SGFA.SETTINGS.DEBUG
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
            for m = 1:obj.M
                % p
                view = obj.views(m);
                elbo = elbo + view.getExpectationLnPX() + view.getExpectationLnW() ... % p(.)
                    + view.alpha.E_LnP + view.mu.E_LnP + view.tau.E_LnP ... % p(.)
                    + view.W.H + view.alpha.H + view.mu.H + view.tau.H; % q(.)
            end

            elbo = elbo + obj.Z.H + obj.Z.E_LnP;
        end





        %% Additional methods
        function obj = removeFactors(obj, it, threshold)
            if nargin < 3
                threshold = SGFA.SETTINGS.LATENT_FACTORS_THRESHOLD;
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