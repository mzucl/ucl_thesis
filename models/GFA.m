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

            % [obj.M, obj.N] = obj.validateSources(...)
            obj.M = length(data);
            obj.N = 12;

            obj.K = K;

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
            % Initialize the model - set random values for the 'mu'
            % This means we will run the update equation for W first and
            % that we should set some values for all the moments that are
            % in those update equations.
            initZMu = randn(obj.K, 1);
            zPrior = GaussianDistribution(initZMu, eye(obj.K));
            obj.Z = GaussianDistributionContainer(obj.N, zPrior, true);

            % Initialize views
            obj.views = GFAGroup.empty(obj.M, 0);

            for i = 1:obj.M
                obj.views(i) = GFAGroup(data{i}, obj.Z, obj.K, true);
            end
        end



        %% Update methods
        % obj.Z is GaussianDistributionContainer(cols = true)
        function obj = qZUpdate(obj)
            % Update covariance
            % All latent variables have the same covariance
            covNew = zeros(obj.K);
            for i = 1:obj.M
                covNew = covNew + obj.views(i).W.E_Ct * ...
                    obj.views(i).T.E_Diag * obj.views(i).W.EC;
            end
            covNew = Utility.matrixInverse(eye(obj.K) + covNew);
            obj.Z.updateAllDistributionsCovariance(covNew);

            % Update mu
            for n = 1:obj.N
                newMu = zeros(obj.K, 1);
                for m = 1:obj.M
                    newMu = newMu + obj.views(i).W.E_Ct * ...
                        obj.views(i).T.E_Diag * obj.views(i).X.getObservation(n);
                end
                obj.Z.updateDistributionMu(n, muNew);
            end
        end

        function obj = qWUpdate(obj)
            for i = 1:obj.M
                obj.views(i).qWUpdate();
            end
        end

        function obj = qAlphaUpdate(obj)
            for i = 1:obj.M
                obj.views(i).qWUpdate();
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
                % obj.qWUpdate(it);
                obj.qZUpdate();
                % obj.qAlphaUpdate();
                % obj.qTauUpdate();

                % [currElbo, res] = obj.computeELBO();
                % 
                % resArr{it} = res;
                % 
                % % CHECK: ELBO has to increase from iteration to iteration
                % if it ~= 1 && currElbo < elboVals(it - 1)
                %     fprintf(2, 'ELBO decreased in iteration %d\n', it);
                % end 
                % 
                % elboVals(it) = currElbo;
                % 
                % if it ~= 1
                %     disp(['======= ELBO increased by: ', num2str(currElbo - elboVals(it - 1))]);
                % end
                % 
                % % Check for convergence
                % if it ~= 1 && abs(currElbo - elboVals(it - 1)) / abs(currElbo) < obj.tol
                %     disp(['Convergence at iteration: ', num2str(it)]);
                %     elboVals = elboVals(1:it); % cut the -Inf values at the end
                %     resArr = resArr(1:it);
                %     break;
                % end
            end
        end
    end
end