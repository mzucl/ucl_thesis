classdef GFAGroup < handle
    properties         
        X               % ViewHandler instance to hold the data

        W               % [D x K] GaussianDistributionContainer      
                        %       --- [size: D; for each row in W matrix]

                        % Prior over W is defined per columns (each column
                        % has its own precision parameter, but update
                        % equations are defined by rows, so we are
                        % representing W as a size D container in a row
                        % format.

        alpha           % [K x 1] GammaDistributionContainer         
                        %       --- [size: K]

        T               % [D x 1] GammaDistributionContainer         
                        %       --- [size: D]

        % TODO (high): Check if there is a better way to define these
        % inside this class - they are initialized in the constructor and
        % shouldn't be changed!
        Z
        K
    end

    % properties (Dependent, SetAccess = private)
    % end

    properties (Dependent)
        D
        N
    end

    methods
        %% Constructors
        function obj = GFAGroup(data, Z, K, featuresInCols)
            if nargin < 3
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
            elseif nargin < 4
                featuresInCols = true;
            end

            obj.X = ViewHandler(data, featuresInCols);
            obj.Z = Z;
            obj.K = K; % TODO(low): This can be a reference to avoid copying


            %% Model setup and initialization
            % alpha
            alphaPrior = GammaDistribution(Constants.DEFAULT_GAMMA_A, Constants.DEFAULT_GAMMA_B);
            obj.alpha = GammaDistributionContainer(repmat(alphaPrior, obj.K.Val, 1));

            % tau
            tauPrior = GammaDistribution(Constants.DEFAULT_GAMMA_A, Constants.DEFAULT_GAMMA_B);
            obj.T = GammaDistributionContainer(repmat(tauPrior, obj.D, 1));

            % W; sample from obj.alpha for the prior
            %       Should we do this? The values for alpha are so small!!!
            % wPrior = GaussianDistribution(0, diag(1./obj.alpha.Value));
            wPrior = GaussianDistribution(0, eye(obj.K.Val));
            obj.W = GaussianDistributionContainer(obj.D, wPrior, false);

            % Model initialization - second part
            % The first update equation is for W, so we need to initialize
            % everything that is used in those two equations and those
            % initilizations are given below.
            %   obj.T.expInit
            %   obj.alpha.expCInit
            % ----------------------------------------------------------------
            obj.T.setExpCInit(1000 * ones(obj.D, 1));        
            obj.alpha.setExpCInit(repmat(1e-1, obj.K.Val, 1));
        end



        %% Update methods
        % obj.alpha is GammaDistributionContainer
        function obj = qAlphaUpdate(obj)
            newAVal = obj.alpha.ds(1).prior.a + obj.D/2; % All 'a' values are the same
            newBVals = obj.alpha.ds(1).prior.b + 1/2 * obj.W.getExpectationOfColumnsNormSq();

            obj.alpha.updateAllDistributionsParams(newAVal, newBVals);
        end

        % obj.T is GammaDistributionContainer
        function obj = qTauUpdate(obj)
            newAVal = obj.T.ds(1).prior.a + obj.N/2; % All 'a' values are the same
            newBVals = obj.T.ds(1).prior.b + 1/2 * diag( ...
                obj.X.XXt ...
                - 2 * obj.W.EC * obj.Z.EC * obj.X.X' ...
                + obj.W.EC * obj.Z.E_CCt * obj.W.EC');

            obj.T.updateAllDistributionsParams(newAVal, newBVals);
        end
    
        function obj = qWUpdate(obj, it)
            % In the first iteration we perform the update based on the
            % initialized moments of T and alpha, and in every
            % subsequent iteration we use the 'normal' update equations

            % TODO (medium): DRY this code: define the vars for the values
            % that are different based on the value of 'it' before the
            % update equations. Use Utility.ternary(it == 1, ...)
            %   obj.alpha.expCInit
            %   obj.T.expInit
            disp(['Min W value: ', num2str(min(obj.W.EC, [], 'all'))]);
            disp(['Max W value: ', num2str(max(obj.W.EC, [], 'all'))]);
            if it > 1
                for d = 1:obj.D
                    covNew = Utility.matrixInverse(obj.T.E{d} * obj.Z.Tr_CtC * eye(obj.K.Val) + ...
                        diag(obj.alpha.EC));
                
                    muNew = covNew * obj.T.E{d} * obj.Z.EC * obj.X.getRow(d, true);
    
                    % TODO (low): Can be done in one call
                    obj.W.updateDistributionMu(d, muNew);
                    obj.W.updateDistributionCovariance(d, covNew);
                end
            else
                expInitT = obj.T.getExpCInit();
                expInitAlpha = obj.alpha.getExpCInit();

                for d = 1:obj.D
                    covNew = Utility.matrixInverse(expInitT(d) * obj.Z.Tr_CtC * eye(obj.K.Val) + ...
                        diag(expInitAlpha));
                
                    muNew = covNew * expInitT(d) * obj.Z.EC * obj.X.getRow(d, true);
    
                    obj.W.updateDistributionMu(d, muNew);
                    obj.W.updateDistributionCovariance(d, covNew);
                end
            end
        end

        function value = getExpectationLnW(obj)
            value = 0;
            colsNormSq = obj.W.getExpectationOfColumnsNormSq();
            for k = 1:obj.K.Val % TODO (high): This can be implemented as a dot product
                value = value + obj.alpha.E{k} * colsNormSq(k);
            end
            value = -1/2 * value + obj.D/2 * (obj.alpha.E_LnC - obj.K.Val * log(2*pi));
        end

        function value = getExpectationLnPX(obj)
            value = obj.N/2 * obj.T.E_LnC - obj.N * obj.D/2 * log(2 * pi) - 1/2 * ...
                obj.T.EC' * diag( ...
                obj.X.XXt ...
                - 2 * obj.W.EC * obj.Z.EC * obj.X.X' ...
                + obj.W.EC * obj.Z.E_CCt * obj.W.EC');
            % ans1 == ans2
            % ans1 = obj.T.EC' * diag( ...
            %     obj.X.XXt ...
            %     - 2 * obj.W.EC * obj.Z.EC * obj.X.X' ...
            %     + obj.W.EC * obj.Z.E_CCt * obj.W.EC');
            % ans2 = trace(obj.T.E_Diag * (obj.X.XXt ...
            %     - 2 * obj.W.EC * obj.Z.EC * obj.X.X' ...
            %     + obj.W.EC * obj.Z.E_CCt * obj.W.EC'));
        end
        

        %% Getters
        function value = get.D(obj)
            value = obj.X.D;
        end

        function value = get.N(obj)
            value = obj.X.N;
        end
    end
end