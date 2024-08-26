classdef GFAGroup < handle
    properties         
        X               % ViewHandler instance to hold the data

        W               % [D x K] GaussianContainer      
                        %       --- [size: D; for each row in W matrix]

                        % Prior over W is defined per columns (each column
                        % has its own precision parameter, but update
                        % equations are defined by rows, so we are
                        % representing W as a size D container in a row
                        % format.

        alpha           % [K x 1] GammaContainer         
                        %       --- [size: K]

        T               % [D x 1] GammaContainer         
                        %       --- [size: D]

        Z
        K
    end

    % TODO: ternaryOpt use!

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
            obj.K = K;
            

            %% Model setup and initialization
            %                          type, size_, a, b, prior
            obj.alpha = GammaContainer("SD", obj.K.Val, Constants.DEFAULT_GAMMA_A, Constants.DEFAULT_GAMMA_B);

            %                      type, size_, a, b, prior
            obj.T = GammaContainer("SD", obj.D, Constants.DEFAULT_GAMMA_A, Constants.DEFAULT_GAMMA_B);

            %                         type, size_, cols,   dim,   mu, cov, priorPrec
            obj.W = GaussianContainer("DD", obj.D, false, obj.K.Val);

            % Model initialization - second part
            % The first update equation is for W, so we need to initialize
            % everything that is used in those two equations and those
            % initilizations are given below.
            %   obj.T.expInit
            %   obj.alpha.expCInit
            % ----------------------------------------------------------------
            obj.T.setExpCInit(1000 * ones(obj.D, 1));        
            obj.alpha.setExpCInit(repmat(1e-1, obj.K.Val, 1));

            % Performed only once
            obj.qConstantUpdates();
        end



        %% Update methods
        function obj = qConstantUpdates(obj)
            % alpha.a
            obj.alpha.updateAllDistributionsA(obj.alpha.prior.a + obj.D/2);
            % T.a
            obj.T.updateAllDistributionsA(obj.T.prior.a + obj.N/2);
        end


        function obj = qAlphaUpdate(obj)
            newBVals = obj.alpha.prior.b + 1/2 * obj.W.E_SNC;
            obj.alpha.updateAllDistributionsB(newBVals);
        end


        function obj = qTauUpdate(obj)
            newBVals = obj.T.prior.b + 1/2 * diag( ...
                obj.X.XXt - 2 * obj.W.E * obj.Z.E * obj.X.X' ...
                + obj.W.E * obj.Z.E_XXt * obj.W.E');

            obj.T.updateAllDistributionsB(newBVals);
        end
    

        function obj = qWUpdate(obj, it)
            alphaExp = Utility.ternary(it == 1, obj.alpha.getExpCInit(), obj.alpha.E);
            TExp = Utility.ternary(it == 1, obj.T.getExpInit(), obj.T.E);

            for d = 1:obj.D
                covNew = Utility.matrixInverse(TExp{d} * obj.Z.E_XXt * eye(obj.K.Val) + ...
                    diag(alphaExp));
            
                muNew = covNew * TExp{d} * obj.Z.E * obj.X.getRow(d, true);

                obj.W.updateDistributionMu(d, muNew);
                obj.W.updateDistributionCovariance(d, covNew);
            end
        end

        function value = getExpectationLnW(obj)
            % The sum in the eq implemented as a dot product
            value = obj.W.getExpectationOfColumnsNormSq()' * obj.alpha.E;
            value = -1/2 * value + obj.D/2 * (obj.alpha.E_LnC - obj.K.Val * log(2*pi));
        end

        function value = getExpectationLnPX(obj)
            value = obj.N/2 * obj.T.E_LnC - obj.N * obj.D/2 * log(2 * pi) - 1/2 * ...
                obj.T.E' * diag( ...
                obj.X.XXt ...
                - 2 * obj.W.E * obj.Z.E * obj.X.X' ...
                + obj.W.E * obj.Z.E_XXt * obj.W.E');
            % ans1 == ans2
            % ans1 = obj.T.E' * diag( ...
            %     obj.X.XXt ...
            %     - 2 * obj.W.E * obj.Z.E * obj.X.X' ...
            %     + obj.W.E * obj.Z.E_XXt * obj.W.E');
            % ans2 = trace(obj.T.E_Diag * (obj.X.XXt ...
            %     - 2 * obj.W.E * obj.Z.E * obj.X.X' ...
            %     + obj.W.E * obj.Z.E_XXt * obj.W.E'));
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