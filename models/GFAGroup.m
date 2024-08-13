classdef GFAGroup < handle
    properties         
        view            % ViewHandler instance to hold the data

        W               % [D x K] GaussianDistributionContainer      
                        %       --- [size: D; for each row in W matrix]

                        % Prior over W is defined per columns (each column
                        % has its own precision parameter, but update
                        % equations are defined by rows, so we are
                        % representing W as a size D container in a row
                        % format.

        alpha           % [K x 1] GammaDistributionContainer         
                        %       --- [size: K]

        tau             % [K x 1] GammaDistributionContainer         
                        %       --- [size: K]
    end

    properties (Dependent, SetAccess = private)
        K               % Initialized in the constructor and can't be changed
        Z               % Initialized in the constructor and can't be changed
                        % (not sure how this works with references in MATLAB, check this!)
    end

    
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

            obj.view = ViewHandler(data, featuresInCols);
            obj.K = K;
            obj.Z = Z; % TODO(low): This can be a reference to avoid copying


            %% Model setup and initialization
            % alpha
            alphaPrior = GammaDistribution(Constants.DEFAULT_GAMMA_A, Constants.DEFAULT_GAMMA_B);
            obj.alpha = GammaDistributionContainer(repmat(alphaPrior, K, 1));

            % tau
            tauPrior = GammaDistribution(Constants.DEFAULT_GAMMA_A, Constants.DEFAULT_GAMMA_B);
            obj.tau = GammaDistributionContainer(repmat(tauPrior, K, 1));
        end



        %% Update methods
        % obj.alpha is GammaDistributionContainer
        function obj = qAlphaUpdate(obj)
            newAVal = obj.alpha.distributions(1).prior.a + obj.D/2; % All 'a' values are the same
            newBVals = obj.alpha.distributions(1).prior.b + 1/2 * obj.W.getExpectationOfColumnsNormSq();

            obj.alpha.updateAllDistributionsParams(newAVal, newBVals);
        end

        % obj.tau is GammaDistributionContainer
        function obj = qTauUpdate(obj)
            newAVal = obj.tau.distributions(1).prior.a + obj.N/2; % All 'a' values are the same
            newBVals = 1/2 * diag( ...
                obj.view.XXt ...
                - 2 * obj.W.ExpectationC * obj.Z.ExpectationC * obj.view.X' ...
                + obj.W.ExpectationC * obj.Z.ExpectationCCt * obj.W.ExpectationC');

            obj.tau.updateAllDistributionsParams(newAVal, newBVals);
        end

        

        %% Getters
        function value = get.D(obj)
            value = obj.view.D;
        end

        function value = get.N(obj)
            value = obj.view.N;
        end

        function value = get.K(obj)
            value = obj.K;
        end
    end
end