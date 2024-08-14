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

        T               % [K x 1] GammaDistributionContainer         
                        %       --- [size: K]

        % TODO (high): Check if there is a better way to define these
        % inside this class - they are initialized in the constructor and
        % shouldn't be changed!
        K
        Z
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
                return;
                % error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
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
            obj.T = GammaDistributionContainer(repmat(tauPrior, K, 1));
        end



        %% Update methods
        % obj.alpha is GammaDistributionContainer
        function obj = qAlphaUpdate(obj)
            newAVal = obj.alpha.ds(1).prior.a + obj.D/2; % All 'a' values are the same
            newBVals = obj.alpha.ds(1).prior.b + 1/2 * obj.W.getExpectationOfColumnsNormSq();

            obj.alpha.updateAllDistributionsParams(newAVal, newBVals);
        end

        % obj.tau is GammaDistributionContainer
        function obj = qTauUpdate(obj)
            newAVal = obj.T.ds(1).prior.a + obj.N/2; % All 'a' values are the same
            newBVals = obj.T.ds(1).prior.b + 1/2 * diag( ...
                obj.view.XXt ...
                - 2 * obj.W.EC * obj.Z.EC * obj.view.X' ...
                + obj.W.EC * obj.Z.E_CCt * obj.W.EC');

            obj.T.updateAllDistributionsParams(newAVal, newBVals);
        end

        

        %% Getters
        function value = get.D(obj)
            value = obj.view.D;
        end

        function value = get.N(obj)
            value = obj.view.N;
        end
    end
end