% TODO: View should be a base class and this should inherit from it, 
% DRY the code!

classdef BinaryView < BaseView
    properties         
        W               % [D x K] GaussianContainer      
                        %       --- [size: D; for each row in W matrix]

                        % Prior over W is defined per columns (each column
                        % has its own precision parameter, but update
                        % equations are defined by rows, so we are
                        % representing W as a size D container in a row
                        % format.

        alpha           % [K x 1] GammaContainer         
                        %       --- [size: K]

        mu              % [D x 1] Gaussian

        xi              % [D x N] Variational parameters

        bound           % Jaakkola or Bohning bound
    

    end
    

    
    methods
        %% Constructors
        function obj = BinaryView(data, Z, K, featuresInCols, bound)
            obj@BaseView(data, Z, K, false);
            % disp('BinaryView');
            % TODO: Deal with the default values in a proper way
            if nargin < 3
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
            elseif nargin < 4
                featuresInCols = true;
            end


            if nargin < 5 || bound == 'B'
                obj.bound = BohningBound(randn(obj.D, obj.N));
            elseif bound == 'J'
                obj.bound = JaakkolaBound(randn(obj.D, obj.N));
            end

            %% Model setup and initialization
            % type, size_, a, b, prior
            obj.alpha = GammaContainer( ...
                "SD", ...
                obj.K.Val, ...
                Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_A'), ...
                Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_B') ...
                );
            
            %                  dim, mu,    cov,  priorPrec
            obj.mu = Gaussian(obj.D, 0, eye(obj.D), 10^3);

            if isa(obj.bound, 'BohningBound')
                %                         type, size_, cols,   dim,        mu, cov, priorPrec
                obj.W = GaussianContainer("DS", obj.D, false, obj.K.Val, randn(obj.K.Val, obj.D));
            elseif isa(obj.bound, 'JaakkolaBound')
                obj.W = GaussianContainer("DD", obj.D, false, obj.K.Val, randn(obj.K.Val, obj.D));
            end

            % Model initialization - second part
            % The first update equation is for W, so we need to initialize
            % everything that is used in those two equations and those
            % initilizations are given below.
            %   obj.T.expInit
            %   obj.alpha.ExpInit
            % ----------------------------------------------------------------      
            obj.alpha.setExpInit(repmat(1e-1, obj.K.Val, 1));

            % Performed only once
            obj.qConstantUpdates();
        end





        %% Update methods
        function obj = qConstantUpdates(obj)
            % alpha.a
            obj.alpha.updateAllDistributionsA(obj.alpha.prior.a + obj.D/2);
        end

        function obj = qWUpdate(obj, it)
            alphaExp = Utility.ternary(it == 1, obj.alpha.getExpInit(true), obj.alpha.E_Diag);
            Rt = (obj.X.X + obj.bound.T + obj.bound.H.* obj.mu.E)';
            
            if isa(obj.bound, 'BohningBound')
                covNew = Utility.matrixInverse(1/4 * obj.Z.E_XXt + alphaExp);
                muNew = covNew * obj.Z.E * Rt;
            elseif isa(obj.bound, 'JaakkolaBound')
                Ht = obj.bound.H'; % Access dth col of Ht instead of dth row of H
                % Non-vectorized code
                covNew = repmat(eye(obj.K.Val), 1, 1, obj.D);
                muNew = zeros(obj.K.Val, obj.D);
                for d = 1:obj.D
                    covNew(:, :, d) = Utility.matrixInverse(obj.Z.getExpXtDX(diag(Ht(:, d)), false) + alphaExp); % XtDX = false
                end

                % Vectorized code
                V = reshape(obj.Z.E * Rt, obj.K.Val, 1, obj.D); % Columns of the matrix will be in the third dimension
                muNew = squeeze(pagemtimes(covNew, V));
            end
           
            obj.W.updateDistributionsParameters(muNew, covNew);
        end

        function obj = qAlphaUpdate(obj)
            bNew = obj.alpha.prior.b + 1/2 * obj.W.E_SNC;
            obj.alpha.updateAllDistributionsB(bNew);
        end




        function obj = qMuUpdate(obj)
            if isa(obj.bound, 'BohningBound')
                covNew = (1 / (obj.N / 4 + obj.mu.priorPrec)) * eye(obj.D);
                muNew = covNew * ((obj.X.X + obj.bound.T - 1/4 * obj.W.E * obj.Z.E) * ...
                    ones(obj.N, 1));

            elseif isa(obj.bound, 'JaakkolaBound')
                covNew = Utility.matrixInverse(diag(obj.bound.H * ones(obj.N, 1)) + ...
                    obj.mu.priorPrec * eye(obj.D));
                muNew = covNew * ((obj.X.X + obj.bound.T - obj.bound.H .* (obj.W.E * obj.Z.E)) * ...
                    ones(obj.N, 1));
            end
            
            obj.mu.updateParameters(muNew, covNew);
        end


        function obj = qXiUpdate(obj)
            if isa(obj.bound, 'BohningBound')
                xiNew = obj.W.E * obj.Z.E + obj.mu.E;
            elseif isa(obj.bound, 'JaakkolaBound')
                xiNew = 0; % Not implemented
            end

            obj.bound.updateXi(xiNew);
        end



        % Same as for continuous view -> move to the base class View
        function value = getExpectationLnW(obj)
            value = obj.W.E_SNC' * obj.alpha.E;
            value = -1/2 * value + obj.D/2 * (obj.alpha.E_LnX - obj.K.Val * log(2*pi));
        end

        function value = getExpectationLnPX(obj)
            % TODO: optimize trace in multiple places
            WZ_tr = (obj.W.E * obj.Z.E)';
            term1 = trace((WZ_tr + obj.mu.E_Xt) * (obj.X.X + obj.bound.T));
            if isa(obj.bound, 'BohningBound')
                term2 = 1/4 * trace(obj.W.E_XtX * obj.Z.E_XXt) + 1/2 * ones(1, obj.N) * WZ_tr * obj.mu.E ...
                    + obj.N/4 * obj.mu.E_XtX;

                const = sum(sum(-obj.bound.C + obj.bound.xi .* obj.bound.G ...
                    -1/8 * obj.bound.xi.^2));
                
            elseif isa(obj.bound, 'JaakkolaBound')
                term2 = 0; % Not implemented
                const = sum(sum(-obj.bound.C + obj.bound.xi .* obj.bound.G ...
                    -1/2 * obj.bound.xi.^2 .* obj.bound.H));
            end

            value = term1 - 1/2 * term2 + const;
        end
    end
end