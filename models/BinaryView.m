% TODO: View should be a base class and this should inherit from it, 
% DRY the code!

classdef BinaryView < handle
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

        mu              % [D x 1] Gaussian

        xi              % [D x N] Variational parameters

        bound           % Jaakkola or Bohning bound
    

        Z
        K

        % CONSTANT (don't change after initialization) dependent properties
        D
        N
    end
    

    properties(Access = private, Constant)
        SETTINGS = ModelSettings.getInstance();
    end

    properties(Access = private)
        expInit
        cache = struct(...
            'C', NaN, ...
            'G', NaN, ...
            'H', NaN);

        cacheFlags = false(1, 3); % I invalidate all entries in the cache at the same time
                                  % but I don't store valid values at the
                                  % same time, so I need a flag for each
                                  % entry.
                                  % Hardcoded for optimization purposes!
    end

    properties (Dependent)
        % Matrices that depend on \xi
        % TODO: Check if I need them??? I can access them directly through
        % the bound
        C
        G
        H
    end


    methods
        %% Constructors
        function obj = BinaryView(data, Z, K, featuresInCols, bound)
            % disp('BinaryView');
            % TODO: Deal with the default values in a proper way
            if nargin < 3
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
            elseif nargin < 4
                featuresInCols = true;
            end

            obj.X = ViewHandler(data, featuresInCols);
            obj.Z = Z;
            obj.K = K;

            % Dependent properties
            obj.D = obj.X.D;
            obj.N = obj.X.N;

            if nargin < 5 || bound == 'B'
                obj.bound = BohningBound(randn(obj.D, obj.N));
            elseif bound == 'J'
                obj.bound = JaakkolaBound(randn(obj.D, obj.N));
            end

            %% Model setup and initialization
            %                          type, size_, a, b, prior
            obj.alpha = GammaContainer("SD", obj.K.Val, BinaryView.SETTINGS.DEFAULT_GAMMA_A, BinaryView.SETTINGS.DEFAULT_GAMMA_B);
            
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
            
            if isa(obj.bound, 'BohningBound')
                covNew = Utility.matrixInverse(1/4 * obj.Z.E_XXt + alphaExp);
                muNew = covNew * obj.Z.E * (obj.X.X + obj.bound.t() + obj.bound.h().* obj.mu.E);
            elseif isa(obj.bound, 'JaakkolaBound')
            end
           
            obj.W.updateDistributionsParameters(muNew, covNew);
        end

        function obj = qAlphaUpdate(obj)
            bNew = obj.alpha.prior.b + 1/2 * obj.W.E_SNC;
            obj.alpha.updateAllDistributionsB(bNew);
        end




        function obj = qMuUpdate(obj)
            if isa(obj.bound, 'BohningBound')
                covNew = (4 / obj.N + obj.mu.priorPrec) * eye(obj.D);
                muNew = covNew * ((obj.X.X + obj.bound.t() - 1/4 * obj.W.E * obj.Z.E) * ...
                    ones(obj.N, 1));

            elseif isa(obj.bound, 'JaakkolaBound')
                covNew = Utility.matrixInverse(diag(obj.bound.h() * ones(obj.N, 1)) + ...
                    obj.mu.priorPrec * eye(obj.D));
                muNew = covNew * ((obj.X.X + obj.bound.t() - obj.bound.h() .* (obj.W.E * obj.Z.E)) * ...
                    ones(obj.N, 1));
            end
            
            obj.mu.updateParameters(muNew, covNew);
        end


        function obj = xiUpdate(obj)
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
            term1 = trace((WZ_tr + obj.mu.E_Xt) * (obj.X.X + obj.bound.t()));
            if isa(obj.bound, 'BohningBound')
                term2 = 1/4 * trace(obj.W.E_XtX * obj.Z.E_XXt) + 1/2 * WZ_tr * obj.mu.E ...
                    + obj.N/4 * obj.mu.E_XtX;

                const = sum(sum(-obj.bound.c() + obj.bound.xi .* obj.bound.g() ...
                    -1/8 * obj.bound.xi.^2));
                
            elseif isa(obj.bound, 'JaakkolaBound')
                term2 = 0; % Not implemented
                const = sum(sum(-obj.bound.c() + obj.bound.xi .* obj.bound.g() ...
                    -1/2 * obj.bound.xi.^2 .* obj.bound.h()));
            end

            value = term1 - 1/2 * term2 + const;
        end










    end
end