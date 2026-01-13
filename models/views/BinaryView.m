classdef BinaryView < BaseView
    properties         
        xi              % [D x N] Variational parameters
        bound           % Jaakkola ('J') or Bohning ('B') bound
    end
    

    
    methods
        function obj = BinaryView(data, Z, K, featuresInCols, bound)
            CustomError.validateNumberOfParameters(nargin, 3, 5);
            
            % Set default values
            if nargin < 5, bound = 'B'; end
            if nargin < 4, featuresInCols = true; end

            obj@BaseView(data, Z, K, featuresInCols);

            %% Model setup and initialization
            if bound == 'B'
                obj.bound = BohningBound(randn(obj.D, obj.N));
            elseif bound == 'J'
                obj.bound = JaakkolaBound(randn(obj.D, obj.N));
            end
                
            if isa(obj.bound, 'JaakkolaBound')
                obj.W = GaussianContainer("DD", obj.D, false, obj.K.Val, randn(obj.K.Val, obj.D)); % Override!
            end


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
            alphaExp = LogicUtils.ternary(it == 1, obj.alpha.getExpInit(true), obj.alpha.E_Diag);
            Rt = (obj.X.X + obj.bound.T - obj.bound.H.* obj.mu.E)';
            
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