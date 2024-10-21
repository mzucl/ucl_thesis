classdef (Abstract) Bound < handle
    properties
        xi
    end


    properties(Access = private)
        cache = struct(...
            'C', NaN, ...
            'G', NaN, ...
            'H', NaN, ...
            'T', NaN, ...
            'ElboConst', NaN);

        cacheFlags = false(1, 5);
    end


    properties (Dependent)   
        C
        G
        H
        T
        ElboConst
    end
    

    % These methods are defined as 'Static' because they are general utility 
    % functions hat are not tied to any specific instance of the 'Bound'
    % object, but still are related to the bounds.
    methods (Static)
        function value = sigma(x)
            value = exp(x)./(exp(x) + 1);
        end

        function value = lse(x)
            value = log(exp(x) + 1);
        end
    end


    methods(Access = private)
        function clearCache(obj)
            obj.cacheFlags = false(1, 5);
        end
    end


    methods(Abstract)
        h(obj)
    end


    methods
        function obj = Bound(xi)
            % Default constructor
            if nargin == 0
                obj.xi = 0;
            else
                obj.xi = xi;
            end
        end


        %% Update methods
        function obj = updateXi(obj, newXi)
            obj.xi = newXi;

            % Clear cache
            obj.clearCache();
        end


        %% Getters
        function value = c(obj)
            value = Bound.lse(obj.xi);
        end
    
        function value = g(obj)
            % TODO (medium): The second term is same as c, maybe this can 
            % be used for improving performance?
            value = exp(obj.xi - obj.C);
        end

        function value = t(obj)
            value = obj.h() .* obj.xi - obj.g();
        end

        function value = getElboConst(obj)
            value = sum(sum(-obj.c() + obj.xi .* obj.g() - 1/2 * obj.xi.^2 .* obj.h()));
        end

        % Dependent properties
        function value = get.C(obj)
            if ~obj.cacheFlags(1)
                obj.cache.C = Bound.lse(obj.xi);
                obj.cacheFlags(1) = true;
            end
            value = obj.cache.C;
        end

        function value = get.G(obj)
            if ~obj.cacheFlags(2)
                obj.cache.G = exp(obj.xi - obj.C);
                obj.cacheFlags(2) = true;
            end
            value = obj.cache.G;
        end

        function value = get.H(obj)
            if ~obj.cacheFlags(3)
                obj.cache.H = obj.h();
                obj.cacheFlags(3) = true;
            end
            value = obj.cache.H;
        end

        function value = get.T(obj)
            if ~obj.cacheFlags(4)
                obj.cache.T = obj.h() .* obj.xi - obj.g();
                obj.cacheFlags(4) = true;
            end
            value = obj.cache.T;
        end

        function value = get.ElboConst(obj)
            if ~obj.cacheFlags(5)
                obj.cache.ElboConst = sum(sum(-obj.C + obj.xi .* obj.G - 1/2 * obj.xi.^2 .* obj.H));
                obj.cacheFlags(5) = true;
            end
            value = obj.cache.ElboConst;
        end
    end
end