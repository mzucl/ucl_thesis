%% [NOTE]
% Components inside the container are always independent and they
% appear as factors inside the product, e.g. p(Z) is product of p(zn) for
% each zn.
%
%
%%
classdef GammaContainer < handle
    properties
        % [NOTE]: Some types will affect obj.Size, check implementation if
        % we add more types
        type   
        % "SS": shared 'a', shared 'b'
        % "SD": shared 'a', different 'b'   -> IMPLEMENTED!
        % "DS": different 'a', shared 'b'   
        % "DD": different 'a', different 'b'

        % [NOTE]: For now only "SD" type is supported, thus 'a' is a
        % scalar, and 'b' is a vector
        a
        b

        % [NOTE]: They all have the same prior, thus we keep the
        % information here and not in the components distributions
        prior

        % It is dependent, but it changes infrequently
        Size                % Number of distributions in the container
    end

    properties(Access = private)
        expInit
        cache = struct(...
            'E', NaN, ...
            'E_Diag', NaN, ...
            'H', NaN, ...
            'E_LnP', NaN, ...
            'E_LnX', NaN);

        cacheFlags = false(1, 5); % Hardcoded for optimization purposes!

        % cacheSize = 0; % Cached value for 'Size' is invalidated only when 'removeDimensions'
        %                % is called, so it make sense for it to have a separate cache! 
        %                % 0 is used because == 0 is much faster than isnan().
    end

    
    properties (Dependent)   
        % Size                % Number of distributions in the container
        E                   % Array (column vector) of expectations of all components
        E_Diag              % Diagonal matrix of expectations of all components
        H                   % Entropy of the collection (sum of components entropies)
        E_LnP               % Sum of 'E_LnP' of all components
        E_LnX               % Sum of 'E_LnX' of all components

        Val           
    end



    methods(Access = private)
        function isValid = validateIndex(obj, idx)
            if idx < 1 || idx > obj.Size 
                error(['##### ERROR IN THE CLASS ' class(obj) ': Index out of range.']); 
            end
            isValid = true;
        end

        function isValid = validateIndices(obj, indices)
            isValid = true;
            for i = 1:length(indices)
                if ~obj.validateIndex(indices(i))
                    isValid = false;
                    break;
                end
            end
        end

        function clearCache(obj)
            obj.cacheFlags = false(1, 5);
        end
    end

    methods (Static)

        function obj = loadobj(s)
            % Reconstruct the object from struct
            obj = GammaContainer(s.type, s.Size, s.a, s.b);

            if isfield(s, 'priorClass') && strcmp(s.priorClass, 'Gamma')
                obj.prior = Gamma.loadobj(s.prior);  % recursively load
            end
        end
    end

    methods

        function s = saveobj(obj)
            % Custom save logic: convert handle object to struct
            s.type = obj.type;
            s.a = obj.a;
            s.b = obj.b;
            s.Size = obj.Size;

            if isa(obj.prior, 'Gamma')
                s.prior = saveobj(obj.prior);  % recursively save prior
                s.priorClass = 'Gamma';
            end
        end

        %% Options for the constructor GammaContainer
        % ZERO PARAMETERS
        % -> error!
        %
        % 1 PARAMETER: type
        % -> container has 1 Gamma with default parameters and default prior
        % 
        % 2 PARAMETERS: type, size_
        % -> container has 'size' Gamma with default parameters and default
        % prior
        %
        % 3 PARAMETERS: type, size_, a
        % -> container has 'size' Gamma with set 'a' and default 'b'
        %
        % 4 PARAMETERS: type, size_, a, b
        % -> container has 'size' Gamma with set 'a' and set 'b' if 'b' is
        % scalar; if b is an array of values then those values are used,
        % but the length of the array must match the size
        %
        % 5 PARAMETERS: type, size_, a, b, prior
        % -> same as previous constructor, but the prior is set
        %
        %%
        function obj = GammaContainer(type, size_, a, b, prior)
            if RunConfig.getInstance().validateInput
                if nargin == 0
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Too few argumentes passed in.']);
                end
                if type == "DS" || type == "SS" || type == "DD"
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Not implemented yet.']);
                end
            end

            obj.type = type;

            % Default values
            obj.a = Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_A');
            obj.b = Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_B'); % scalar
            obj.prior = Gamma();

            switch nargin
                case 2 % type, size
                    obj.b = repmat(Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_B'), size_, 1);

                case 3 % type, size, a
                    obj.a = a;
                    obj.b = repmat(Utility.getConfigValue('Distribution', 'DEFAULT_GAMMA_B'), size_, 1);
                   
                case {4, 5} % type, size, a, b
                    obj.a = a;
                    if Utility.isSingleNumber(b)
                        obj.b = repmat(b, size_, 1);
                    else
                        if RunConfig.getInstance().validateInput && size(b, 1) ~= size_ % 'b' must be a column vector
                            error(['##### ERROR IN THE CLASS ' class(obj) ': Length of b doesn''t match the size.']);
                        end
                        obj.b = b;
                    end
                    
                    if nargin > 4 % prior
                        if RunConfig.getInstance().validateInput && ~Utility.isNaNOrInstanceOf(prior, 'Gamma')
                            error(['##### ERROR IN THE CLASS ' class(obj) ': Invalid prior parameter.']);
                        end
                        obj.prior = prior.copy();
                    end
            end

            % Set obj.Size
            % It is a dependent property, but it only changes when
            % 'removeDimensions' is called
            obj.Size = length(obj.b);
            
            % Set initial expectation to the actual expectation
            obj.setExpInit(obj.E);
        end





        %% Update methods
        % [NOTE]: 'type' dependent
        function obj = updateAllDistributionsA(obj, a)
            if RunConfig.getInstance().validateInput
                if nargin < 2
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
                end
                if ~Utility.isSingleNumber(a)
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Parameter must be a numerical value.']);
                end
                if (a <= 0)
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Parameter a is a strictly positive number.']);
                end
            end

            obj.a = a;

            % Clear cache
            obj.clearCache();
        end

        % [NOTE]: 'type' dependent
        function obj = updateAllDistributionsB(obj, b)
            if RunConfig.getInstance().validateInput
                if nargin < 2
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
                end
                if ~(Utility.isSingleNumber(b) || Utility.isArray(b) && size(b, 1) == obj.Size) % 'b' must be a column vector
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Dimensions do not match.']);
                end
                if ~all(b > 0)
                    fprintf('b: %f\n', b);
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Parameter b is a strictly positive number.']);
                    
                end
            end

            if Utility.isSingleNumber(b)
                obj.b = repmat(b, obj.Size, 1);
            else
                obj.b = b;
            end

            % Clear cache
            obj.clearCache();
        end


        % Remove distributions from the container with idx in 'indices'
        function obj = removeDimensions(obj, indices)
            if nargin < 2 || isempty(indices)
                return; % No change
            end
            if RunConfig.getInstance().validateInput && ~obj.validateIndices(indices)
                error(['##### ERROR IN THE CLASS ' class(obj) ': Index out of range.']); 
            end
            obj.b(indices) = [];

            % Clear cache
            obj.clearCache();
            % obj.cacheSize = 0;

            obj.Size = length(obj.b);
        end





        %% Setters
        function obj = setExpInit(obj, value)
            if RunConfig.getInstance().validateInput
                if ~all(value > 0)
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Expectation is a strictly positive number.']);
                end
                if size(value, 1) ~= obj.Size % Expectation is a column vector
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Number of elements in the expectation must be equal to the number of ' ...
                        'components. ']);
                end
            end
            obj.expInit = value;
        end





        %% Getters
        function value = getExpInit(obj, d)
            if nargin < 2
                d = false;
            end
            value = Utility.ternary(d, diag(obj.expInit), obj.expInit);
        end




        
        %% Dependent properties
        % [NOTE]: 'type' dependent
        % function value = get.Size(obj)
        %     if obj.cacheSize == 0
        %         obj.cacheSize = length(obj.b);
        %     end
        %     value = obj.cacheSize;
        % end

        function value = get.E(obj)
            if ~obj.cacheFlags(1)
                obj.cache.E = obj.a ./ obj.b;
                obj.cacheFlags(1) = true;
            end
            value = obj.cache.E;
        end

        function value = get.E_Diag(obj)
            if ~obj.cacheFlags(2)
                obj.cache.E_Diag = diag(obj.E);
                obj.cacheFlags(2) = true;
            end
            value = obj.cache.E_Diag;
        end

        % [NOTE]: 'type' dependent
        function value = get.H(obj)
            if ~obj.cacheFlags(3)
                obj.cache.H = obj.Size * (gammaln(obj.a) - (obj.a - 1) * psi(obj.a) ...
                    + obj.a) - sum(log(obj.b));
                obj.cacheFlags(3) = true;
            end
            value = obj.cache.H;
        end

        % [NOTE]: 'type' dependent
        function value = get.E_LnP(obj)
            if ~obj.cacheFlags(4)
                if isa(obj.prior, 'Gamma') % The value is set only when prior is defined
                    obj.cache.E_LnP = obj.Size * (-gammaln(obj.prior.a) + obj.prior.a * log(obj.prior.b) + ...
                        (obj.prior.a - 1) * psi(obj.a)) - (obj.prior.a - 1) * sum(log(obj.b)) - ...
                        obj.prior.b * obj.a * sum(1 ./ obj.b);
                else 
                    obj.cache.E_LnP = NaN;
                end
                obj.cacheFlags(4) = true;
            end
            value = obj.cache.E_LnP;
        end

        function value = get.E_LnX(obj)
            if ~obj.cacheFlags(5)
                obj.cache.E_LnX = obj.Size * psi(obj.a) - sum(log(obj.b));
                obj.cacheFlags(5) = true;
            end
            value = obj.cache.E_LnX;
        end
        
        function value = get.Val(obj)
            value = zeros(obj.Size, 1);
            for i = 1:obj.Size
                value(i) = Gamma(obj.a, obj.b(i)).Val;
            end
        end
    end
end