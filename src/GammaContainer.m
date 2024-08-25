%% [NOTE]
% Components inside the container are always independent and they
% appear as factors inside the product, e.g. p(Z) is product of p(zn) for
% each zn.
%
%
%%
classdef GammaContainer < handle
    properties
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
    end

    properties (Constant)
        VALIDATE = Constants.VALIDATE;
    end

    properties(Access = private)
        expCInit
        cache = struct(...
            'Size', NaN, ...
            'E', NaN, ...
            'E_Diag', NaN, ...
            'H', NaN, ...
            'E_LnP', NaN, ...
            'E_LnX', NaN);
    end
    
    properties (Dependent)   
        Size                % Number of distributions in the container
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
            fields = fieldnames(obj.cache);
            
            for i = 1:length(fields)
                obj.cache.(fields{i}) = NaN;
            end
        end
    end

    methods
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
            if obj.VALIDATE
                if nargin == 0
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Too few argumentes passed in.']);
                end
                if type == "DS" || type == "SS" || type == "DD"
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Not implemented yet.']);
                end
            end

            obj.type = type;

            % Default values
            obj.a = Constants.DEFAULT_GAMMA_A;
            obj.b = Constants.DEFAULT_GAMMA_B; % scalar
            obj.prior = Gamma();

            switch nargin
                case 2 % type, size
                    obj.b = repmat(Constants.DEFAULT_GAMMA_B, size_, 1);

                case 3 % type, size, a
                    obj.a = a;
                    obj.b = repmat(Constants.DEFAULT_GAMMA_B, size_, 1);
                   
                case {4, 5} % type, size, a, b
                    obj.a = a;
                    if Utility.isSingleNumber(b)
                        obj.b = repmat(b, size_, 1);
                    else
                        if obj.VALIDATE && size(b, 1) ~= size_ % 'b' must be a column vector
                            error(['##### ERROR IN THE CLASS ' class(obj) ': Length of b doesn''t match the size.']);
                        end
                        obj.b = b;
                    end
                    
                    if nargin > 4 % prior
                        if obj.VALIDATE && ~Utility.isNaNOrInstanceOf(prior, 'Gamma')
                            error(['##### ERROR IN THE CLASS ' class(obj) ': Invalid prior parameter.']);
                        end
                        obj.prior = prior.copy();
                    end
            end
            
            % Set initial expectation to the actual expectation
            obj.setExpCInit(obj.E);
        end





        %% Update methods
        % [NOTE]: 'type' dependent
        function obj = updateAllDistributionsA(obj, a)
            if obj.VALIDATE
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
            if obj.VALIDATE
                if nargin < 2
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
                end
                if ~(Utility.isSingleNumber(b) || Utility.isArray(b) && size(b, 1) == obj.Size) % 'b' must be a column vector
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Dimensions do not match.']);
                end
                if ~all(b > 0)
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
            if obj.VALIDATE && ~obj.validateIndices(indices)
                error(['##### ERROR IN THE CLASS ' class(obj) ': Index out of range.']); 
            end
            obj.b(indices) = [];

            % Clear cache
            obj.clearCache();
        end





        %% Setters
        function obj = setExpCInit(obj, value)
            if obj.VALIDATE
                if ~all(value > 0)
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Expectation is a strictly positive number.']);
                end
                if size(value, 1) ~= obj.Size % Expectation is a column vector
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Number of elements in the expectation must be equal to the number of ' ...
                        'components. ']);
                end
            end
            obj.expCInit = value;
        end





        %% Getters
        % 'expInit' is a private property -> needs a getter for the access
        function value = getExpCInit(obj)
            value = obj.expCInit;
        end




        
        %% Dependent properties
        % [NOTE]: 'type' dependent
        function value = get.Size(obj)
            if isnan(obj.cache.Size)
                obj.cache.Size = length(obj.b);
            end
            value = obj.cache.Size;
        end

        function value = get.E(obj)
            if isnan(obj.cache.E)
                obj.cache.E = obj.a ./ obj.b;
            end
            value = obj.cache.E;
        end

        function value = get.E_Diag(obj)
            if isnan(obj.cache.E_Diag)
                obj.cache.E_Diag = diag(obj.E);
            end
            value = obj.cache.E_Diag;
        end

        % [NOTE]: 'type' dependent
        function value = get.H(obj)
            if isnan(obj.cache.H)
                obj.cache.H = obj.Size * (gammaln(obj.a) - (obj.a - 1) * psi(obj.a) + obj.a) - sum(log(obj.b));
            end
            value = obj.cache.H;
        end

        % [NOTE]: 'type' dependent
        function value = get.E_LnP(obj)
            if isnan(obj.cache.E_LnP)
                if isa(obj.prior, 'Gamma') % The value is set only when prior is defined
                    obj.cache.E_LnP = obj.Size * (-gammaln(obj.prior.a) + obj.prior.a * log(obj.prior.b) + ...
                        (obj.prior.a - 1) * psi(obj.a)) - (obj.prior.a - 1) * sum(log(obj.b)) - ...
                        obj.prior.b * obj.a * sum(1 ./ obj.b);
                else 
                    obj.cache.E_LnP = NaN;
                end
            end
            value = obj.cache.E_LnP;
        end

        function value = get.E_LnX(obj)
            if isnan(obj.cache.E_LnX)
                obj.cache.E_LnX = obj.Size * psi(obj.a) - sum(log(obj.b));
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