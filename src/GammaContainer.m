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
        % 'SS': shared 'a', shared 'b'
        % 'SD': shared 'a', different 'b'
        % 'DS': different 'a', shared 'b'
        % 'DD': different 'a', different 'b'

        % [NOTE]: For now only 'SD' type is supported, thus 'a' is a
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
            'Distributions', NaN, ...
            'Size', NaN, ...
            'E', NaN, ...
            'E_Diag', NaN, ...
            'H', NaN, ...
            'E_LnP', NaN);
    end
    
    properties (Dependent)   
        Size                % Number of distributions in the container
        E                   % Array (column vector) of expectations of all components
        E_Diag              % Diagonal matrix of expectations of all components
        H                   % Entropy of the collection (sum of components entropies)
        E_LnP               % Sum of 'E_LnP' of all components

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
    end

    methods
        %% Constructor
        function obj = GammaContainer(type, size, a, b, prior)
            if nargin == 0
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few argumentes passed in.']);
            end
            if type == "DS" || type == "SS" || type == "DD"
                error(['##### ERROR IN THE CLASS ' class(obj) ': Not implemented yet.']);
            end

            obj.type = type;

            % Default values
            obj.prior = Gamma();
            obj.a = Constants.DEFAULT_GAMMA_A;
            obj.b = Constants.DEFAULT_GAMMA_B; % scalar

            if nargin > 1 % 'size'
                obj.b = repmat(Constants.DEFAULT_GAMMA_B, size, 1);
                if nargin > 2 % 'a'
                    obj.a = a;
                    if nargin > 3 % 'b'
                        if obj.VALIDATE && length(b) ~= size
                            error(['##### ERROR IN THE CLASS ' class(obj) ': Length of b doesn''t match the size.']);
                        end
                        obj.b = b;
                        if nargin > 4 % 'prior'
                            if obj.VALIDATE && ~Utility.isNaNOrInstanceOf('Gamma')
                                error(['##### ERROR IN THE CLASS ' class(obj) ': Invalid argument for prior.']);
                            end
                            obj.prior = prior;
                        end
                    end
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
        end

        % [NOTE]: 'type' dependent
        function obj = updateAllDistributionsB(obj, b)
            if obj.VALIDATE
                if nargin < 2
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
                end
                if ~(isscalar(b) || Utility.isArray(b) && length(b) == obj.Size)
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Dimensions do not match.']);
                end
                if ~all(b > 0)
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Parameter b is a strictly positive number.']);
                end
            end

            obj.b = b;
        end

        % Remove distributions from the container with idx in 'indices'
        function obj = removeDistributions(obj, indices)
            if nargin < 2 || isempty(indices)
                return; % No change
            end
            if obj.VALIDATE && ~obj.validateIndices(indices)
                error(['##### ERROR IN THE CLASS ' class(obj) ': Index out of range.']); 
            end
            obj.b(indices) = [];
        end





        %% Setters
        function obj = setExpCInit(obj, value)
            if obj.VALIDATE
                if ~all(value > 0)
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Expectation is a strictly positive number.']);
                end
                if length(value) ~= obj.Size
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
            value = length(obj.b);
        end

        function value = get.E(obj)
            value = obj.a ./ obj.b;
        end

        function value = get.E_Diag(obj)
            value = diag(obj.E);
        end

        % [NOTE]: 'type' dependent
        function value = get.H(obj)
            value = obj.Size * (gammaln(obj.a) - (obj.a - 1) * psi(obj.a) + obj.a) - sum(log(obj.b));
        end
       
        % [NOTE]: 'type' dependent
        function value = get.E_LnP(obj)
            if Gamma.VALIDATE && ~isa(obj.prior, 'Gamma')
                error(['##### ERROR IN THE CLASS ' class(obj) ': Prior must be defined.']);
            end
            value = obj.Size * (-gammaln(obj.prior.a) + obj.prior.a * log(obj.prior.b) + ...
                (obj.prior.a - 1) * psi(obj.a)) - (obj.prior.a - 1) * sum(log(obj.b)) - ...
                obj.prior.b * obj.prior.a * sum(1 ./ obj.b);
        end
        
        function value = get.Val(obj)
            value = zeros(obj.Size, 1);
            for i = 1:obj.Size
                value(i) = Gamma(obj.a, obj.b(i)).Val;
            end
        end
    end
end