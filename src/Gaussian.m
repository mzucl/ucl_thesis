classdef Gaussian < handle
    properties
        dim
        mu 
        cov
        priorPrec % scalar, only valid if the prior covariance is spherical (e.g. BPCA mu param) -> E_LnP
    end

    properties (Constant)
        VALIDATE = Constants.VALIDATE;
    end

    properties(Access = private)
        expInit
        cache = struct(...
            'E', NaN, ...
            'H', NaN, ...
            'E_Xt', NaN, ...
            'E_XtX', NaN, ...
            'E_XXt', NaN, ...
            'E_LnP', NaN);
    end
    
    properties (Dependent)
        E
        H                   % Entropy
        E_Xt                % E[x^T]
        E_XtX               % E[x^Tx]
        E_XXt               % E[xx^T]
        E_LnP               % E[ln(p(x))] wrt to the q(x)

        Var
        Val                 % Sample from the distribution
    end





    %% Private methods
    methods(Access = private)
        function isValid = validateDimIndices(obj, indices)
            isValid = true;
            for i = 1:length(indices)
                if indices(i) < 1 || indices(i) > obj.dim
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


   


    %% Static methods
    %% Options for the constructor Gaussian
    % ZERO PARAMETERS
    % -> default values for all parameters; dim = 1
    %
    % 1 PARAMETER: dim
    % OPTION 1: 'dim' is an instance of Gaussian
    %       -> set parameters to the values of that object
    % OPTION 2: 'dim' is a scalar
    %       -> standard normal Gaussian distribution; dim = dim;
    % 
    % 2 PARAMETERS: dim, mu
    % cov: eye(dim)
    % OPTION 1: 'mu' is a scalar
    %       -> all dimensions have the same mu
    % OPTION 2: 'mu' is a vector
    %       -> it must be a column vector and of length equal to 'dim'
    %
    % 3 PARAMETERS: dim, mu, cov
    % -> same as previous constructor regarding 'dim' and 'mu'
    % OPTION 1: 'cov' is a scalar -> spherical covariance
    % OPTION 2: 'cov' is an array -> diagonal covariance
    % OPTION 3: 'cov' is a matrix -> full covariance
    %
    % 4 PARAMETERS: dim, mu, cov, priorPrec
    % -> same as previous constructor + prior precision is set
    %
    %
    %
    methods (Static, Access = private)
        function [dim, mu, cov, priorPrec] = initParameters(varargin)

            % Default values
            dim = Constants.DEFAULT_GAUSS_DIM;
            mu = Constants.DEFAULT_GAUSS_MU;
            cov = 1/Constants.DEFAULT_GAUSS_PRECISION;
            priorPrec = Constants.DEFAULT_GAUSS_PRECISION;

            switch nargin
                case 1 % dim
                    % OPTION 1
                    if isscalar(varargin{1}) && isa(varargin{1}, 'Gaussian')
                        dist = varargin{1};
                        dim = dist.dim;
                        mu = dist.mu;
                        cov = dist.cov;
                        priorPrec = dist.priorPrec;

                    % OPTION 2
                    elseif Utility.isSingleNumber(varargin{1})
                        dim = varargin{1};
                        mu = zeros(dim, 1);
                        cov = eye(dim);
                    end

                case {2, 3, 4} % dim, mu
                    dim = varargin{1};

                    if Utility.isSingleNumber(varargin{2})
                        mu = repmat(varargin{2}, dim, 1);
                    else
                        if Gaussian.VALIDATE && size(varargin{2}, 1) ~= dim || size(varargin{2}, 2) ~= 1 % 'mu' is a column vector
                            error(['##### ERROR IN THE CLASS ' mfilename('class') ': Length of mu doesn''t match dimension.']);
                        end
                        mu = varargin{2};
                    end
                    cov = eye(dim);

                    if nargin > 2 % dim, mu, cov
                        covParam = varargin{3};

                        if Utility.isSingleNumber(covParam)
                            if Gaussian.VALIDATE && covParam <= 0
                                error(['##### ERROR IN THE CLASS ' mfilename ': Covariance parameter must be greater than 0.']);
                            end
                            cov = covParam * eye(dim); % Spherical
                            
                        elseif Utility.isArray(covParam)
                            if Gaussian.VALIDATE && (length(covParam) ~= dim || ~Utility.isValidCovarianceMatrix(diag(covParam)))
                                error(['##### ERROR IN THE CLASS ' mfilename ': Parameter is either not a valid covariance matrix or' ...
                                    ' dimensionality doesn''t match.']);
                            end
                            cov = diag(covParam); % Diagonal

                        elseif Utility.isMatrix(covParam)
                            if Gaussian.VALIDATE && (~isequal(size(covParam), [dim, dim]) || ~Utility.isValidCovarianceMatrix(covParam))
                                error(['##### ERROR IN THE CLASS ' mfilename ': Parameter is either not a valid covariance matrix or' ...
                                    'dimensionality doesn''t match.']);
                            end
                            cov = covParam; % Full
                        end

                        if nargin > 3 % dim, mu, cov, priorPrec
                            if Gaussian.VALIDATE && (~Utility.isSingleNumber(varargin{4}) || varargin{4} <= 0)
                                error(['##### ERROR IN THE CLASS ' mfilename ': Invalid precision parameter.']);
                            end
                            priorPrec = varargin{4};
                        end
                    end
            end % switch(nargin)
        end
    end
    




    methods
        %% Deep copy and operators overloading
        function newObj = copy(obj)
            newObj = Gaussian();
            
            newObj.mu = obj.mu;
            newObj.cov = obj.cov;
            newObj.dim = obj.dim;
            newObj.priorPrec = obj.priorPrec;
        end

        function isEqual = eq(obj1, obj2)
            if ~Utility.isNaNOrInstanceOf(obj1, 'Gaussian') || ...
                    ~Utility.isNaNOrInstanceOf(obj2, 'Gaussian')
                isEqual = false;
                return;
            end

            % Both are NaN
            if Utility.isNaN(obj1) && Utility.isNaN(obj2)
                isEqual = true;
                return;
            % Only one is NaN
            elseif xor(Utility.isNaN(obj1), Utility.isNaN(obj2))
                isEqual = false;
                return;
            end
                
            % Both are set - compare them!    
            isEqual =  obj1.dim == obj2.dim && isequal(obj1.mu, obj2.mu) && isequal(obj1.cov, obj2.cov) && ...
                obj1.priorPrec == obj2.priorPrec;
        end

        function isNotEqual = ne(obj1, obj2)
            isNotEqual = ~eq(obj1, obj2);
        end



        %% Constructor
        function obj = Gaussian(varargin)
            % Parameters: dim, mu, cov, priorPrec - they are all optional
            [dim, mu, cov, priorPrec] = Gaussian.initParameters(varargin{:});

            obj.dim = dim;
            obj.mu = mu;
            obj.cov = cov;
            obj.priorPrec = priorPrec;

            % Set initial expectation to the real expectation
            obj.setExpInit(obj.E);
        end



        

        %% Update methods
        function obj = updateMu(obj, mu)
            if obj.VALIDATE
                if nargin < 2
                    error(['##### ERROR IN THE THE CLASS ' class(obj) ': Too few arguments passed.']);
                end
    
                if size(mu, 1) ~= obj.dim % 'mu' is a column vector
                    error(['##### ERROR IN THE THE CLASS ' class(obj) ': Dimension cannot be changed via update method.']);
                end
            end
           
            obj.mu = mu;

            % Clear cache
            obj.clearCache();
        end

        function obj = updateCovariance(obj, cov)
            if obj.VALIDATE
                if nargin < 2
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
                end
    
                if ~Utility.isValidCovarianceMatrix(cov)
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Argument is not a valid covariance matrix.']);
                end
    
                if  size(cov, 1) ~= obj.dim
                    error(['##### ERROR IN THE CLASS ' class(obj) ': The dimensionality of new covariance doesn' ...
                        't match the distribution dimensionality.']);
                end
            end
            
            obj.cov = cov;

            % Clear cache
            obj.clearCache();
        end
        
        function removeDimensions(obj, indices)
            if nargin < 2 || isempty(indices)
                return; % No change
            end

            if obj.VALIDATE && ~obj.validateDimIndices(indices)
                error(['##### ERROR IN THE CLASS ' class(obj) ': Index out of range.']); 
            end
            
            obj.mu(indices) = [];
            obj.cov(:, indices) = [];
            obj.cov(indices, :) = [];
            obj.dim = obj.dim - length(indices);

            % Clear cache
            obj.clearCache();
        end


        


        %% Setters
        function obj = setExpInit(obj, value)
            obj.expInit = value;
        end





        %% Getters
        function value = getExpInit(obj)
            value = obj.expInit;
        end





        %% Dependent properties
        function value = get.E(obj)
            value = obj.mu;
        end

        function value = get.H(obj)
            if isnan(obj.cache.H)
                obj.cache.H = 1/2 * Utility.logDetUsingCholesky(obj.cov) + obj.dim/2 * (1 + log(2 * pi));
            end
            value = obj.cache.H;
        end

        function value = get.E_Xt(obj)
            if isnan(obj.cache.E_Xt)
                obj.cache.E_Xt = obj.mu';
            end
            value = obj.cache.E_Xt;
        end
       
        function value = get.E_XtX(obj)
            if isnan(obj.cache.E_XtX)
                obj.cache.E_XtX = dot(obj.mu, obj.mu) + trace(obj.cov);
            end
            value = obj.cache.E_XtX;
        end

        function value = get.E_XXt(obj)
            if isnan(obj.cache.E_XXt)
                obj.cache.E_XXt = obj.mu * obj.mu' + obj.cov;
            end
            value = obj.cache.E_XXt;
        end

        % [NOTE]: This formula is only valid if the prior covariance is spherical!
        % TODO (performance): use formula instead of obj.E_XtX vs use
        % cached value for it?
        function value = get.E_LnP(obj)
            if isnan(obj.cache.E_LnP)
                obj.cache.E_LnP = obj.dim/2 * log(obj.priorPrec / (2 * pi)) - ...
                    obj.priorPrec/2 * obj.E_XtX;
            end
            value = obj.cache.E_LnP;
        end

        function value = get.Var(obj)
            value = obj.cov;
        end

        function value = get.Val(obj)
            value = mvnrnd(obj.mu, obj.cov)'; % Transpose to get a column vector
        end
    end
end