classdef GaussianDistribution < handle
    properties
        mu 
        cov
        dim
        prior
    end
    
    properties (Dependent)
        Expectation
        Variance
        PriorPrecision      % Set only if the prior covariance matrix is spherical or diagonal, otherwise NaN
                            %   if precision matrix of the full prior
                            %   covariance matrix is needed ->
                            %   (obj.prior.cov)^-1
        Value               % Sample from the distribution
        H                   % Entropy
        ExpectationXt       % E[x^T]
        ExpectationXtX      % E[x^Tx]
        ExpectationXXt      % E[xx^T]
    end

    %% Options for the constructor
    % SINGLE PARAMETER
    % ()                                -> default GaussianDistribution object
    %                                   (1-dim standard Gaussian distribution)
    % (numeric: scalar or array)        -> value used for the mu; cov: spherical
    % (GaussianDistribution)            -> GaussianDistribution object with the same prior

    % 2 PARAMETERS: 'mu' and 'cov'
    % 'dim' is inferred from 'mu' and 'cov' and they must agree!
    % ([], scalar)                      -> spherical covariance
    % ([], [])                          -> diagonal covariance
    % ([], matrix)                      -> full covariance
    % (scalar, matrix)                  -> full covariance
    % (scalar, [])                      -> diag. covariance
    % (scalar, scalar)                  -> 1-dim Gaussian distribution

    % 3 PARAMETERS: 'mu' and 'cov', 'prior'
    % Same as for the previous constructor, but with the addition of the
    % 'prior', which can be set to:
    %       - NaN
    %       - a GaussianDistribution with the same 'dim' parameter as the one
    %       inferred from the 'mu' and 'cov' params
            
    % 4 PARAMETERS
    % (mu, cov, prior, dim)             -> spherical multivariate (dim)
    %                                   GaussianDistribution object with 
    %                                   mu = mu * ones(dim, 1) and set prior


    
    %% Static methods
    methods (Static, Access = private)
        function [mu, cov, prior, dim] = getDefaultDistributionParams()
            mu = Constants.DEFAULT_GAUSS_MU;
            cov = 1 / Constants.DEFAULT_GAUSS_PRECISION;
            prior = NaN; % By default prior for the distribution is not defined
            dim = Constants.DEFAULT_GAUSS_DIM;
        end

        % [NOTE] % Left 'dim' as the last parameter because it can be
        % inferred from other parameters in some cases. 
        function [mu, cov, prior, dim] = initParameters(varargin)
            prior = NaN; % Default value for prior

            switch nargin
                % One dimensional Gaussian
                case 0
                    [mu, cov, prior, dim] = GaussianDistribution.getDefaultDistributionParams();

                case 1
                    mu = varargin{1};
                    % If the ONLY parameter is of a type GaussianDistribution that is
                    % the prior and we initialize the 'obj' and its 'prior'
                    % using that value
                    if isscalar(mu) && isa(mu, 'GaussianDistribution')
                        priorDist = mu; % Shallow copy is okay here!
                        mu = priorDist.mu;
                        cov = priorDist.cov;
                        dim = priorDist.dim;
                        prior = priorDist.copy();

                    elseif isnumeric(mu)
                        if size(mu, 1) == 1
                            mu = mu'; % 'mu' is a column vector
                        end 
    
                        cov = eye(length(mu));
                        dim = length(mu);
                    else
                        error(['##### ERROR IN THE CLASS ' class(obj) ': Invalid arguments passed.']);
                    end

                case {2, 3}
                    mu = varargin{1};
                    if size(mu, 1) == 1
                        mu = mu'; % 'mu' is a column vector
                    end

                    cov = varargin{2};

                    if Utility.isMatrix(cov) && ~Utility.isValidCovarianceMatrix(cov) 
                        error(['##### ERROR IN THE CLASS ' mfilename ': Parameter is not a valid covariance matrix.']);
                    end
    
                    % [NOTE]: For 'cov' it is okay to use 'length' because it is passed to the diag() in
                    % case it is an array!
                    if Utility.isArray(mu) && Utility.isArray(cov) && size(mu, 1) ~= length(cov) || ...
                            Utility.isArray(mu) && Utility.isMatrix(cov) && size(mu, 1) ~= size(cov, 1)
                        error(['##### ERROR IN THE CLASS ' mfilename ': Dimensions for mu and cov do not match.']);
                    end

                    % ------------------------------------------------------
                    % When we get this far mu and cov parameters are valid and all
                    % dimensions match
    
                    % 'mu' is an array
                    if Utility.isArray(mu)
                        dim = size(mu, 1);
                        % If it is a matrix it is already set with 'cov = varargin{2};' above
                        if ~Utility.isMatrix(cov)
                            if Utility.isArray(cov)
                                cov = diag(cov);
                            else
                                cov = cov * eye(size(mu, 1));
                            end
                        end

                    % 'mu' is a scalar        
                    else 
                        if Utility.isMatrix(cov)
                            mu  = mu * ones(size(cov, 1), 1);
                            dim = size(cov, 1);
                        elseif Utility.isArray(cov)
                            mu = mu * ones(length(cov), 1);
                            cov = diag(cov);
                            dim = length(cov);
                        else
                            % Both 'mu' and 'cov' are scalars
                            dim = 1;
                        end
                    end
                
                    % Validate 'prior'
                    % I validate it here and not above with mu and cov
                    % because the 'dim' depends on which case we are in, so
                    % we can validate the dimensionality of the prior when
                    % we get 'final' value for 'mu' and 'cov' and thus 'dim'
                    if nargin == 3
                        if ~Utility.isNaNOrInstanceOf(varargin{3}, 'GaussianDistribution')
                            error(['##### ERROR IN THE CLASS ' mfilename ': Prior must be an instance of a class or NaN.']);
                        elseif ~Utility.isNaN(varargin{3}) && varargin{3}.dim ~= size(mu, 1)
                            error(['##### ERROR IN THE CLASS ' mfilename ': Dimension of the prior must match with mu and cov.']);
                        end
                    end 

            
                    if nargin == 3 && ~Utility.isNaN(varargin{3}) % Set the prior only if not NaN is passed in for priorPrec
                        prior = varargin{3}.copy();
                    end

                case 4
                    % [NOTE] 'NaN' can be passed for 'prior'
                    % This is the case where 'dim' is explicitly set, and it
                    % can ONLY be used for the case where 'mu' and 'cov' are
                    % scalars. This means 'cov' is just sigma^2. 
                    % In case where we want to infer the 'dim' from other params,
                    % constructors with two or three parameters should be used.
                    mu = varargin{1};
                    cov = varargin{2};
                    prior = varargin{3};
                    dim = varargin{4};

                    if size(mu, 1) == 1
                        mu = mu'; % 'mu' is a column vector
                    end 
                    
                    if ~Utility.isSingleNumber(mu) || ~Utility.isSingleNumber(cov) || ~isscalar(prior) || ...
                            ~Utility.isNaNOrInstanceOf(prior, 'GaussianDistribution') || ~Utility.isSingleNumber(dim)
                        error(['##### ERROR IN THE CLASS ' mfilename ': Method with four parameters can be used only with scalar values for mu, cov, and prior.']);
                    end
    
                    if ~Utility.isNaN(prior) && prior.dim ~= dim
                        error(['##### ERROR IN THE CLASS ' mfilename ': dim parameter and dimension of the prior do not align']);
                    end

                    % Spherical multivariate Gaussian distribution
                    mu = mu * ones(dim, 1); % Column vector
                    cov = diag(cov * ones(dim, 1));
    
                otherwise
                    error(['##### ERROR IN THE CLASS ' mfilename ': Too many arguments passed.']);
            end
        end

    end
    
    methods
        %% Deep copy and operators overloading
        function newObj = copy(obj)
            newObj = GaussianDistribution();
            
            newObj.mu = obj.mu;
            newObj.cov = obj.cov;
            newObj.dim = obj.dim;
            
            % Copy the prior (manually)
            if ~Utility.isNaN(obj.prior)
                newObj.prior = GaussianDistribution();

                newObj.prior.mu = obj.prior.mu;
                newObj.prior.cov = obj.prior.cov;
                newObj.prior.dim = obj.prior.dim;
            end
        end

        function isEqual = eq(obj1, obj2)
            if ~Utility.isNaNOrInstanceOf(obj1, 'GaussianDistribution') || ...
                    ~Utility.isNaNOrInstanceOf(obj2, 'GaussianDistribution')
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
            isEqual = all(obj1.mu == obj2.mu, 'all') && ...
                all(obj1.cov == obj2.cov, 'all') && obj1.dim == obj2.dim;

            % If parameters mu, cov and dim are not equal, they are different objects
            % regardless of the prior!
            %   Second part of the condition is added because obj1.prior ==
            %   obj2.prior won't call 'isEqual' if both are NaN, and we are
            %   relaying on that call for comparing priors.
            if ~isEqual || Utility.isNaN(obj1.prior) && Utility.isNaN(obj2.prior)
                return;
            end

            isEqual = obj1.prior == obj2.prior;
        end

        function isNotEqual = ne(obj1, obj2)
            isNotEqual = ~eq(obj1, obj2);
        end



        %% Constructor
        function obj = GaussianDistribution(varargin)
            % [NOTE] We pass in priorPrec not the prior (scalar not the
            % GaussianDistribution object)
            % Parameters: mu, cov, priorPrec, dim - they are all optional
            [mu, cov, prior, dim] = GaussianDistribution.initParameters(varargin{:});

            obj.mu = mu;
            obj.cov = cov;
            obj.prior = prior;
            obj.dim = dim;
        end



        %% Update methods
        % [NOTE] Gaussian parameters are correlated (via 'dim'), don't allow setting
        % them independently (if we want to implement this properly we should make attributes private)! 
        % Also, 'dim' will never change in these updates, so for now they
        % can be set separately as long as the update doesn't change the
        % 'dim'. 
        % Dimension and prior cannot be changed in these update
        % methods.

        % Updates both 'mu' and 'cov'
        function updateParameters(obj, varargin)
            if nargin ~= 3
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
            end
            
            % We can ignore 'prior' returned, but not 'dim', because 'dim'
            % can be inferred from new 'mu' and 'cov' params, thus can be
            % changed, unlike the 'prior'
            [m, c, ~, d] = GaussianDistribution.initParameters(varargin{:});

            if obj.dim ~= d
                error(['##### ERROR IN THE CLASS ' class(obj) ': Dimension cannot be changed via update method.']);
            end

            if size(m, 1) == 1
                m = m'; % Column vector
            end
            
            % Update parameters
            obj.mu = m;
            obj.cov = c;
        end

        function obj = updateCovariance(obj, cov)
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
            
            % Update covariance
            obj.cov = cov;
        end
        
        function obj = updateMu(obj, mu)
            if nargin < 2
                error(['##### ERROR IN THE THE CLASS ' class(obj) ': Too few arguments passed.']);
            end

            if length(mu) ~= obj.dim
                error(['##### ERROR IN THE THE CLASS ' class(obj) ': Dimension cannot be changed via update method.']);
            end
            
            % 'mu' is a column vector
            if size(mu, 1) == 1
                mu = mu';
            end

            obj.mu = mu;
        end



        %% Getters
        function value = get.Expectation(obj)
            value = obj.mu;
        end

        function value = get.ExpectationXt(obj)
            value = obj.mu';
        end
        
        function value = get.Variance(obj)
            value = obj.cov;
        end

        function value = get.H(obj)
            value = 1/2 * log(det(obj.cov)) + obj.dim/2 * (1 + log(2 * pi)); 
        end
        
        function value = get.PriorPrecision(obj)
            if Utility.isNaN(obj.prior)
                value = NaN;
            else
                [isDiag, diagEl] = Utility.checkAndExtractDiagonal(obj.prior.cov);
                value = Utility.ternary(isDiag, 1./diagEl, NaN);
            end
        end

        function value = get.Value(obj)
            value = mvnrnd(obj.mu, obj.cov)'; % Transpose to get a column vector
        end

        function value = get.ExpectationXtX(obj)
            value = obj.mu' * obj.mu + trace(obj.cov);
        end

        function value = get.ExpectationXXt(obj)
            value = obj.mu * obj.mu' + obj.cov;
        end
    end
end