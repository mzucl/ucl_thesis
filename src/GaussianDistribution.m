classdef GaussianDistribution < handle
    properties
        mu 
        cov
        dim
    end
    
    properties (Dependent)
        Expectation
        Variance
        Value
        ExpectationXt       % E[x^T]
        ExpectationXtX      % E[x^Tx]
        ExpectationXXt      % E[xx^T]
    end
  
    methods (Static, Access = private)
        function [mu, cov, dim] = initParameters(varargin)
            % One dimensional Gaussian
            switch nargin
                case 0
                    mu = Constants.DEFAULT_GAUSS_MU;
                    cov = Constants.DEFAULT_GAUSS_COV;
                    dim = Constants.DEFAULT_GAUSS_DIM;
    
                case 1
                    mu = varargin{1};

                    cov = eye(length(mu));
                    dim = length(mu);
    
                case 2
                    mu = varargin{1};
                    cov = varargin{2};

                    if Utility.isMatrix(cov) && ~Utility.isValidCovarianceMatrix(cov) 
                        error(['Error in class ' mfilename ': Parameter is not a valid covariance matrix.']);
                    end
    
                    if Utility.isArray(mu) && Utility.isArray(cov) && length(mu) ~= length(cov) || ...
                            Utility.isArray(mu) && Utility.isMatrix(cov) && length(mu) ~= size(cov, 1)
                        error(['Error in ' mfilename ': Dimensions do not match.']);
                    end
    
                    % 'mu' is an array
                    if Utility.isArray(mu)
                        dim = length(mu);
                        if ~Utility.isMatrix(cov)
                            if Utility.isArray(cov)
                                cov = diag(cov);
                            else
                                cov = cov * eye(length(mu));
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

                case 3
                    % This is the case where 'dim' is explicitly set, and it
                    % can only be used for the case where 'mu' and 'cov' are
                    % scalars; If we want to infer the 'dim' from other params,
                    % the previous constructor should be used.
                    mu = varargin{1};
                    cov = varargin{2};
                    dim = varargin{3};

                    if ~isscalar(mu) || ~isscalar(cov)
                        error(['Error in class ' mfilename ': Constructor with 3 parameters can be used only with scalar values for mu and cov.']);
                    end
    
                    % Spherical multivariate Gaussian distribution
                    mu = mu * ones(dim, 1);
                    cov = diag(cov * ones(dim, 1));
    
                otherwise
                    error(['Error in class ' mfilename ': Too many arguments passed into the constructor.']);
            end
        end

    end
    
    methods
        function obj = GaussianDistribution(varargin)
            [mu, cov, dim] = GaussianDistribution.initParameters(varargin{:});

            obj.mu = mu;
            obj.cov = cov;
            obj.dim = dim;
        end



        %% Methods
        % [NOTE] Gaussian parameters are correlated, don't allow setting
        % them independently (if we want to implement this properly we should make attributes private)! 
        % Also, 'dim' will never change in these updates, so for now they
        % can be set separately as long as the update doesn't change the
        % 'dim'.
        function updateParameters(obj, varargin)
            if nargin < 3
                error(['Error in class ' class(obj) ': Too few arguments passed.']);
            end
            
            [m, c, d] = GaussianDistribution.initParameters(varargin{:});

            if obj.dim ~= d
                error(['Error in class ' class(obj) ': Dimension cannot be changed via update method.']);
            end

            obj.mu = m;
            obj.cov = c;
            obj.dim = d;
        end

        function obj = updateCovariance(obj, cov)
            if nargin < 1
                error(['Error in class ' class(obj) ': Too few arguments passed.']);
            end

            if ~Utility.isValidCovarianceMatrix(cov)
                error(['Error in class ' class(obj) ': Argument is not a valid covariance matrix.']);
            end

            if obj.dim ~= size(cov, 1)
                error(['Error in class ' class(obj) ': Dimension cannot be changed via update method.']);
            end
            
            % Update covariance
            obj.cov = cov;
        end
        
        function obj = updateMu(obj, mu)
            if nargin < 1
                error(['Error in class ' class(obj) ': Too few arguments passed.']);
            end

            if obj.dim ~= length(mu)
                error(['Error in class ' class(obj) ': Dimension cannot be changed via update method.']);
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

        function value = get.Value(obj)
            value = mvnrnd(obj.mu, obj.cov)';
        end

        function value = get.ExpectationXtX(obj)
            value = obj.mu' * obj.mu + trace(obj.cov);
        end

        function value = get.ExpectationXXt(obj)
            value = obj.mu * obj.mu' + obj.cov;
        end
    end
end