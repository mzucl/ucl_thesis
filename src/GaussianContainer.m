classdef GaussianContainer < handle
    properties
        type   
        % "SS": shared 'mu', shared 'cov'
        % "SD": shared 'mu', different 'cov'        
        % "DS": different 'mu', shared 'cov'        -> IMPLEMENTED
        % "DD": different 'mu', different 'cov'     -> IMPLEMENTED

        cols                %   - 'true' when distributions describe columns of a matrix
                            %   - 'false' when they describe rows. This is important
                            % when the dependent properties are calculated

        dim                 % Dimensionality of each Gaussian component

        mu                  % Matrix that stores mu of each component in its COLUMNS (for both values of 'cols')

        cov                 % Matrix for the "DS" type
                            % Multidimensional array for the "DD" type
          
        priorPrec           % scalar or array, only valid if the prior covariance for each covariance is spherical
                            % -> E_LnP

        Size                % Number of distributions in the container
    end


    properties(Access = private)
        expInit
        cache = struct(...
            'E', NaN, ...
            'H', NaN, ...
            'E_Xt', NaN, ...
            'E_XXt', NaN, ...
            'E_XtX', NaN, ...
            'E_TrXtX', NaN, ...
            'E_LnP', NaN, ...
            'Cov_Tr', NaN, ...
            'E_SNC', NaN);

        cacheFlags = false(1, 9); % I invalidate all entries in the cache at the same time
                                  % but I don't store valid values at the
                                  % same time, so I need a flag for each
                                  % entry.
                                  % Hardcoded for optimization purposes!

        % cacheSize = 0; % Cached value for 'Size' is invalidated only when 'removeDimensions'
        %                % is called, so it make sense for it to have a separate cache! 
        %                % 0 is used because == 0 is much faster than isnan().
    end
    
    properties (Dependent)
        % Size                % Number of distributions in the container
        E                   % Matrix: expectation of the whole container
        H                   % Entropy of the collection
        E_Xt
        E_XXt               
        E_XtX   
        E_TrXtX             % Tr(E[X^TX]), e.g. Tr(Z^TZ)
        E_LnP              
        E_SNC               % E[squared norm of a column] for each column
        Cov_Tr              % Trace for the covariance
    end

    methods(Access = private)
        % Helper method that throws an error if index is not valid
        function validateIndex(obj, idx)
            if idx < 1 || idx > obj.Size 
                error(['##### ERROR IN THE CLASS ' class(obj) ': Index out of range.']); 
            end
        end

        function clearCache(obj)
            obj.cacheFlags = false(1, 9);
        end

        function isValid = validateDimIndices(obj, indices)
            isValid = true;
            for i = 1:length(indices)
                if indices(i) < 1 || indices(i) > obj.dim
                    isValid = false;
                    break;
                end
            end
        end

        % [NOTE] E[Z^T * Z] is a matrix where each entry is e.g. E[z1^T *
        % z2] (z1 and z2 are columns of the matrix Z and there
        % are independent random variable). Given the
        % independence this is E[z1^T]E[z2], but when i == j
        % then the independence doesn't hold, so for diagonal elements we
        % need to correct with trace(covariance).
        function value = getContainerExpectations(obj, XtX)
            if (obj.cols && XtX) || (~obj.cols && ~XtX)
                value = obj.mu' * obj.mu;
                if obj.type == "DS"
                    value = value + trace(obj.cov) * eye(obj.Size); % diag(repmat(trace(obj.cov), obj.Size, 1));
                elseif obj.type == "DD"
                    value = value + diag(obj.Cov_Tr); % returns trace for each covariance matrix
                end
            else
                value = obj.mu * obj.mu';
                if obj.type == "DS"
                    value = value + obj.Size * obj.cov;
                elseif obj.type == "DD"
                    value = value + sum(obj.cov, 3);
                end
            end
        end
    end

    


    methods (Static)

        function obj = loadobj(s)
            % Reconstruct the object from struct
            obj = GaussianContainer(s.type, s.Size, s.cols, s.dim, s.mu, s.cov, s.priorPrec);
        end
    end


    methods
        % TODO: Write tests for this and the loadObj
        function s = saveobj(obj)
            % Custom save logic: convert handle object to struct
            s.type = obj.type;
            s.cols = obj.cols;
            s.dim = obj.dim;
            s.mu = obj.mu;
            s.cov = obj.cov;
            s.priorPrec = obj.priorPrec;
            s.Size = obj.Size;
    
            % s.E = obj.E;
            % s.H = obj.H;
            % s.E_Xt = obj.E_Xt;
            % s.E_XXt = obj.E_XXt;
            % s.E_XtX = obj.E_XtX;
            % s.E_TrXtX = obj.E_TrXtX;
            % s.E_LnP = obj.E_LnP;
            % s.E_SNC = obj.E_SNC;
            % s.Cov_Tr = obj.Cov_Tr;

        end
        %% Options for the constructor GaussianContainer
        % 4 PARAMETERS
        % -> obj.mu is set to all zeros; covariance for all dimensions is
        % identity matrix
        %
        %
        % 5 PARAMETER: mu
        % OPTION 1: scalar
        %       -> all entries in obj.mu are set to that value
        % OPTION 2: array
        %       -> all columns in obj.mu are set to that value
        % OPTION 3: matrix
        %       -> obj.mu is set to that matrix
        %
        % [NOTE] cov is the same for all distributions!
        % 6 PARAMETERS: mu, cov
        % OPTION 1: scalar
        %       -> all covariance matrices are spherical
        % OPTION 2: array
        %       -> all covariance matrices are diagonal
        % OPTION 3: matrix
        %       -> all covariance matrices set to that value
        %
        % 7 PARAMETERS: mu, cov, priorPrec
        % -> same as previous constructor + prior precision is set
        %
        %
        function obj = GaussianContainer(type, size_, cols, dim, mu, cov, priorPrec)
            % Non-optional parameters: type, size_, cols, dim: 'size_' is the
            % size of the container, and 'dim' is the Gaussian dimension 
            if RunConfig.getInstance().validateInput
                if nargin < 4
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Too few argumentes passed in.']);
                end
                if type == "SD" || type == "SS"
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Not implemented yet.']);
                end
            end

            obj.type = type;
            obj.cols = cols;
            obj.dim = dim;
            obj.priorPrec = Utility.getConfigValue('Distribution', 'DEFAULT_GAUSS_PRECISION');

            % Preallocate + default values
            obj.mu = zeros(obj.dim, size_);
            if obj.type == "DS"
                obj.cov = eye(obj.dim);
            elseif obj.type == "DD"
                obj.cov = repmat(eye(obj.dim), 1, 1, size_);
            end
            
            if nargin > 4 % mu
                if Utility.isSingleNumber(mu)
                    obj.mu = mu * ones(obj.dim, size_);
                 
                elseif MatrixValidation.isNumericVector(mu)
                    if RunConfig.getInstance().validateInput && ~isequal(size(mu), [obj.dim, 1])
                        error(['##### ERROR IN THE CLASS ' mfilename ': Parameter mu dimensionality doesn''t match.']);
                    end
                    obj.mu = repmat(mu, 1, size_);

                elseif MatrixValidation.isNumeric2DMatrix(mu)
                    if RunConfig.getInstance().validateInput && ~isequal(size(mu), [obj.dim, size_])
                        error(['##### ERROR IN THE CLASS ' mfilename ': Parameter mu dimensionality doesn''t match.']);
                    end
                    obj.mu = mu;
                end
                if nargin > 5 % cov
                    % 'cov' is the same for all distributions
                    compCov = zeros(obj.dim, obj.dim); % Preallocate

                    if Utility.isSingleNumber(cov)
                        if RunConfig.getInstance().validateInput && cov <= 0
                            error(['##### ERROR IN THE CLASS ' mfilename ': Covariance parameter must be greater than 0.']);
                        end
                        compCov = cov * eye(dim); % Spherical
                        
                    elseif MatrixValidation.isNumericVector(cov)
                        if RunConfig.getInstance().validateInput && (length(cov) ~= dim || ~MatrixValidation.isCovarianceMatrix(diag(cov)))
                            error(['##### ERROR IN THE CLASS ' mfilename ': Parameter is either not a valid covariance matrix or' ...
                                ' dimensionality doesn''t match.']);
                        end
                        compCov = diag(cov); % Diagonal

                    elseif MatrixValidation.isNumeric2DMatrix(cov)
                        if RunConfig.getInstance().validateInput && (~isequal(size(cov), [dim, dim]) || ~MatrixValidation.isCovarianceMatrix(cov))
                            error(['##### ERROR IN THE CLASS ' mfilename ': Parameter is either not a valid covariance matrix or' ...
                                'dimensionality doesn''t match.']);
                        end
                        compCov = cov; % Full
                    end

                    if type == "DS"
                        obj.cov = compCov;
                    elseif type == "DD"
                        obj.cov = repmat(compCov, 1, 1, size_);
                    end

                    if nargin > 6 % priorPrec
                        if Utility.isSingleNumber(priorPrec)
                            if RunConfig.getInstance().validateInput && (priorPrec <= 0 || obj.type ~= "DS")
                                error(['##### ERROR IN THE CLASS ' mfilename ': Invalid precision parameter.']);
                            end
                            obj.priorPrec = priorPrec; % Spherical shared covariance
                            
                        elseif MatrixValidation.isNumericVector(priorPrec)
                            if RunConfig.getInstance().validateInput && (length(priorPrec) ~= size_ || any(priorPrec < 0) || obj.type ~= "DD")
                                error(['##### ERROR IN THE CLASS ' mfilename ': Parameter is either not a valid precision or' ...
                                    ' dimensionality doesn''t match.']);
                            end
                            obj.priorPrec = priorPrec;
                        end
                    end
                end
            end

            % Set obj.Size
            obj.Size = size(obj.mu, 2);
            
            % Set initial expectation to the actual expectation
            obj.setExpInit(obj.E);
        end
       




        %% Setters
        function obj = setExpInit(obj, value)
            if RunConfig.getInstance().validateInput
                if ~(obj.cols && isequal(size(value), [obj.dim, obj.Size]) || ...
                        ~obj.cols && isequal(size(value), [obj.Size, obj.dim]))
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Number of elements in the expectation must be equal to the number of ' ...
                        'components. ']);
                end
            end
            
            obj.expInit = value;
        end





        %% Getters
        % 'expInit' is a private property -> needs a getter for the access
        function value = getExpInit(obj)
            value = obj.expInit;
        end

        % Returns E[X^TDX] or E[XDX^T] where D is a diagonal matrix
        % TODO: Check if trace of diagonal is computed efficiently
        % TODO: A special case of this method is the method 'getContainerExpectations'
        % potentially merge them in code refactoring.
        % TODO: Add tests for this method
        function value = getExpXtDX(obj, diagMatrix, XtDX)
            if (obj.cols && XtDX) || (~obj.cols && ~XtDX)
                value = obj.mu' * diagMatrix * obj.mu;
                if obj.type == "DS"
                    value = value + trace(obj.cov * diagMatrix) * eye(obj.Size); 
                elseif obj.type == "DD"
                    % TODO: Haven't test this, not sure if it returns
                    % correct value
                    sigmaDiagTr = squeeze(sum((obj.cov * diagMatrix) .* eye(obj.dim), [1, 2])); % returns trace for each covariance matrix * diagonalMatrix
                    value = value + diag(sigmaDiagTr); 
                end
            else
                value = obj.mu * diagMatrix * obj.mu';
                if obj.type == "DS"
                    value = value + trace(diagMatrix) * obj.cov;
                elseif obj.type == "DD"
                    diagonal = reshape(diag(diagMatrix), 1, 1, []);
                    value = value + sum(obj.cov .* diagonal, 3);
                end
            end
        end




        
        %% Update methods
        % TODO (medium): Validate inputs!
        % Independent of the type (currently implemented)
        function obj = updateDistributionsMu(obj, mu)
            if RunConfig.getInstance().validateInput && nargin < 2
                CustomError.raiseError('InputCheck', CustomError.ERR_TOO_FEW_ARGUMENTS);
            end

            obj.mu = mu;

            % Clear cache
            obj.clearCache();
        end

        % For "DD" type 'cov' parameter must be a multidimensional array
        function obj = updateDistributionsCovariance(obj, cov)
            if RunConfig.getInstance().validateInput
                if nargin < 2
                    CustomError.raiseError('InputCheck', CustomError.ERR_TOO_FEW_ARGUMENTS);
                end
                if ~MatrixValidation.isCovarianceMatrix(cov)
                    CustomError.raiseError('InputCheck', ...
                    'Invalid update parameter. Ensure the sizes are correct and the covariance matrix is positive definite.');
                end
            end
           
            obj.cov = cov;

            % Clear cache
            obj.clearCache();
        end

        function obj = updateDistributionsParameters(obj, mu, cov)
            if RunConfig.getInstance().validateInput
                if nargin < 3
                    CustomError.raiseError('InputCheck', CustomError.ERR_TOO_FEW_ARGUMENTS);
                end
                if ~MatrixValidation.isCovarianceMatrix(cov)
                    CustomError.raiseError('InputCheck', ...
                    'Invalid update parameter. Ensure the sizes are correct and the covariance matrix is positive definite.');
                end
            end
            
            % Update params
            obj.mu = mu;
            obj.cov = cov;

            % Clear cache
            obj.clearCache();
        end

        function removeDimensions(obj, indices)
            if nargin < 2 || isempty(indices)
                return; % No change
            end

            if RunConfig.getInstance().validateInput && ~obj.validateDimIndices(indices)
                error(['##### ERROR IN THE CLASS ' class(obj) ': Index out of range.']); 
            end

            % Update dim
            obj.dim = obj.dim - length(indices);

            % Remove rows from mu
            obj.mu(indices, :) = [];
            
            if obj.type == "DS"
                % Remove rows and cols from the cov matrix
                obj.cov(:, indices) = [];
                obj.cov(indices, :) = [];
            elseif obj.type == "DD"
                newCov = zeros(obj.dim, obj.dim, obj.Size);

                for k = 1:obj.Size
                    newCov(:, :, k) = obj.cov(setdiff(1:end, indices), setdiff(1:end, indices), k);
                end
                obj.cov = newCov;
            end

            % Clear cache
            obj.clearCache();
            % obj.cacheSize = 0;

            obj.Size = size(obj.mu, 2);
        end




        
        %% Dependent properties
        % function value = get.Size(obj)
        %     if obj.cacheSize == 0
        %         obj.cacheSize = size(obj.mu, 2);
        %     end
        %     value = obj.cacheSize;
        % end        

        function value = get.E(obj)
            if ~obj.cacheFlags(1)
                obj.cache.E = LogicUtils.ternary(obj.cols, obj.mu, obj.mu'); % 'mu' is always stored in columns format
                obj.cacheFlags(1) = true;
            end
            value = obj.cache.E;
        end

        function value = get.E_Xt(obj)
            if ~obj.cacheFlags(2)
                obj.cache.E_Xt = LogicUtils.ternary(obj.cols, obj.mu', obj.mu);
                obj.cacheFlags(2) = true;
            end
            value = obj.cache.E_Xt;
        end

        function value = get.H(obj)
            if ~obj.cacheFlags(3)
                value = obj.Size * obj.dim/2 * (1 + log(2 * pi));
                if obj.type == "DS"
                    obj.cache.H = value + obj.Size/2 * LinearAlgebra.logDetCholesky(obj.cov);
                elseif obj.type == "DD"
                    for k = 1:obj.Size
                        value = value + 1/2 * LinearAlgebra.logDetCholesky(obj.cov(:, :, k));
                    end
                    obj.cache.H = value;
                end
                obj.cacheFlags(3) = true;
            end
            value = obj.cache.H;
        end

        function value = get.E_SNC(obj)
            if ~obj.cacheFlags(4)
                % cols = true
                if obj.cols == true
                    value = (sum(obj.mu.^2, 1))';
                    if obj.type == "DS"
                        obj.cache.E_SNC = value + trace(obj.cov);
                    elseif obj.type == "DD"
                        obj.cache.E_SNC = value + obj.Cov_Tr; % returns trace for each covariance matrix
                    end
    
                % cols = false
                else
                    value = sum(obj.mu.^2, 2);
                    if obj.type == "DS"
                        obj.cache.E_SNC = value + obj.Size * diag(obj.cov);
                    
                    % Sum diagonals of different cov matrices
                    elseif obj.type == "DD"
                        diagIdx = 1:size(obj.cov, 1) + 1:size(obj.cov, 1)^2;
                        covReshaped = reshape(obj.cov, size(obj.cov, 1) * size(obj.cov, 2), size(obj.cov, 3));
                        diags = covReshaped(diagIdx, :);
                        obj.cache.E_SNC = value + sum(diags, 2);
                    end
                end
                obj.cacheFlags(4) = true;
            end
            value = obj.cache.E_SNC;
        end

        function value = get.E_XXt(obj)
            if ~obj.cacheFlags(5)
                obj.cache.E_XXt = obj.getContainerExpectations(false);
                obj.cacheFlags(5) = true;
            end
            value = obj.cache.E_XXt;
        end

        function value = get.E_XtX(obj)
            if ~obj.cacheFlags(6)
                obj.cache.E_XtX = obj.getContainerExpectations(true);
                obj.cacheFlags(6) = true;
            end
            value = obj.cache.E_XtX;
        end

        function value = get.E_TrXtX(obj)
            if ~obj.cacheFlags(7)
                if obj.cols
                    value = sum(sum(obj.mu.^2, 2));
                    if obj.type == "DD"
                        obj.cache.E_TrXtX = value + sum(obj.Cov_Tr);
                    elseif obj.type == "DS"
                        obj.cache.E_TrXtX = value + obj.Size * trace(obj.cov);
                    end
                else % Optimization not possible
                    obj.cache.E_TrXtX = trace(obj.E_XtX);
                end
                obj.cacheFlags(7) = true;
            end
            value = obj.cache.E_TrXtX;
        end


        % [NOTE] Only implemented for the "DS" type that was needed for the
        % Z;
        function value = get.E_LnP(obj)
            if ~obj.cacheFlags(8)
                obj.cache.E_LnP = obj.Size * obj.dim / 2 * log (obj.priorPrec/(2 * pi)) ...
                    - 1/2 * obj.E_TrXtX;
                obj.cacheFlags(8) = true;
            end
            value = obj.cache.E_LnP;
        end

        function value = get.Cov_Tr(obj)
            if ~obj.cacheFlags(9)
                obj.cache.Cov_Tr = squeeze(sum(obj.cov .* eye(obj.dim), [1, 2]));
                obj.cacheFlags(9) = true;
            end
            value = obj.cache.Cov_Tr;    
        end

        function value = getCacheFlags(obj)
            value = obj.cacheFlags;
        end
    end
end