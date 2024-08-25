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
    end

    properties (Constant)
        VALIDATE = Constants.VALIDATE;
    end

    properties(Access = private)
        expCInit
        cache = struct(...
            'Size', NaN, ...
            'E', NaN, ...
            'H', NaN, ...
            'E_Xt', NaN, ...
            'E_XXt', NaN, ...
            'E_XtX', NaN, ...
            'E_TrXtX', NaN, ...
            'E_LnP', NaN, ...
            'E_SNC', NaN);
    end
    
    properties (Dependent)
        Size                % Number of distributions in the container
        E                   % Matrix: expectation of the whole container
        H                   % Entropy of the collection
        E_Xt
        E_XXt               
        E_XtX   
        E_TrXtX             % Tr(E[X^TX]), e.g. Tr(Z^TZ)
        E_LnP              
        E_SNC               % E[squared norm of a column] for each column
    end

    methods(Access = private)
        % Helper method that throws an error if index is not valid
        function validateIndex(obj, idx)
            if idx < 1 || idx > obj.Size 
                error(['##### ERROR IN THE CLASS ' class(obj) ': Index out of range.']); 
            end
        end

        function clearCache(obj)
            fields = fieldnames(obj.cache);
            
            for i = 1:length(fields)
                obj.cache.(fields{i}) = NaN;
            end
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

        % [!!!] Dependent on 'cols'
        function value = getContainerExpectations(obj, CtC)
            if (obj.cols && CtC) || (~obj.cols && ~CtC)
                value = zeros(obj.Size, obj.Size);
                % MATLAB stores data in column-major order
                for j = 1:obj.Size
                    for i = 1:obj.Size
                        % 
                        % [NOTE] E[Wt * W] is a matrix where each entry is e.g. E[w1^T *
                        % w2] (w1 and w2 are columns of the matrix W and there
                        % are independent random variable). Given the
                        % independence this is E[w1^T]E[w2], but when i == j
                        % then the independence doesn't hold and in that case
                        % we should use E_XtX from the GaussianDistribution class.
                        %
                        if i < j
                            value(i, j) = value(j, i); % Matrix is symmetric
                        elseif i ~= j
                            value(i, j) = obj.ds(i).E_Xt * ...
                            obj.ds(j).E;
                        else 
                            value(i, j) = obj.ds(i).E_XtX;
                        end
                    end
                end
            else
                dim = obj.ds(1).dim; % They are all of the same dimension
                value = zeros(dim, dim);
                for i = 1:obj.Size
                    value = value + obj.ds(i).E_XXt;
                end
            end
        end
    end

    



    methods
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
            if obj.VALIDATE
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
            obj.priorPrec = Constants.DEFAULT_GAUSS_PRECISION;

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
                 
                elseif Utility.isArray(mu)
                    if obj.VALIDATE && ~isequal(size(mu), [obj.dim, 1])
                        error(['##### ERROR IN THE CLASS ' mfilename ': Parameter mu dimensionality doesn''t match.']);
                    end
                    obj.mu = repmat(mu, 1, size_);

                elseif Utility.isMatrix(mu)
                    if Gaussian.VALIDATE && ~isequal(size(mu), [obj.dim, size_])
                        error(['##### ERROR IN THE CLASS ' mfilename ': Parameter mu dimensionality doesn''t match.']);
                    end
                    obj.mu = mu;
                end
                if nargin > 5 % cov
                    % 'cov' is the same for all distributions
                    compCov = zeros(obj.dim, obj.dim); % Preallocate

                    if Utility.isSingleNumber(cov)
                        if obj.VALIDATE && cov <= 0
                            error(['##### ERROR IN THE CLASS ' mfilename ': Covariance parameter must be greater than 0.']);
                        end
                        compCov = cov * eye(dim); % Spherical
                        
                    elseif Utility.isArray(cov)
                        if obj.VALIDATE && (length(cov) ~= dim || ~Utility.isValidCovarianceMatrix(diag(cov)))
                            error(['##### ERROR IN THE CLASS ' mfilename ': Parameter is either not a valid covariance matrix or' ...
                                ' dimensionality doesn''t match.']);
                        end
                        compCov = diag(cov); % Diagonal

                    elseif Utility.isMatrix(cov)
                        if obj.VALIDATE && (~isequal(size(cov), [dim, dim]) || ~Utility.isValidCovarianceMatrix(cov))
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
                            if obj.VALIDATE && (priorPrec <= 0 || obj.type ~= "DS")
                                error(['##### ERROR IN THE CLASS ' mfilename ': Invalid precision parameter.']);
                            end
                            obj.priorPrec = priorPrec; % Spherical shared covariance
                            
                        elseif Utility.isArray(priorPrec)
                            if obj.VALIDATE && (length(priorPrec) ~= size_ || any(priorPrec < 0) || obj.type ~= "DD")
                                error(['##### ERROR IN THE CLASS ' mfilename ': Parameter is either not a valid precision or' ...
                                    ' dimensionality doesn''t match.']);
                            end
                            obj.priorPrec = priorPrec;
                        end
                    end
                end
            end

            % Set initial expectation to the actual expectation
            obj.setExpCInit(obj.E);
        end
       




        %% Setters
        function obj = setExpCInit(obj, value)
            if obj.VALIDATE
                if ~(obj.cols && isequal(size(value), [obj.dim, obj.Size]) || ...
                        ~obj.cols && isequal(size(value), [obj.Size, obj.dim]))
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




        
        %% Update methods
        % Independent of the type (currently implemented)
        function obj = updateDistributionsMu(obj, mu)
            if obj.VALIDATE && nargin < 2
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
            end

            obj.mu = mu;
        end

        % For "DD" type 'cov' parameter must be a multidimensional array
        function obj = updateDistributionsCovariance(obj, cov)
            if obj.VALIDATE && nargin < 2
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
            end
            
            obj.cov = cov;
        end

        function removeDimensions(obj, indices)
            if nargin < 2 || isempty(indices)
                return; % No change
            end

            if obj.VALIDATE && ~obj.validateDimIndices(indices)
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
        end










        %% Dependent properties
        function value = get.Size(obj)
            value = size(obj.mu, 2);
        end        

        function value = get.E(obj)
            value = Utility.ternary(obj.cols, obj.mu, obj.mu'); % 'mu' is always stored in columns format
        end

        function value = get.E_Xt(obj)
            value = Utility.ternary(obj.cols, obj.mu', obj.mu);
        end

        function value = get.H(obj)
            value = obj.Size * obj.dim/2 * (1 + log(2 * pi));
            if obj.type == "DS"
                value = value + obj.Size/2 * Utility.logDetUsingCholesky(obj.cov);

            elseif obj.type == "DD"
                for k = 1:obj.Size
                    value = value + 1/2 * Utility.logDetUsingCholesky(obj.cov(:, :, k));
                end
            end
        end

        function value = get.E_SNC(obj)
            % cols = true
            if obj.cols == true
                value = (sum(obj.mu.^2, 1))';
                if obj.type == "DS"
                    value = value + trace(obj.cov);
                elseif obj.type == "DD"
                    value = value + squeeze(sum(obj.cov .* eye(obj.dim), [1, 2])); % returns trace for each covariance matrix
                end

            % cols = false
            else
                value = sum(obj.mu.^2, 2);
                if obj.type == "DS"
                    value = value + obj.Size * diag(obj.cov);
                
                % Sum diagonals of different cov matrices
                elseif obj.type == "DD"
                    diagIdx = 1:size(obj.cov, 1) + 1:size(obj.cov, 1)^2;
                    covReshaped = reshape(obj.cov, size(obj.cov, 1) * size(obj.cov, 2), size(obj.cov, 3));
                    diags = covReshaped(diagIdx, :);
                    value = value + sum(diags, 2);
                end
            end
        end


        % //// Tr
        % [NOTE] Only implemented for the "DS" type that was needed for the
        % Z;
        function value = get.E_LnP(obj)
            value = obj.Size * obj.dim / 2 * log (obj.priorPrec/(2 * pi)) ...
                - 1/2 * obj.E_TrXtX;
        end




        function value = get.E_XXt(obj)
            value = obj.getContainerExpectations(false); % CtC
        end

        function value = get.E_XtX(obj)
            value = obj.getContainerExpectations(true); % CtC
        end

        % TODO (high): This is implemented just for the cols format as it
        % was needed for the Tr(<Z^TZ>). Implement this properly and add
        % tests; The purpose of it is to optimize the calculation, where we
        % don't calculate all elements of the Z^TZ product, just diagonal
        % function value = get.Tr_CtC(obj)
        %     value = 0;
        %     for i = 1:obj.Size
        %         value = value + obj.ds(i).E_XtX;
        %     end
        % end




        % % TODO: Can the obj.Size change??? What about dim?
        % function value = get.E_LnPC(obj)
        %     value = (-obj.Size * obj.ds(1).dim/2) * log(2 * pi) - 1/2 * obj.Tr_CtC; 
        % end
    end
end