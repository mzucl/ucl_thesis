classdef GaussianDistributionContainer < handle
    properties
        distributions       % A list of distributions inside the container
        cols                % Boolean: true when the distributions are columns of a matrix
                            % false when they are rows. This is important
                            % when the dependent properties are calculated
    end
    
    properties (Dependent)
        Size                % Number of distributions in the container
        Expectation         % Under the assumption of independence of each component; 
                            % This is a cell array where each entry is an
                            % expectation of one component;
        ExpectationXXt      % Similar to above, each entry is E[XXt]
        ExpectationXtX      % Similar to above, each entry is E[XtX], which is equivalent to E[|X|^2]
        ExpectationC        % 'C' stands for container, and this is expectation of the whole container
        ExpectationCt 
        ExpectationCtC
        H                   % Similar to above, each entry is entropy of a single component
        HC                  % Entropy of the collection
    end

    methods(Access = private)
        % Helper method that throws an error if index is not valid
        function validateIndex(obj, idx)
            if idx < 1 || idx > obj.Size 
                error(['Error in ' class(obj) ': Index out of range.']); 
            end
        end
    end

    methods
        function obj = GaussianDistributionContainer(dim, cols, numDistributions, prec)
            % [NOTE] All distibutions in the container have the same 'dim' (for now).
            % Constructor can be extended to handle other cases beside the
            % two described below.
            switch nargin
                case {0, 1, 2}
                    error(['Error in class ' class(obj) ': Too few arguments passed.']);
                
                case 3
                    % Constructor with 3 scalar parameters
                    %   All distributions are standard normal multivariate
                    %   Gaussians with the same dimensionality
                    obj.cols = cols;
                    obj.distributions = repmat(GaussianDistribution(), numDistributions, 1); % Preallocate
                    for i = 1:numDistributions
                        obj.distributions(i) = GaussianDistribution(0, 1, dim);
                    end

                case 4
                    % Constructor with 4 scalar parameters
                    %   Option 1: When 'prec' is an array, it is an array
                    %   of precisions for each of the components in container.
                    %   All distributions are zero-mean normal multivariate
                    %   Gaussians with the same dimensionality and
                    %   different covariances given by precision array
                    %   
                    %   Option 2: When 'prec' is a matrix it is a
                    %   covariance for each of the components
                    obj.cols = cols;
                    obj.distributions = repmat(GaussianDistribution(), numDistributions, 1); % Preallocate

                    if Utility.isArray(prec)
                        if length(prec) ~= numDistributions
                            error(['Error in class ' class(obj) ': Invalid arguments.']);
                        end
    
                        for i = 1:numDistributions
                            obj.distributions(i) = GaussianDistribution(0, 1./prec(i), dim);
                        end
                    elseif Utility.isMatrix(prec)
                        if size(prec, 1) ~= dim || size(prec, 2) ~= dim
                            error(['Error in class ' class(obj) ': Invalid arguments.']);
                        end

                        if ~Utility.isValidCovarianceMatrix(prec)
                            error(['Error in class ' class(obj) ': Invalid arguments: covariance parameter is not a valid covariance matrix.']);
                        end

                        for i = 1:numDistributions
                            obj.distributions(i) = GaussianDistribution(0, prec);
                        end
                            
                    else
                        error(['Error in class ' class(obj) ': Invalid arguments.']);
                    end
            end
        end
       


        %% Methods
        function dist = getDistribution(obj, idx)
            % Returns the distribution at index 'idx'
            obj.validateIndex(idx);

            dist = obj.distributions(idx);
        end

        % [NOTE] This is used in the scenario when W is stored in the rows
        % format, but we need expectation of the squared norm of a column
        function res = getExpectationOfColumnNormSq(obj)
            numOfCols = obj.distributions(1).dim;
            res = zeros(1, numOfCols);
            for idx = 1:numOfCols
                for i = 1:obj.Size
                    res(idx) = res(idx) + obj.distributions(i).mu(idx) ^ 2 + obj.distributions(i).cov(idx, idx) ^ 2;
                end
            end
        end

        function obj = updateDistribution(obj, idx, dist)
            obj.validateIndex(idx);

            % Updates the distribution at index 'idx'
            obj.distributions(idx) = dist;
        end

        function obj = updateDistributionParams(obj, idx, mu, cov)
            if nargin < 4
                error(['Error in class ' class(obj) ': Too few arguments passed.']);
            end
            obj.validateIndex(idx);

            obj.distributions(idx).updateParameters(mu, cov);
        end

        function obj = updateDistributionMu(obj, idx, mu)
            if nargin < 3
                error(['Error in class ' class(obj) ': Too few arguments passed.']);
            end
            obj.validateIndex(idx);

            obj.distributions(idx).updateMu(mu); % Will raise an exception if we try to change dimension
        end

        % [NOTE] The type of these update methods depends on the current
        % need (e.g. in the update equations for BPCA all covariances are set to the
        % same value). In future more general methods can be added (maybe
        % needed in the update equations of more complex models like GFA).
        function obj = updateAllDistributionsCovariance(obj, cov)
            % Update covariance of all distributions
            if nargin < 2
                error(['Error in class ' class(obj) ': Too few arguments passed.']);
            end
            
            for i=1:obj.Size
                obj.distributions(i).updateCovariance(cov);
            end
        end

        %% Getters
        function value = get.Size(obj)
            value = length(obj.distributions);
        end

        % [NOTE] The value for 'cols' is irrelevant for the properties that
        % are of type cell, because they are used with indexing to access
        % components values. And for each component (multivariate Gauss) we
        % use columns for mean and expectations, this 'cols' parameter is
        % not related to that.
        function value = get.Expectation(obj)
            value = cell(1, obj.Size);
            for i = 1:obj.Size
                value{i} = obj.distributions(i).Expectation;
            end
        end

        function value = get.ExpectationXXt(obj)
            value = cell(1, obj.Size);
            for i = 1:obj.Size
                value{i} = obj.distributions(i).ExpectationXXt;
            end
        end
        
        function value = get.ExpectationXtX(obj)
            value = cell(1, obj.Size);
            for i = 1:obj.Size
                value{i} = obj.distributions(i).ExpectationXtX;
            end
        end

        function value = get.ExpectationC(obj)
            value = Utility.ternary(obj.cols, cell2mat(obj.Expectation), cell2mat(obj.Expectation)');
        end

        function value = get.ExpectationCt(obj)
            value = obj.ExpectationC';
        end

        function value = get.ExpectationCtC(obj)
            if obj.cols
                value = zeros(obj.Size, obj.Size);
                for i = 1:obj.Size
                    for j = 1:obj.Size
                        % [NOTE] E[Wt * W] is a matrix where each entry is e.g. E[w1^T *
                        % w2] (w1 and w2 are columns of the matrix W and there
                        % are independent random variable). Given the
                        % independence this is E[w1^T]E[w2], but when i == j
                        % then the independence doesn't hold and in that case
                        % we should use ExpectationXtX from the
                        % GaussianDistribution class.
                        if i ~= j
                            value(i, j) = obj.distributions(i).ExpectationXt * ...
                            obj.distributions(j).Expectation;
                        else 
                            value(i, j) = obj.distributions(i).ExpectationXtX;
                        end
                    end
                end
            else
            dim = obj.distributions(1).dim; % They are all of the same dimension
            value = zeros(dim, dim);
                for i = 1:obj.Size
                    value = value + obj.distributions(i).ExpectationXXt;
                end
            end
        end

        function value = get.H(obj)
            value = cell(1, obj.Size);
            for i = 1:obj.Size
                value{i} = obj.distributions(i).H;
            end
        end

        function value = get.HC(obj)
            value = 0;
            for i = 1:obj.Size
                value = value + obj.distributions(i).H;
            end
        end
    end
end