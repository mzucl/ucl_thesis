classdef GaussianDistributionContainer < handle
    properties
        distributions       % A list of distributions inside the container
        cols                %   - 'true' when distributions describe columns of a matrix
                            %   - 'false' when they describe rows. This is important
                            % when the dependent properties are calculated
    end
    
    properties (Dependent)
        Size                % Number of distributions in the container
        Expectation         % Under the assumption of independence of each component; 
                            % This is a cell array where each entry is an
                            % expectation of one component;
        PriorPrecision      % Similar to above, each entry is precision of the prior of one component
        ExpectationXXt      % Similar to above, each entry is E[XXt]
        ExpectationXtX      % Similar to above, each entry is E[XtX], which is equivalent to E[|X|^2]
        ExpectationC        % 'C' stands for container, and this is expectation of the whole container
        ExpectationCt 
        ExpectationCtC      % This corresponds to e.g. E[W^TW]
        ExpectationCCt      % Sum of E[XXt] of each distribution
        H                   % Similar to above, each entry is entropy of a single component
        HC                  % Entropy of the collection

        ExpectationLnP      % Cell array where each entry is ExpectationLnP

        ExpectationLnPC     % The sum of all entries in ExpectationLnP
    end

    methods(Access = private)
        % Helper method that throws an error if index is not valid
        function validateIndex(obj, idx)
            if idx < 1 || idx > obj.Size 
                error(['##### ERROR IN THE CLASS ' class(obj) ': Index out of range.']); 
            end
        end
    end

    %% Options for the constructor GaussianDistributionContainer
    % 3 PARAMETERS
    % (numDistributions, prior, cols)   -> Creates a container with
    % the specified format ('cols') and specified number of components
    % ('numDistributions') where each component is a GaussianDistribution
    % instance defined by the prior (SINGLE prior used for all components)

    methods
        %% Constructors
        function obj = GaussianDistributionContainer(numDistributions, prior, cols)
            % [NOTE] For now we have only three parameters constructor 
            if nargin < 3
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
            else
                validParams = isnumeric(numDistributions) && isscalar(prior) && ...
                    Utility.isNaNOrInstanceOf(prior, 'GaussianDistribution') && ...
                    ~Utility.isNaN(prior) && islogical(cols);

                if ~validParams
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Invalid parameters.']);
                end

                obj.cols = cols;
                obj.distributions = repmat(GaussianDistribution(), numDistributions, 1); % Preallocate
                for i = 1:numDistributions
                    obj.distributions(i) = GaussianDistribution(prior);
                end
            end
        end
       

        
        %% Update and retrieve methods
        function dist = getDistribution(obj, idx)
            % Returns the distribution at index 'idx'
            obj.validateIndex(idx);

            dist = obj.distributions(idx);
        end

        % [NOTE] This is used in the scenario when e.g. W is stored in the rows
        % format, but we need expectation of the squared norm of a column
        function res = getExpectationOfColumnsNormSq(obj)
            if obj.cols == true
                numOfCols = obj.Size;
                res = zeros(1, numOfCols);
                for idx = 1:numOfCols
                    res(idx) = obj.ExpectationXtX{idx};
                end
            else
                numOfCols = obj.distributions(1).dim;
                res = zeros(1, numOfCols);
                for k = 1:numOfCols
                    for d = 1:obj.Size
                        res(k) = res(k) + obj.distributions(d).mu(k) ^ 2 + ...
                            obj.distributions(d).cov(k, k);
                    end
                end
            end
        end

        function obj = updateDistribution(obj, idx, dist)
            obj.validateIndex(idx);

            % Updates the distribution at index 'idx'
            obj.distributions(idx) = dist.copy();
        end

        function obj = updateDistributionParams(obj, idx, mu, cov)
            if nargin < 4
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
            end
            obj.validateIndex(idx);

            obj.distributions(idx).updateParameters(mu, cov);
        end

        function obj = updateDistributionMu(obj, idx, mu)
            if nargin < 3
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
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
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
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
        % use columns for mean and expectations, 'cols' parameter is
        % not related to that (to the single distr. inside)
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

        % [!!!] Dependent on 'cols'
        function value = get.ExpectationC(obj)
            value = Utility.ternary(obj.cols, cell2mat(obj.Expectation), cell2mat(obj.Expectation)');
        end

        function value = get.ExpectationCt(obj)
            value = obj.ExpectationC';
        end

        % TODO (high): Refactor this and the one below this - the code is
        % the same!
        function value = get.ExpectationCCt(obj)
            if obj.cols
                dim = obj.distributions(1).dim; % They are all of the same dimension
                value = zeros(dim, dim);
                for i = 1:obj.Size
                    value = value + obj.distributions(i).ExpectationXXt;
                end
            else
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
            end
        end

        % [!!!] Dependent on 'cols'
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

        function value = get.PriorPrecision(obj)
            value = cell(1, obj.Size);
            for i = 1:obj.Size
                value{i} = obj.distributions(i).PriorPrecision;
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

        function value = get.ExpectationLnP(obj)
            value = cell(1, obj.Size);
            for i = 1:obj.Size
                value{i} = obj.distributions(i).ExpectationLnP;
            end
        end

        function value = get.ExpectationLnPC(obj)
            value = 0;
            for i = 1:obj.Size
                value = value + obj.distributions(i).ExpectationLnP;
            end
        end
    end
end