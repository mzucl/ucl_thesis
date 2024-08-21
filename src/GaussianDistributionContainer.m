classdef GaussianDistributionContainer < handle
    properties
        ds                  % A list of distributions inside the container
        cols                %   - 'true' when distributions describe columns of a matrix
                            %   - 'false' when they describe rows. This is important
                            % when the dependent properties are calculated
    end
    
    properties (Dependent)
        Size                % Number of distributions in the container
        E                   % Under the assumption of independence of each component; 
                            % This is a cell array where each entry is an
                            % expectation of one component;
        PPrec               % Similar to above, each entry is precision of the prior of one component
        E_XXt               % Similar to above, each entry is E[XXt]
        E_XtX               % Similar to above, each entry is E[XtX], which is equivalent to E[|X|^2]
        EC                  % 'C' stands for container, and this is expectation of the whole container
        E_Ct 
        E_CtC               % This corresponds to e.g. E[W^TW]
        E_CCt               % Sum of E[XXt] of each distribution
        H                   % Similar to above, each entry is entropy of a single component
        HC                  % Entropy of the collection
        E_LnP               % Cell array where each entry is E_LnP
        E_LnPC              % The sum of all entries in E_LnP
        Tr_CtC              % Tr(X^TX), but computed without computing all elements of the matrix
                            % TODO(high): only computed for the cols format
                            % because it is needed for Z^TZ
    end

    properties(Access = private)
        expCInit
    end

    methods(Access = private)
        % Helper method that throws an error if index is not valid
        function validateIndex(obj, idx)
            if idx < 1 || idx > obj.Size 
                error(['##### ERROR IN THE CLASS ' class(obj) ': Index out of range.']); 
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

    %% Options for the constructor GaussianDistributionContainer
    % 3 PARAMETERS
    % (numDistributions, prior, cols)   -> Creates a container with
    % the specified format ('cols') and specified number of components
    % ('numDistributions') where each component is a GaussianDistribution
    % instance defined by the 'priors' (either SINGLE prior used for all components or 
    % each component has its OWN prior, in which case the length(priors) must be equal to
    % the 'numDistributions').
    methods
        %% Constructors
        function obj = GaussianDistributionContainer(numDistributions, priors, cols)
            % [NOTE] 'priors' can't be NaN because it is used to infer the
            % information about the distributions!
            if nargin < 3
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
            else
                validParams = isnumeric(numDistributions) && islogical(cols) && ...
                    Utility.areAllInstancesOf(priors, 'GaussianDistribution') && ...
                    (isscalar(priors) || length(priors) == numDistributions);

                if ~validParams
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Invalid parameters.']);
                end

                obj.cols = cols;
                obj.ds = repmat(GaussianDistribution(), numDistributions, 1); % Preallocate
                for i = 1:numDistributions
                    obj.ds(i) = Utility.ternaryOpt(isscalar(priors), ...
                            @() GaussianDistribution(priors), @() GaussianDistribution(priors(i)));
                end
            end

            % Set initial expectation to the real expectation
            obj.setExpCInit(obj.EC);
        end
       

        
        %% Update and retrieve methods
        function dist = getDistribution(obj, idx)
            % Returns the distribution at index 'idx'
            obj.validateIndex(idx);

            dist = obj.ds(idx);
        end

        % [NOTE] This is used in the scenario when e.g. W is stored in the rows
        % format, but we need expectation of the squared norm of a column
        function res = getExpectationOfColumnsNormSq(obj)
            if obj.cols == true
                numOfCols = obj.Size;
                res = zeros(numOfCols, 1);
                for idx = 1:numOfCols
                    res(idx) = obj.E_XtX{idx};
                end
            else
                numOfCols = obj.ds(1).dim;
                res = zeros(numOfCols, 1);
                for k = 1:numOfCols
                    for d = 1:obj.Size
                        res(k) = res(k) + obj.ds(d).mu(k) ^ 2 + ...
                            obj.ds(d).cov(k, k);
                    end
                end
            end
        end

        function obj = updateDistribution(obj, idx, dist)
            obj.validateIndex(idx);

            % Updates the distribution at index 'idx'
            obj.ds(idx) = dist.copy();
        end

        function obj = updateDistributionParams(obj, idx, mu, cov)
            if nargin < 4
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
            end
            obj.validateIndex(idx);

            obj.ds(idx).updateParameters(mu, cov);
        end

        function obj = updateDistributionMu(obj, idx, mu)
            if nargin < 3
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
            end
            obj.validateIndex(idx);

            obj.ds(idx).updateMu(mu); % Will raise an exception if we try to change dimension
        end

        function obj = updateDistributionCovariance(obj, idx, cov)
            if nargin < 3
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
            end
            obj.validateIndex(idx);

            obj.ds(idx).updateCovariance(cov); % Will raise an exception if we try to change dimension
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
                obj.ds(i).updateCovariance(cov);
            end
        end

        function obj = rotateAllDistributionsMu(obj, R)
            if nargin < 2
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
            end
            for i = 1:obj.Size
                obj.ds(i).rotateMu(R);
            end
        end

        function obj = rotateAllDistributionsCovariance(obj, R, RtCovR)
            % Update covariance of all distributions
            % if(RtCovR = true) -> R' * cov * R
            % else              -> R * cov * R'
            if nargin < 3
                error(['##### ERROR IN THE CLASS ' class(obj) ': Too few arguments passed.']);
            end
            
            for i=1:obj.Size
                obj.ds(i).rotateCovariance(R, RtCovR);
            end
        end

        function removeDimensions(obj, indices)
            if nargin < 2 || isempty(indices)
                return; % No change
            end

            for i = 1:obj.Size
                obj.ds(i).removeDimensions(indices);
            end
        end


        %% Setters
        function obj = setExpCInit(obj, value)
            validValue = obj.cols == true && isequal(size(value), [obj.ds(1).dim, obj.Size]) || ...
                obj.cols == false && isequal(size(value), [obj.Size, obj.ds(1).dim]);

            if ~validValue
                error(['##### ERROR IN THE CLASS ' class(obj) ': Number of elements in the expectation must be equal to the number of ' ...
                    'components. ']);
            end
            obj.expCInit = value;
        end



        %% Getters
        % 'expInit' is a private property -> needs a getter for the access
        function value = getExpCInit(obj)
            value = obj.expCInit;
        end

        function value = get.Size(obj)
            value = length(obj.ds);
        end

        % [NOTE] The value for 'cols' is irrelevant for the properties that
        % are of type cell, because they are used with indexing to access
        % components values. And for each component (multivariate Gauss) we
        % use columns for mean and expectations, 'cols' parameter is
        % not related to that (to the single distr. inside)
        function value = get.E(obj)
            value = cell(1, obj.Size);
            for i = 1:obj.Size
                value{i} = obj.ds(i).E;
            end
        end

        function value = get.E_XXt(obj)
            value = cell(1, obj.Size);
            for i = 1:obj.Size
                value{i} = obj.ds(i).E_XXt;
            end
        end
        
        function value = get.E_XtX(obj)
            value = cell(1, obj.Size);
            for i = 1:obj.Size
                value{i} = obj.ds(i).E_XtX;
            end
        end

        % [!!!] Dependent on 'cols'
        function value = get.EC(obj)
            value = Utility.ternary(obj.cols, cell2mat(obj.E), cell2mat(obj.E)');
        end

        function value = get.E_Ct(obj)
            value = obj.EC';
        end

        function value = get.E_CCt(obj)
            value = obj.getContainerExpectations(false); % CtC
        end

        function value = get.E_CtC(obj)
            value = obj.getContainerExpectations(true); % CtC
        end

        % TODO (high): This is implemented just for the cols format as it
        % was needed for the Tr(<Z^TZ>). Implement this properly and add
        % tests; The purpose of it is to optimize the calculation, where we
        % don't calculate all elements of the Z^TZ product, just diagonal
        function value = get.Tr_CtC(obj)
            value = 0;
            for i = 1:obj.Size
                value = value + obj.ds(i).E_XtX;
            end
        end

        function value = get.PPrec(obj)
            value = cell(1, obj.Size);
            for i = 1:obj.Size
                value{i} = obj.ds(i).PPrec;
            end
        end

        function value = get.H(obj)
            value = cell(1, obj.Size);
            for i = 1:obj.Size
                value{i} = obj.ds(i).H;
            end
        end

        function value = get.HC(obj)
            value = 0;
            for i = 1:obj.Size
                value = value + obj.ds(i).H;
            end
        end

        function value = get.E_LnP(obj)
            value = cell(1, obj.Size);
            for i = 1:obj.Size
                value{i} = obj.ds(i).E_LnP;
            end
        end

        function value = get.E_LnPC(obj)
            value = 0;
            for i = 1:obj.Size
                value = value + obj.ds(i).E_LnP;
            end
        end
    end
end