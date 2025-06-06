function [W, sigmaSq] = PPCA(X, K)
    %PPCA - Computes the Maximum Likelihood solution of Probabilistic PCA.
    %
    % Description:
    %   This function computes the loading matrix and isotropic noise variance 
    %   for a Probabilistic PCA (PPCA) model using the closed-form solution 
    %   described in Bishop's "Pattern Recognition and Machine Learning" 
    %   (Equations 12.45 and 12.46). The input data is first centered, and then 
    %   eigendecomposition is performed on the sample covariance matrix.
    %
    % Input:
    %   X  - A [D x N] data matrix where each column represents an observation.
    %   K  - (Optional) The number of principal components to retain. If not 
    %        specified, the default is D - 1.
    %
    % Output:
    %   W        - A [D x K] matrix of principal component directions (loadings).
    %   sigmaSq  - A scalar representing the estimated isotropic noise variance.
    CustomError.validateNumberOfParameters(nargin, 1, 2);

    [D, ~] = size(X);
    if nargin < 2
        K = D - 1; % Maximum value for K
    end

    % Center the data
    X_centered = X - mean(X, 2);
    
    % [NOTE] For matrices, where each row is an observation, and each column a variable, 
    % cov(X) is the covariance matrix.
    covarianceMatrix = cov(X_centered');
    
    % Eigenvalue decomposition
    [eigVectors, eigValues] = eig(covarianceMatrix);
    
    eigValues = diag(eigValues);
    [sortedEigValues, idx] = sort(eigValues, 'descend');
    
    % ML solution for sigma^2 (Bishop: 12.46)
    sigmaSq = 1./(D - K) * sum(sortedEigValues(K + 1:end));

    % ML solution for W (Bishop: 12.45); R = I;
    W = eigVectors(:, idx(1:K)) * ( diag(sortedEigValues(1:K)) - ...
        sigmaSq * eye(K) )^(1./2);
end