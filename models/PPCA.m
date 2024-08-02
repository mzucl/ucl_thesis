function W = PPCA(X, k)
    [~, d] = size(X);

    if nargin < 2
        k = d - 1; % Maximum value for k
    end

    % Center the data
    X_centered = X - mean(X, 1);
    
    covarianceMatrix = cov(X_centered);
    
    % Eigenvalue decomposition
    [eigVectors, eigValues] = eig(covarianceMatrix);
    
    eigValues = diag(eigValues);
    [sortedEigValues, idx] = sort(eigValues, 'descend');
    
    % ML solution for sigma^2 (Bishop: 12.46)
    sigmaSq = 1./(d - k) * sum(sortedEigValues(k + 1:end));

    % ML solution for W (Bishop: 12.45); R = I;
    W = eigVectors(:, idx(1:k)) * ( diag(sortedEigValues(1:k)) - ...
        sigmaSq * eye(k) )^(1./2);
end