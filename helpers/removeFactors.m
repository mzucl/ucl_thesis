function [Z_new, W_new, alpha_new] = removeFactors(Z, W, alpha, threshold)
    if nargin < 4
        threshold = Constants.LATENT_FACTORS_THRESHOLD;
    end
    % Calculate the average of the square of elements for each row of Z
    avgSquare = mean(Z.^2, 2);

    % Find indices of rows where the average of squares is below a threshold 
    removeIdx = find(avgSquare < threshold);

    % Remove those rows from Z, corresponding columns from W, and elements from alpha
    Z_new = Z;
    W_new = W;
    alpha_new = alpha;

    Z_new(removeIdx, :) = [];
    W_new(:, removeIdx) = [];
    alpha_new(removeIdx) = [];
end
