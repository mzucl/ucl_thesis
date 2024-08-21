classdef Datasets
    methods (Static)
        % X is [N x D]
        function [X, D] = generateBPCA(N, D, stdDevs)
            if nargin == 0
                N = 100;
                D = 10;
                % Standard deviations along orthogonal directions
                % stdDevs = [5, 4, 3, 2, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
                stdDevs = [5, 5, 5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
            elseif nargin ~= 3
                error(['##### ERROR IN THE CLASS Datasets' ': Pass either none or all three arguments.']);
            end

            if length(stdDevs) ~= D
                error(['##### ERROR IN THE CLASS Utility' ': stdDevs array must have D elements.']);
            end
            
            cov = diag(stdDevs.^2);
            
            % Generate random orthogonal matrix for rotation
            [U, ~] = qr(randn(D));
            
            % The true covariance matrix after rotation
            TrueCovMatrix = U * cov * U';
            
            % Mean vector (centered at origin)
            mu = zeros(1, D);
            
            % Generate dataset from multivariate Gaussian distribution
            X = mvnrnd(mu, TrueCovMatrix, N);
        end
    
        

        %% Standard scaler
        function X_scaled = standardScaler(X)
            % X: The input matrix of size [N x D]
            mean_X = mean(X);

            std_X = std(X);
            
            % Standardize the data
            X_scaled = (X - mean_X) ./ std_X;
        end
    end
end
