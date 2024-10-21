classdef Datasets
    methods (Static)
        % X is [N x D]
        function [X, D] = generateBPCA(N, D, stdDevs)
            if nargin == 0
                N = 300;
                D = 10;
                % Standard deviations along orthogonal directions
                % stdDevs = [5, 4, 3, 2, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
                stdDevs = [5, 5, 5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
            elseif nargin ~= 3
                error(['##### ERROR IN THE CLASS ' mfilename('class') ': Pass either none or all three arguments.']);
            end

            if length(stdDevs) ~= D
                error(['##### ERROR IN THE CLASS ' mfilename('class') ': stdDevs array must have D elements.']);
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
    
        

        % Standard scaler
        function X_scaled = standardScaler(X)
            % X: The input matrix of size [N x D]
            mean_X = mean(X);

            std_X = std(X);
            
            % Standardize the data
            X_scaled = (X - mean_X) ./ std_X;
        end
    

        % Returns N of one digit and N of other with labels
        % One digit is in the top rows the other is in the bottom rows
        function [X, y] = getMNISTData(folderName, binaryLabels, N)
            labelsFileName = Utility.ternary(binaryLabels, '/binaryLabels.tsv', '/continuousLabels.tsv');

            % featuresInCols = true for both X1 and X2
            X1 = readmatrix(['datasets/', folderName, '/pixels.tsv'], 'FileType', 'text'); % [N x D1];
            X2 = readmatrix(['datasets/', folderName, labelsFileName], 'FileType', 'text'); % [N x D2], D2 = 1;
            
            if size(X1, 1) ~= size(X2, 1)
                error(['##### ERROR IN THE CLASS ' mfilename('class') ': X1 and X2 don''t have the same number of observations']);
            end

            % Combine top and bottom rows to form the test set
            X = [X1(1:N, :); X1(end-N+1:end, :)];
            y = [X2(1:N, :); X2(end-N+1:end, :)];
        end


        % This function is used to generate a randomized train/test split
        % for the experiments involving MNIST datasets given in 'mnist18'
        % and 'mnist38'
        function [X1_train, X2_train, X1_test, X2_test] = trainTestSplitMNIST(folderName, binaryLabels, testPerc)
            % Optional parameter: testPerc
            if nargin < 3
                testPerc = 0.3;
            end

            labelsFileName = Utility.ternary(binaryLabels, '/binaryLabels.tsv', '/continuousLabels.tsv');

            % featuresInCols = true for both X1 and X2
            X1 = readmatrix(['datasets/', folderName, '/pixels.tsv'], 'FileType', 'text'); % [N x D1];
            X2 = readmatrix(['datasets/', folderName, labelsFileName], 'FileType', 'text'); % [N x D2], D2 = 1;
            
            if size(X1, 1) ~= size(X2, 1)
                error(['##### ERROR IN THE CLASS ' mfilename('class') ': X1 and X2 don''t have the same number of observations']);
            end
            
            testPerc = testPerc / 2; % The dataset is organize so one digit is on top and other in bottom rows
            N = size(X1, 1);

            % Calculate the number of rows for the test set
            testRowsCnt = round(testPerc * N);
            
            % Combine top and bottom rows to form the test set
            X1_test = [X1(1:testRowsCnt, :); X1(end-testRowsCnt+1:end, :)];
            
            % Get the remaining rows for the train set
            X1_train = X1(testRowsCnt+1:end-testRowsCnt, :);
            
            % Combine top and bottom rows to form the test set
            X2_test = [X2(1:testRowsCnt, :); X2(end-testRowsCnt+1:end, :)];
            
            % Get the remaining rows for the train set
            X2_train = X2(testRowsCnt+1:end-testRowsCnt, :);
            
            shuffle_indices = randperm(size(X1_train, 1));
            
            % Shuffle rows of both matrices using the same permutation
            X1_train = X1_train(shuffle_indices, :);
            X2_train = X2_train(shuffle_indices, :);
        end



        function data = generateGFA_2G()
            % Spherical noise - check the definition for tau!
            %
            %
            N_tr = 400; 
            N_te = 100; % For predicting views
            N = N_tr + N_te;
        
            M = 2;
            D = [50, 30]; % Dimensions of the views
        
            % Latent factors
            K = 4;
            Z = zeros(K, N);
        
            % Generate latent factors - row by row (in total K rows)
            n = 1:N;
            Z(1, :) = sin((n)/(N/20));
            Z(2, :) = cos((n)/(N/20));
            Z(3, :) = 2 * ((n)/N-0.5); 
            Z(4, :) = normrnd(0, 1, [N, 1]);
        
            % Noise precisions (tau) for each view; Spherical noise;
            tau = {5; 10};
            
            % Alpha (alpha) for each view K values
            alpha = zeros(K, M);
            alpha(:, 1) = [1e6, 1, 1e6, 1];
            alpha(:, 2) = [1, 1, 1, 1e6];
            
            % W matrix for each view
            W = cell(M, 1);
        
            % Train and test observarions: for all views M
            X_tr = cell(M, 1);
            X_te = cell(M, 1);
            
            for m = 1:M
                Dm = D(m);
                W{m} = zeros(Dm, K);
                for k = 1:K
                    % Generate w_k (kth column of W) from p(W | alpha)
                    alpha_m_k = alpha(k, m); % alpha (precision) for the view 'm' and column 'k'
        
                    % normrnd(MU,SIGMA) returns an array of random numbers chosen from a
                    % normal distribution with mean MU and standard deviation SIGMA.
                    W{m}(:, k) = normrnd(0, 1/sqrt(alpha_m_k), [Dm, 1]);
                end
        
                % Generate X for view 'm'
                X = W{m} * Z + normrnd(0, 1/sqrt(tau{m}), [Dm, N]);
        
                % Get training and test data
                X_tr{m} = X(:, 1:N_tr);
                X_te{m} = X(:, N_tr+1:end);
            end
            
            % Latent variables for training the model    
            Z = Z(:, 1:N_tr);
            
            % Store data and model parameters            
            data.X_tr = X_tr;
            data.X_te = X_te;
            data.W = W;
            data.Z = Z;
            data.tau = tau;
            data.alpha = alpha;
            data.K = K;
        end

    

        function data = generateGFA_3G()
            % Spherical noise - check the definition for tau!
            %
            %
            N_tr = 400; 
            N_te = 100; % For predicting views
            N = N_tr + N_te;
        
            M = 3;
            D = [50, 30, 20]; % Dimensions of the views
        
            % Latent factors
            K = 4;
            Z = zeros(K, N);
        
            % Generate latent factors - row by row (in total K rows)
            n = 1:N;
            Z(1, :) = sin((n)/(N/20));
            Z(2, :) = cos((n)/(N/20));
            Z(3, :) = 2 * ((n)/N-0.5); 
            Z(4, :) = normrnd(0, 1, [N, 1]);
        
            % Noise precisions (tau) for each view; Spherical noise;
            tau = {5; 10; 8};
            
            % Alpha (alpha) for each view K values
            alpha = zeros(K, M);
            alpha(:, 1) = [1, 1, 1e6, 1];
            alpha(:, 2) = [1, 1, 1, 1e6];
            alpha(:, 3) = [1, 1e6, 1, 1e6];
            
            % W matrix for each view
            W = cell(M, 1);
        
            % Train and test observarions: for all views M
            X_tr = cell(M, 1);
            X_te = cell(M, 1);
            
            for m = 1:M
                Dm = D(m);
                W{m} = zeros(Dm, K);
                for k = 1:K
                    % Generate w_k (kth column of W) from p(W | alpha)
                    alpha_m_k = alpha(k, m); % alpha (precision) for the view 'm' and column 'k'
        
                    % normrnd(MU,SIGMA) returns an array of random numbers chosen from a
                    % normal distribution with mean MU and standard deviation SIGMA.
                    W{m}(:, k) = normrnd(0, 1/sqrt(alpha_m_k), [Dm, 1]);
                end
        
                % Generate X for view 'm'
                X = W{m} * Z + normrnd(0, 1/sqrt(tau{m}), [Dm, N]);
        
                % Get training and test data
                X_tr{m} = X(:, 1:N_tr);
                X_te{m} = X(:, N_tr+1:end);
            end
            
            % Latent variables for training the model    
            Z = Z(:, 1:N_tr);
            
            % Store data and model parameters            
            data.X_tr = X_tr;
            data.X_te = X_te;
            data.W = W;
            data.Z = Z;
            data.tau = tau;
            data.alpha = alpha;
            data.K = K;
        end
        
        % Used for repeated train-test splits
        function [X_tr, y_tr, X_te, y_te] = trainTestSplit(X, y, verbose, testPerc)
            % % Set the random seed for reproducibility
            % rng(1);
            if nargin < 4
                testPerc = 0.2;
                if nargin < 3
                    verbose = false;
                end
            end
            cv = cvpartition(y, 'HoldOut', testPerc);

            % Train and test split
            X_tr = X(training(cv), :);
            y_tr = y(training(cv), :);
            X_te = X(test(cv), :);
            y_te = y(test(cv), :);

            if verbose
                disp('Training class distribution:');
                tabulate(y_tr)
                disp('Test class distribution:');
                tabulate(y_te)
            end
        end
    end
end
