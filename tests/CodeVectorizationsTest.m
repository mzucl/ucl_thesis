classdef CodeVectorizationsTest < matlab.unittest.TestCase
    methods (Static)
        function [X, W, T] = generateViewMatrices(D, N, K)
            X = Utility.generateRandomIntMatrix(D, N);
            W = Utility.generateRandomIntMatrix(D, K);
            T = diag(diag(Utility.generateRandomIntMatrix(D, D)));
        end

        % Generate view for the new GFA model: mu included, noise
        % spherical;
        function [X, W, mu, tau] = generateViewMatricesGFANew(D, N, K)
            X = Utility.generateRandomIntMatrix(D, N);
            W = Utility.generateRandomIntMatrix(D, K);
            mu = Utility.generateRandomIntMatrix(D, 1);
            tau = randi(15); % Drawn from Unif[1:15]
        end
    end




    
    methods (Test)
        %% Code vectorization BPCA model
        % Vectorization for qZUpdate
        function testIdentityBPCAVectorization1(testCase)
            D = 20;
            K = 10;
            N = 50;
            % <tau> * \Sigma_z * <W>^T
            X = Utility.generateRandomIntMatrix(D, N);
            A = Utility.generateRandomIntMatrix(K, D);
            mu = Utility.generateRandomIntMatrix(D, 1);

            % Vectorized code
            MU = A * (X - mu);

            % Non-vectorized code
            for n = 1:N
                mu_n = A * (X(:, n) - mu);
                testCase.verifyEqual(mu_n, MU(:, n));
            end
        end


        % Vectorization for qMuUpdate
        function testIdentityBPCAVectorization2(testCase)
            % Setup
            D = 20;
            K = 10;
            N = 50;
            X = Utility.generateRandomIntMatrix(D, N);
            W = Utility.generateRandomIntMatrix(D, K);
            Z = Utility.generateRandomIntMatrix(K, N);

            sum = 0;

            % Non-vectorized code
            for n = 1:N
                sum = sum + (X(:, n) - W * Z(:, n));
            end

            res = (X - W * Z) * ones(N, 1);

            testCase.verifyEqual(sum, res);
        end


        % Vectorization for qWUpdate
        function testIdentityBPCAVectorization3(testCase)
            % Setup
            D = 20;
            K = 10;
            N = 50;
            X = Utility.generateRandomIntMatrix(D, N);
            Z = Utility.generateRandomIntMatrix(K, N);
            mu = Utility.generateRandomIntMatrix(D, 1);

            % Vectorized code
            MU = Z * (X' - mu');

            for d = 1:size(MU, 2) % D
                sum = 0;
                % Non-vectorized code
                for n = 1:N
                    sum = sum + Z(:, n) * (X(d, n) - mu(d));
                end
    
                testCase.verifyEqual(sum, MU(:, d));
            end
        end


        % Vectorization for qTauUpdate - first step
        function testIdentityBPCAVectorization4(testCase)
            % Setup
            D = 20;
            K = 10;
            N = 50;
            X = Utility.generateRandomIntMatrix(D, N);
            Z = Utility.generateRandomIntMatrix(K, N);
            mu = Utility.generateRandomIntMatrix(D, 1);
            W = Utility.generateRandomIntMatrix(D, K);

            sum = 0;

            for n = 1:N
                sum = sum + X(:, n)' * X(:, n) + mu' * mu + trace(W' * W * Z(:, n) * Z(:, n)') + ...
                    2 * mu' * W * Z(:, n) - 2 * X(:, n)' * W * Z(:, n) - 2 * X(:, n)' * mu;
            end

            sum = 1/2 * sum;

            res = 1/2 * trace(X' * X) + N/2 * (mu' * mu) + 1/2 * trace((W' * W) * (Z * Z')) + ...
                mu' * W * Z * ones(N, 1) - trace(W * Z * X') - mu' * X * ones(N, 1);

            testCase.verifyEqual(sum, res);
        end

        
        % Vectorization for qTauUpdate - second step
        function testIdentityBPCAVectorization5(testCase)
            % Setup
            D = 20;
            K = 10;
            N = 50;
            X = Utility.generateRandomIntMatrix(D, N);
            Z = Utility.generateRandomIntMatrix(K, N);
            mu = Utility.generateRandomIntMatrix(D, 1);
            W = Utility.generateRandomIntMatrix(D, K);

            % Vectorized, but not optimized
            res1 = 1/2 * trace(X' * X) + N/2 * dot(mu, mu) + 1/2 * trace((W' * W) * (Z * Z')) + ...
                mu' * W * Z * ones(N, 1) - trace(W * Z * X') - mu' * X * ones(N, 1);
            
            WtW = (W' * W)'; % This can be optimized, but the purpose here is to test identity
            ZZt = Z * Z';
            WZ = (W * Z);

            res2 = 1/2 * dot(X(:), X(:)) + N/2 * dot(mu, mu) + 1/2 * dot(WtW(:), ZZt(:)) + ...
                - dot(X(:), WZ(:)) +  mu' * (W * Z - X) * ones(N, 1);

            testCase.verifyEqual(res1, res2);
        end





        %% Code vectorization for the GFA model
        % Vectorization for qZUpdate: mu
        function testIdentityGFAVectorization1(testCase)
            % Setup
            Ds = [5, 10, 12]; % View dimensions
            K = 10;
            N = 100;

            [X1, W1, T1] = CodeVectorizationsTest.generateViewMatrices(Ds(1), N, K);
            [X2, W2, T2] = CodeVectorizationsTest.generateViewMatrices(Ds(2), N, K);
            [X3, W3, T3] = CodeVectorizationsTest.generateViewMatrices(Ds(3), N, K);

            sigmaZ = Utility.generateRandomIntMatrix(K, K);
            
            % Vectorized code
            MU = sigmaZ * (W1' * T1 * X1 + W2' * T2 * X2 + W3' * T3 * X3);

            % Non-vectorized code
            for n = 1:N
                mu = W1' * T1 * X1(:, n) + W2' * T2 * X2(:, n) + W3' * T3 * X3(:, n);
                mu = sigmaZ * mu;

                % Compare
                testCase.verifyEqual(MU(:, n), mu);
            end
        end


        % Vectorization for qZUpdate: mu - full vectorization over the
        % views as well!
        function testIdentityGFAVectorization1_2(testCase)
            % Setup
            Ds = [5, 10, 12]; % View dimensions
            K = 10;
            N = 100;
            M = length(Ds);

            [X1, W1, T1] = CodeVectorizationsTest.generateViewMatrices(Ds(1), N, K);
            [X2, W2, T2] = CodeVectorizationsTest.generateViewMatrices(Ds(2), N, K);
            [X3, W3, T3] = CodeVectorizationsTest.generateViewMatrices(Ds(3), N, K);

            X = [X1; X2; X3];
            W = blkdiag(W1, W2, W3);
            T = blkdiag(T1, T2, T3);

            sigmaZ = Utility.generateRandomIntMatrix(K, K);
            
            % Vectorized code
            MU = W' * T * X; % blkdiag(W1' * T1, W2' * T2, W3' * T3)
            MU = reshape(MU, K, M, []);

            MU = sigmaZ * squeeze(sum(MU, 2));

            % Non-vectorized code
            for n = 1:N
                mu = W1' * T1 * X1(:, n) + W2' * T2 * X2(:, n) + W3' * T3 * X3(:, n);
                mu = sigmaZ * mu;

                % Compare
                testCase.verifyEqual(MU(:, n), mu);
            end
        end


        % Vectorization for qWUpdate: mu
        function testIdentityGFAVectorization2_1(testCase)
            % Setup
            D = 5; % View dimensions
            K = 10;
            N = 100;

            [X, ~, T] = CodeVectorizationsTest.generateViewMatrices(D, N, K);
            Z = Utility.generateRandomIntMatrix(K, N);

            sigma = repmat(Utility.generateRandomSPDMatrix(K), 1, 1, D);
            
            % Vectorized code
            V = reshape(Z * X' * T, K, 1, D); % Columns of the matrix will be in the third dimension

            MU = squeeze(pagemtimes(sigma, V));

            % Non-vectorized code
            for d = 1:D
                mu = sigma(:, :, d) * T(d, d) * Z * X(d, :)'; 

                % Compare
                diffSqNorm = norm(MU(:, d) - mu);
                testCase.verifyTrue(diffSqNorm < 1e-9);
            end
        end


        % Vectorization for qWUpdate: sigma
        % Step 1: Vectorize the equation with the inverse
        % Step 2: Inverse each of the matrices inside the multidimensional
        % array
        function testIdentityGFAVectorization2_2(testCase)
            % Setup
            D = 5; % View dimensions
            K = 3;
            N = 10;

            [~, ~, T] = CodeVectorizationsTest.generateViewMatrices(D, N, K);
            Z = Utility.generateRandomIntMatrix(K, N);
            alpha = diag(Utility.generateRandomIntMatrix(K, K));

            T_3D = reshape(diag(T), 1, 1, []);

            ZZt_3D = repmat(Z * Z', 1, 1, D);
            
            alpha_3D = repmat(alpha, 1, 1, D);

            SIGMA = pagemtimes(ZZt_3D, T_3D) + alpha_3D;

            % Matrix inverse
            SIGMA_inv = arrayfun(@(i) Utility.matrixInverse(SIGMA(:,:,i)), ...
                1:size(SIGMA, 3), 'UniformOutput', false);

            % Convert cell array to multidimensional matrix
            SIGMA_inv = cat(3, SIGMA_inv{:});

            % Non-vectorized code
            for d = 1:D
                sigma = T(d, d) * (Z * Z') + alpha;

                % Compare without inverse
                diffSqNorm = norm(SIGMA(:, :, d) - sigma);
                testCase.verifyTrue(diffSqNorm < 1e-9);

                % Compare inverse
                diffSqNorm = norm(SIGMA_inv(:, :, d) - Utility.matrixInverse(sigma));
                testCase.verifyTrue(diffSqNorm < 1e-9);
            end
        end
    

        % Vectorization for binary - Jaakkola bound (SSHIBA approach)
        function testIdentityGFAVectorization3(testCase)
            % Setup
            D = 10; % View dimensions
            K = 10;
            N = 100;

            [~, W, T] = CodeVectorizationsTest.generateViewMatrices(D, N, K);
            Z = Utility.generateRandomIntMatrix(K, N);
            Sigma = Utility.generateRandomSPDMatrix(D);

            T_obs = Utility.generateRandomBinaryMatrix(D, N);

            % Vectorized code
            MU = Sigma * (T * W * Z + T_obs - 1/2 * ones(D, N));
            
            % Non-vectorized code
            for n = 1:N
                mu = Sigma * (T * W * Z(:, n) + T_obs(:, n) - 1/2 * ones(D, 1));
                
                % Compare
                diffSqNorm = norm(MU(:, n) - mu);
                testCase.verifyTrue(diffSqNorm < 1e-9);
            end
        end





        %% Code vectorization for the new GFA model
        % Vectorization for qZUpdate: mu
        function testIdentityGFANewVectorization1(testCase)
            % Setup
            Ds = [5, 10, 12]; % View dimensions
            K = 10;
            N = 100;

            
            [X1, W1, mu1, tau1] = CodeVectorizationsTest.generateViewMatricesGFANew(Ds(1), N, K);
            [X2, W2, mu2, tau2] = CodeVectorizationsTest.generateViewMatricesGFANew(Ds(2), N, K);
            [X3, W3, mu3, tau3] = CodeVectorizationsTest.generateViewMatricesGFANew(Ds(3), N, K);

            sigmaZ = Utility.generateRandomIntMatrix(K, K);
            
            % Vectorized code
            MU = sigmaZ * (tau1 * W1' * (X1 - mu1) + tau2 * W2' * (X2 - mu2) + ...
                tau3 * W3' * (X3 - mu3));

            % Non-vectorized code
            for n = 1:N
                mu = tau1 * W1' * (X1(:, n) - mu1) + tau2 * W2' * (X2(:, n) - mu2) + ...
                    tau3 * W3' * (X3(:, n) - mu3);
                mu = sigmaZ * mu;

                % Compare
                testCase.verifyEqual(MU(:, n), mu);
            end
        end
    
        % Vectorization for qWUpdate: mu
        function testIdentityGFANewVectorization2(testCase)
            % Setup
            D = 5; % View dimension
            K = 10;
            N = 100;

            [X, ~, mu, tau] = CodeVectorizationsTest.generateViewMatricesGFANew(D, N, K);
            Z = Utility.generateRandomIntMatrix(K, N);
            
            sigmaW = Utility.generateRandomIntMatrix(K, K);

            % Vectorized code
            MU = sigmaW * tau * Z * (X - mu)';

            % Non-vectorized code
            for d = 1:D
                mu_Wd = sigmaW * tau * Z * (X(d, :)' - mu(d));

                % Compare
                testCase.verifyEqual(MU(:, d), mu_Wd);
            end
        end
    
    

    
    
        %% Code vectorization for the binary extension
        % Vectorization for qZUpdate: mu (part 2)
        function testIdentityBinary1(testCase)
            D = 20;
            N = 100;
            mu = Utility.generateRandomIntMatrix(D, 1);
            H = Utility.generateRandomIntMatrix(D, N);

            % Vectorized approach
            res = H .* mu;

            for n = 1:N
                resSingle = diag(H(:, n)) * mu;
                testCase.verifyEqual(resSingle, res(:, n));
            end
        end

        % Vectorization for qMuUpdate: term between mu^T and mu
        function testIdentityBinary2(testCase)
            D = 20;
            N = 100;
            H = Utility.generateRandomIntMatrix(D, N);

            % Vectorized approach
            res = diag(H * ones(N, 1));

            sum = zeros(D, D);
            for n = 1:N
                sum = sum + diag(H(:, n)); 
            end
            testCase.verifyEqual(res, sum);
        end

        % Vectorization for qMuUpdate: term with mu^T
        function testIdentityBinary3(testCase)
            D = 20;
            N = 100;
            K = 8;
            H = Utility.generateRandomIntMatrix(D, N);
            W = Utility.generateRandomIntMatrix(D, K);
            Z = Utility.generateRandomIntMatrix(K, N);

            % Vectorized approach
            res = H .* (W * Z) * ones(N, 1);

            sum = zeros(D, 1);

            for n = 1:N
                sum = sum + diag(H(:, n)) * W * Z(:, n);
            end
            testCase.verifyEqual(sum, res);
        end

        % Vectorization for qWUpdate: term with zn^Tzn
        function testIdentityBinary4_1(testCase)
            D = 20;
            N = 100;
            K = 8;

            H = Utility.generateRandomIntMatrix(D, N);
            W = Utility.generateRandomIntMatrix(D, K);
            Z = Utility.generateRandomIntMatrix(K, N);

            % Sum over n
            sum1 = 0;
            for n = 1:N
                sum1 = sum1 + Z(:, n)' * W' * diag(H(:, n)) * W * Z(:, n);
            end
            
            % Sum over d
            sum2 = 0;
            for d = 1:D
                sum2 = sum2 + W(d, :) * Z * diag(H(d, :)) * Z' * W(d, :)';
            end

            testCase.verifyEqual(sum1, sum2);
        end

        % Vectorization for qWUpdate: term with zn^T
        function testIdentityBinary4_2(testCase)
            D = 20;
            N = 100;
            K = 8;

            W = Utility.generateRandomIntMatrix(D, K);
            Z = Utility.generateRandomIntMatrix(K, N);
            R = Utility.generateRandomIntMatrix(D, N);

            % Sum over n
            sum1 = 0;
            for n = 1:N
                sum1 = sum1 + Z(:, n)' * W' * R(:, n);
            end

            % Sum over d
            sum2 = 0;
            for d = 1:D
                sum2 = sum2 + W(d, :) * Z * R(d, :)';
            end

            testCase.verifyEqual(sum1, sum2);
        end

        % Vectorization for the const term in ELBO
        function testIdentityBinary5_1(testCase)
            D = 20;
            N = 50;
            Xi = Utility.generateRandomIntMatrix(D, N);
            C = Utility.generateRandomIntMatrix(D, N);
            G = Utility.generateRandomIntMatrix(D, N);
            H = Utility.generateRandomIntMatrix(D, N);

            sum1 = 0;
            for n = 1:N
                for d = 1:D
                    sum1 = sum1 - C(d, n) + Xi(d, n) * G(d, n) - 1/2 * Xi(d, n)^2 * H(d, n);
                end
            end

            res = sum(sum(-C + Xi .* G - 1/2 * Xi.^2 .* H));

            testCase.verifyEqual(sum1, res);
        end

        % Vectorization for the phi_1 term in ELBO
        function testIdentityBinary5_2(testCase)
            D = 20;
            N = 100;
            K = 8;

            W = Utility.generateRandomIntMatrix(D, K);
            Z = Utility.generateRandomIntMatrix(K, N);
            T = Utility.generateRandomIntMatrix(D, N);
            X = Utility.generateRandomIntMatrix(D, N);
            mu = Utility.generateRandomIntMatrix(D, 1);

            sum = 0;
            for n = 1:N
                sum = sum + (Z(:, n)' * W'  + mu') * (X(:, n) + T(:, n));
            end

            % Vectorized form
            res = trace((Z' * W' + mu') * (X + T));

            testCase.verifyEqual(sum, res);
        end

        % Vectorization for the phi_2 term in ELBO: term independent of
        % zn^T
        function testIdentityBinary5_3_1(testCase)
            D = 20;
            N = 100;

            H = Utility.generateRandomIntMatrix(D, N);
            mu = Utility.generateRandomIntMatrix(D, 1);

            sum_ = 0;
            for n = 1:N
                sum_ = sum_ + mu' * diag(H(:, n)) * mu;
            end

            % Vectorized form
            res = mu' * diag(H * ones(N, 1)) * mu;

            testCase.verifyEqual(sum_, res);
        end

        % Vectorization for the phi_2 term in ELBO: linear term
        function testIdentityBinary5_3_2(testCase)
            D = 20;
            N = 100;
            K = 8;

            W = Utility.generateRandomIntMatrix(D, K);
            Z = Utility.generateRandomIntMatrix(K, N);
            H = Utility.generateRandomIntMatrix(D, N);
            mu = Utility.generateRandomIntMatrix(D, 1);

            sum_ = 0;
            for n = 1:N
                sum_ = sum_ + Z(:, n)' * W' * diag(H(:, n)) * mu;
            end

            % Vectorized form
            res = ones(1, N) * ((Z' * W') .* H') * mu;

            testCase.verifyEqual(sum_, res);
        end

        % Vectorization for the phi_2 term in ELBO: square term
        function testIdentityBinary5_3_3(testCase)
            D = 20;
            N = 100;
            K = 8;

            W = Utility.generateRandomIntMatrix(D, K);
            Z = Utility.generateRandomIntMatrix(K, N);
            H = Utility.generateRandomIntMatrix(D, N);
            mu = Utility.generateRandomIntMatrix(D, 1);

            sum_ = 0;
            for n = 1:N
                sum_ = sum_ + Z(:, n)' * W' * diag(H(:, n)) * W * Z(:, n);
            end

            % Vectorized form
            res = trace((W * Z)' * (H .* (W * Z)));

            testCase.verifyEqual(sum_, res);
        end


        % Vectorization for q(W) update, vectorization for mu of W under
        % Bohning bound
        function testIdentityBinary6(testCase)
            D = 20;
            N = 100;
            K = 8;

            Sigma = Utility.generateRandomIntMatrix(K, K);
            R = Utility.generateRandomIntMatrix(D, N);
            Z = Utility.generateRandomIntMatrix(K, N);

            % Vectorized approach
            MU = Sigma * Z * R';

            for d = 1:D
                mu_d = Sigma * Z * R(d, :)';
                testCase.verifyEqual(mu_d, MU(:, d));
            end
        end
    end
end
