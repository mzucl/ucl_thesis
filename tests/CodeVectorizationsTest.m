classdef CodeVectorizationsTest < matlab.unittest.TestCase
    methods (Static)
        function [X, W, T] = generateViewMatrices(D, N, K)
            X = Utility.generateRandomIntMatrix(D, N);
            W = Utility.generateRandomIntMatrix(D, K);
            T = diag(diag(Utility.generateRandomIntMatrix(D, D)));
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





        %% Code vectorization GFA model


        % Vectorization for qZUpdate: mu
        function testIdentityGFAVectorization1(testCase)
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

        % Vectorization for qZUpdate: sigma
        function testIdentityGFAVectorization2(testCase)
            % Setup
            Ds = [5, 10, 12]; % View dimensions
            K = 80;
            N = 10000;
            M = length(Ds);

            [X1, W1, T1] = CodeVectorizationsTest.generateViewMatrices(Ds(1), N, K);
            [X2, W2, T2] = CodeVectorizationsTest.generateViewMatrices(Ds(2), N, K);
            [X3, W3, T3] = CodeVectorizationsTest.generateViewMatrices(Ds(3), N, K);

            W = blkdiag(W1, W2, W3, W1, W2, W3, W1, W2, W3, W1, W2, W3, W1, W2, W3, W1, W2, W3);
            T = blkdiag(T1, T2, T3, T1, T2, T3, T1, T2, T3, T1, T2, T3, T1, T2, T3, T1, T2, T3);

            % Vectorized code
            tic;
            temp = W' * T * W;
            % sigmaZ_Vect = Utility.matrixInverse(eye(K) + temp(1:K, 1:K) + temp(K + 1:2 * K, K + 1:  2 * K) + ...
            %     temp(2 * K + 1:end, 2 * K + 1:end));

            elapsedTime = toc;
            fprintf('Elapsed time : %.4f seconds\n', elapsedTime);

            % Non-vectorized code
            tic;
            sigmaZ = Utility.matrixInverse(eye(K) + W1' * T1 * W1 + W2' * T2 * W2 + W3' * T3 * W3 + ...
                W1' * T1 * W1 + W2' * T2 * W2 + W3' * T3 * W3 + W1' * T1 * W1 + W2' * T2 * W2 + W3' * T3 * W3 + ...
                W1' * T1 * W1 + W2' * T2 * W2 + W3' * T3 * W3 + W1' * T1 * W1 + W2' * T2 * W2 + W3' * T3 * W3 + ...
                W1' * T1 * W1 + W2' * T2 * W2 + W3' * T3 * W3);
            elapsedTime = toc;
            fprintf('Elapsed time : %.4f seconds\n', elapsedTime);
            
            % testCase.verifyEqual(sigmaZ, sigmaZ_Vect);
        end
    end
end