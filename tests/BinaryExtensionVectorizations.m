classdef BinaryExtensionVectorizations < matlab.unittest.TestCase    
    methods (Test)
        %% Block 4.1
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
        function testIdentityBinary2_1(testCase)
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
        function testIdentityBinary2_2(testCase)
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

        function testIdentityBinary2_3(testCase)
            D = 20;
            N = 100;
            H = 1/4 * ones(D, N);

            testCase.verifyEqual(diag(H * ones(N, 1)), N / 4 * eye(D, D));
        end



        % Transformation for qWUpdate: term with zn^Tzn
        function testIdentityBinary3_1(testCase)
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

        % Transformation for qWUpdate: term with zn^T
        function testIdentityBinary3_2(testCase)
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

        % Vectorization in Block 4.6
        function testIdentityBinary3_3(testCase)
            D = 20;
            N = 100;
            K = 8;

            sigma = Utility.generateRandomIntMatrix(K, K);
            Z = Utility.generateRandomIntMatrix(K, N);
            R = Utility.generateRandomIntMatrix(D, N);

            % Vectorized approach
            res = sigma * Z * R';

            for d = 1:D
                resSingle = sigma * Z * R(d, :)';
                testCase.verifyEqual(resSingle, res(:, d));
            end
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
