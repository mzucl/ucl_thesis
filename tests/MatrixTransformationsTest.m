classdef MatrixTransformationsTest < matlab.unittest.TestCase
    methods (Test)
        %% Code transformations for BPCA model
        % Sum in E[ln p(X | Z, W, mu, tau)]
        function testIdentityBPCA1(testCase)
            % Setup
            D = 20;
            K = 10;
            N = 50;
            X = Utility.generateRandomIntMatrix(D, N);
            W = Utility.generateRandomIntMatrix(D, K);
            Z = Utility.generateRandomIntMatrix(K, N);
            mu = Utility.generateRandomIntMatrix(D, 1);

            sum = 0;
            for n = 1:N
                el = X(:, n) - W * Z(:, n) - mu;
                sum = sum + el' * el;
            end

            expr = trace(X' * X) - 2 * trace(W * Z * X') - 2 * mu' * X * ones(N, 1) ...
                + 2 * mu' * W * Z * ones(N, 1) + trace(W' * W * (Z * Z')) + N * (mu' * mu);

            testCase.verifyEqual(sum, expr);
        end





        %% Code transformations for GFA model 
        % ln q(W) terms when we expand the quadratic form
        function testIdentityGFA1_1(testCase)
            % Setup
            D = 20;
            K = 10;
            N = 50;

            X = Utility.generateRandomIntMatrix(D, N);
            W = Utility.generateRandomIntMatrix(D, K);
            T = diag(diag(Utility.generateRandomIntMatrix(D, D)));
            Z = Utility.generateRandomIntMatrix(K, N);

            sum1 = 0;
            for n = 1:N
                sum1 = sum1 + X(:, n)' * T * W * Z(:, n);
            end

            sum2 = 0;
            for d = 1:D
                sum2 = sum2 + W(d, :) * T(d, d) * Z * X(d, :)';
            end

            testCase.verifyEqual(sum1, sum2);
        end

        function testIdentityGFA1_2(testCase)
            % Setup
            D = 20;
            K = 10;
            N = 50;

            X = Utility.generateRandomIntMatrix(D, N);
            W = Utility.generateRandomIntMatrix(D, K);
            T = diag(diag(Utility.generateRandomIntMatrix(D, D)));
            Z = Utility.generateRandomIntMatrix(K, N);

            sum1 = 0;
            for n = 1:N
                sum1 = sum1 + Z(:, n)' * W' * T * W * Z(:, n);
            end

            sum2 = 0;
            for d = 1:D
                sum2 = sum2 + W(d, :) * T(d, d) * (Z * Z') * W(d, :)';
            end

            testCase.verifyEqual(sum1, sum2);
        end

        function testIdentityGFA1_3(testCase)
            % Setup
            D = 20;
            K = 10;
            N = 50;
            W = Utility.generateRandomIntMatrix(D, K);
            alpha = Utility.generateRandomIntMatrix(K, 1);

            sum1 = 0;
            for k = 1:K
                sum1 = sum1 + alpha(k) * W(:, k)' * W(:, k);
            end

            sum2 = 0;
            for d = 1:D
                sum2 = sum2 + W(d, :) * diag(alpha) * W(d, :)';
            end

            testCase.verifyEqual(sum1, sum2);
        end


        % Expectation in q(tau) - part 1
        function testIdentityGFA2_1(testCase)
            % Setup
            D = 20;
            K = 10;
            N = 50;
            X = Utility.generateRandomIntMatrix(D, N);
            W = Utility.generateRandomIntMatrix(D, K);
            Z = Utility.generateRandomIntMatrix(K, N);

            sum = 0;
            for n = 1:N
                sum = sum + X(:, n) * X(:, n)' - 2 * W * Z(:, n) * X(:, n)' + ...
                W * Z(:, n) * Z(:, n)' * W';
            end

            expr = X * X' - 2 * W * Z * X' + W * (Z * Z') * W';

            testCase.verifyEqual(sum, expr);
        end


        % Expectation in q(tau) - part 2
        function testIdentityGFA2_2(testCase)
            % Setup
            D = 20;
            K = 10;
            N = 50;
            X = Utility.generateRandomIntMatrix(D, N);
            W = Utility.generateRandomIntMatrix(D, K);
            T = diag(diag(Utility.generateRandomIntMatrix(D, D)));
            Z = Utility.generateRandomIntMatrix(K, N);

            sum1 = 0;
            for n = 1:N
                sum1 = sum1 + (X(:, n) - W * Z(:, n))' * T * (X(:, n) - W * Z(:, n));
            end

            sum2 = 0;
            temp = diag(X * X' - 2 * W * Z * X' + W * (Z * Z') * W');
            for d = 1:D
                sum2 = sum2 + T(d, d) * temp(d);
            end

            testCase.verifyEqual(sum1, sum2);
        end     
    end
end
