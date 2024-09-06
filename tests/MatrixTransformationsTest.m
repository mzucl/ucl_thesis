classdef MatrixTransformationsTest < matlab.unittest.TestCase
    methods (Test)
        %% Code transformations for the BPCA model
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





        %% Code transformations for the GFA model 
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
    
    



        %% Code transformations for the new GFA model (spherical noise)
        % ln q(W) terms when we expand the quadratic form (term 3 is the
        % same as in the diagonal noise GFA model)
        function testIdentityGFANew1_1(testCase)
            % Setup
            D = 20;
            K = 10;
            N = 50;

            X = Utility.generateRandomIntMatrix(D, N);
            W = Utility.generateRandomIntMatrix(D, K);
            Z = Utility.generateRandomIntMatrix(K, N);
            mu = Utility.generateRandomIntMatrix(D, 1);

            sum1 = 0;
            for n = 1:N
                sum1 = sum1 + (X(:, n) - mu)' * W * Z(:, n);
            end

            sum2 = 0;
            for d = 1:D
                sum2 = sum2 + W(d, :) * Z * (X(d, :)' - mu(d) * ones(N, 1));
            end

            testCase.verifyEqual(sum1, sum2);
        end

        function testIdentityGFANew1_2(testCase)
            % Setup
            D = 20;
            K = 10;
            N = 50;

            W = Utility.generateRandomIntMatrix(D, K);
            Z = Utility.generateRandomIntMatrix(K, N);

            sum1 = 0;
            for n = 1:N
                sum1 = sum1 + Z(:, n)' * (W' * W) * Z(:, n);
            end

            sum2 = 0;
            for d = 1:D
                sum2 = sum2 + W(d, :) * (Z * Z') * W(d, :)';
            end

            testCase.verifyEqual(sum1, sum2);
        end
    
        % ln q(tau) first term when we expand the quadratic form
        function testIdentityGFANew2_1(testCase)
            % Setup
            N = 50;
            D = 20;

            X = Utility.generateRandomIntMatrix(D, N);
            mu = Utility.generateRandomIntMatrix(D, 1);

            sum1 = 0;
            for n = 1:N
                sum1 = sum1 + (X(:, n) - mu)' * (X(:, n) - mu);
            end

            vect = trace(X' * X) - 2 * mu' * X * ones(N, 1) + N * (mu' * mu);

            testCase.verifyEqual(sum1, vect);
        end
    
        % ln q(tau) the whole expression
        function testIdentityGFANew2_2(testCase)
            % Setup
            K = 10;
            N = 50;
            D = 20;

            X = Utility.generateRandomIntMatrix(D, N);
            W = Utility.generateRandomIntMatrix(D, K);
            Z = Utility.generateRandomIntMatrix(K, N);
            mu = Utility.generateRandomIntMatrix(D, 1);

            sum1 = 0;
            for n = 1:N
                sum1 = sum1 + (X(:, n) - mu - W * Z(:, n))' * (X(:, n) - mu - W * Z(:, n));
            end

            vect = trace(X' * X) - 2 * mu' * X * ones(N, 1) + N * (mu' * mu) ...
                - 2 * trace(W * Z * (X - mu)') + trace( ...
                (W' * W) * (Z * Z'));

            testCase.verifyEqual(sum1, vect);
        end

        % ln q(mu) the whole expression
        function testIdentityGFANew3(testCase)
            % Setup
            K = 10;
            N = 50;
            D = 20;

            X = Utility.generateRandomIntMatrix(D, N);
            W = Utility.generateRandomIntMatrix(D, K);
            Z = Utility.generateRandomIntMatrix(K, N);
            mu = Utility.generateRandomIntMatrix(D, 1);
            beta = randi(15);
            tau = randi(100);

            sum = 0;
            for n = 1:N
                sum = sum - 2 * mu' * X(:, n) + 2 * mu' * W * Z(:, n) + mu' * mu;
            end

            startEq = -tau/2 * sum - beta/2 * (mu' * mu);

            endEq = -1/2 * mu' * (beta + N * tau) * eye(D) * mu + ...
                mu' * (tau * X * ones(N,1) - tau * W * Z * ones(N, 1));

            testCase.verifyEqual(startEq, endEq);
        end

        % q(tau) update
        %% trace(W * Z * (X - mu)') == sum(sum((W * Z) .* (X - mu)))
        %% X * ones(N, 1) = sum(X, 2)
        function testIdentityGFANew4(testCase)
            % Setup
            K = 10;
            N = 50;
            D = 20;

            X = Utility.generateRandomIntMatrix(D, N);
            W = Utility.generateRandomIntMatrix(D, K);
            Z = Utility.generateRandomIntMatrix(K, N);
            mu = Utility.generateRandomIntMatrix(D, 1);

            tr1 = trace(W * Z * (X - mu)');
            
            tr2 = sum(sum((W * Z) .* (X - mu)));

            testCase.verifyEqual(tr1, tr2);
            testCase.verifyEqual(X * ones(N, 1), sum(X, 2));
        end
    end
end
