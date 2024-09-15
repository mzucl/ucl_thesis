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





        %% Code transformations for the binary extension
        % ln q(Z) - second line
        function testIdentityBinary1_1(testCase)
            D = 10;
            a = Utility.generateRandomIntMatrix(D, 1);
            h = Utility.generateRandomIntMatrix(D, 1);

            sum = 0;
            for d = 1:D
                sum = sum + h(d) * (a(d))^2;
            end

            altRes = a' * diag(h) * a;

            testCase.verifyEqual(sum, altRes);
        end

        % For vectorization in q(Z) for Bohning bound
        function testIdentityBinary1_2(testCase)
            D = 10;
            K = 8;
            W = Utility.generateRandomIntMatrix(D, K);
            D = 1/4 * eye(D);

            testCase.verifyEqual(W' * D * W, 1/4 * (W' * W));
        end

        % For vectorization in q(mu) for Bohning bound
        function testIdentityBinary1_2_1(testCase)
            D = 10;
            N = 100;
            H = 1/4 * ones(D, N);
            
            res1 = diag(H * ones(N, 1));
            res2 = N/4 * eye(D);

            testCase.verifyEqual(res1, res2);
            testCase.verifyEqual(Utility.matrixInverse(res1), 4/N * eye(D));
        end

        % For vectorization in q(mu) for Bohning bound
        function testIdentityBinary1_2_2(testCase)
            D = 10;
            N = 100;
            H = 1/4 * ones(D, N);
            WZ = Utility.generateRandomIntMatrix(D, N);
            
            res1 = H .* WZ;
            res2 = 1/4 * WZ;

            testCase.verifyEqual(res1, res2);
        end

        % For vectorization in q(Z) for Bohning bound
        function testIdentityBinary1_3(testCase)
            D = 10;
            N = 20;
            T = Utility.generateRandomIntMatrix(D, N);
            H = 1/4 * ones(D, N);
            mu = Utility.generateRandomIntMatrix(D, 1);

            testCase.verifyEqual(T - H .* mu, T - 1/4 * mu);
        end

        % E[Y^TDY] -> there is a section in the thesis with these kind of
        % expectations;
        % Part 1
        function testIdentityBinary2_1(testCase)
            D = 10;
            K = 8;
            MU = Utility.generateRandomIntMatrix(K, D);
            diagMatrix = diag(Utility.generateRandomIntMatrix(D, 1));
            Sigma = Utility.generateRandomIntMatrix(K, K);

            sum = zeros(K, K);
            for d = 1:D
                sum = sum + diagMatrix(d, d) * (MU(:, d) * MU(:, d)' + Sigma);
            end

            altRes =  MU * diagMatrix * MU' + trace(diagMatrix) * Sigma;

            testCase.verifyEqual(sum, altRes);
        end

        % Part 2
        function testIdentityBinary2_2(testCase)
            K = 8;
            mu_1 = Utility.generateRandomIntMatrix(K, 1);
            Sigma_1 = Utility.generateRandomIntMatrix(K, K);
            D = diag(Utility.generateRandomIntMatrix(K, 1));

            res1 = trace((mu_1 * mu_1' + Sigma_1) * D);
            res2 = mu_1' * D * mu_1 + trace(Sigma_1 * D);

            testCase.verifyEqual(res1, res2);
        end

        % Part 3
        function testIdentityBinary2_3(testCase)
            K = 8;
            N = 100;
            MU = Utility.generateRandomIntMatrix(K, N);
            Sigma = Utility.generateRandomIntMatrix(K, K);
            D = diag(Utility.generateRandomIntMatrix(K, 1));

            % Test 2: D is diagonal
            res1 = zeros(N, N);
            for i = 1:N
                for j = 1:N
                    if i == j
                        res1(i, j) = MU(:, i)' * D * MU(:, j) + trace(Sigma * D);
                    else
                        res1(i, j) = MU(:, i)' * D * MU(:, j);
                    end
                end
            end

            res2 = MU' * D * MU + trace(Sigma * D) * eye(N);

            testCase.verifyEqual(res1, res2);
        end




        % 'Monster' expectation in q(W) for the binary extension - part 1
        % Tr(Y' * (H .* Y)) = sum(sum(...))
        function testIdentityBinary3_1(testCase)
            D = 10;
            N = 100;
            Y = Utility.generateRandomIntMatrix(D, N); % Y = WZ
            H = Utility.generateRandomIntMatrix(D, N);

            res1 = trace(Y' * (H .* Y));
            res2 = 0;

            for n = 1:N
                for d = 1:D
                    res2 = res2 + Y(d, n)^2 * H(d, n);
                end
            end
            testCase.verifyEqual(res1, res2);
        end

        % 'Monster' expectation in q(W) for the binary extension - part 1
        % Tr(Y' * (H .* Y)) = sum(sum(...))
        function testIdentityBinary3_1_1(testCase)
            D = 10;
            N = 100;
            H = Utility.generateRandomIntMatrix(D, N); % A = WZ
            Tr = Utility.generateRandomIntMatrix(D, N); % Matrix of traces

            res1 = sum(sum(H .* Tr));

            res2 = 0;
            for n = 1:N
                for d = 1:D
                    res2 = res2 + H(d, n) * Tr(d, n);
                end
            end

            testCase.verifyEqual(res1, res2);
        end

        % 'Monster' expectation in q(W) for the binary extension - part 2
        % Computing Tr(mu_zn * mu_zn^T * Sigma_wd) as a DxN matrix, where
        % each element is the value of the trace for specific d and n
        function testIdentityBinary3_2(testCase)
            D = 20;
            N = 100;
            K = 12;
            H = Utility.generateRandomIntMatrix(D, N);
            MU_Z = Utility.generateRandomIntMatrix(K, N); % mu_zn are the cols

            % randi(10, [3, 3, 5]) -> generate a 3x3x5 array of random integers between 1 and 10
            Sigma_W = randi(10, [K, K, D]);

            % Non-vectorized approach
            res1 = zeros(D, N);
            for n = 1:N
                for d = 1:D
                    res1(d, n) = trace(MU_Z(:, n) * MU_Z(:, n)' * Sigma_W(:, :, d));
                end
            end
            
            % Vectorized approach
            res2 = Utility.flatten3DTo2D(Sigma_W)' * Utility.flatten3DTo2D(Utility.outerProduct3D(MU_Z));
            
            testCase.verifyEqual(res1, res2);
        end

        % 'Monster' expectation in q(W) for the binary extension - part 3
        % Computing Tr(Sigma_zn * mu_wd * mu_wd^T) as a DxN matrix, where
        % each element is the value of the trace for specific d and n
        function testIdentityBinary3_3(testCase)
            D = 20;
            N = 100;
            K = 12;
            H = Utility.generateRandomIntMatrix(D, N);
            MU_Z = Utility.generateRandomIntMatrix(K, N); % mu_zn are the cols

            % randi(10, [3, 3, 5]) -> generate a 3x3x5 array of random integers between 1 and 10
            Sigma_W = randi(10, [K, K, D]);

            % Non-vectorized approach
            res1 = zeros(D, N);
            for n = 1:N
                for d = 1:D
                    res1(d, n) = trace(MU_Z(:, n) * MU_Z(:, n)' * Sigma_W(:, :, d));
                end
            end
            
            % Vectorized approach
            res2 = Utility.flatten3DTo2D(Sigma_W)' * Utility.flatten3DTo2D(Utility.outerProduct3D(MU_Z));
            
            testCase.verifyEqual(res1, res2);
        end
    end
end
