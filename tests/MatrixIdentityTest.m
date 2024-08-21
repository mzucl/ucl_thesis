classdef MatrixIdentityTest < matlab.unittest.TestCase
    methods (Test)
        function testIdentity1(testCase)
            %% sum(an^T * an) = Tr(A^T * A)

            % Test matrices of different shapes and sizes
            mVals = [5, 5, 2, 10];
            nVals = [5, 3, 7, 10];
            for i = 1:length(mVals)
                A = Utility.generateRandomIntMatrix(mVals(i), nVals(i));

                sum = 0;
                for col = 1:nVals(i)
                    sum = sum + A(:, col)' * A(:, col);
                end
                testCase.verifyEqual(sum, trace(A' * A));
            end 
        end



        function testIdentity2(testCase)
            %% sum(an^T * v) = vT * A * 1 (vector of 1s)

            % Test matrices of different shapes and sizes
            mVals = [2, 5, 5, 2, 10];
            nVals = [2, 5, 3, 7, 10];
            for i = 1:length(mVals)
                A = Utility.generateRandomIntMatrix(mVals(i), nVals(i));

                v = Utility.generateRandomIntMatrix(mVals(i), 1); % Mx1
                o = ones(nVals(i), 1); % Nx1 vector of '1'

                sum = 0;
                for col = 1:nVals(i)
                    sum = sum + A(:, col)' * v;
                end

                testCase.verifyEqual(sum, v' * A * o);
            end 
        end



        function testIdentity3(testCase)
            %% sum(an^T * W * bn) = Tr(WBA^T)

            % Test matrices of different shapes and sizes
            mVals = [2, 5, 5, 2, 10];
            nVals = [2, 5, 3, 7, 10];
            kVals = [2, 4, 4, 9, 12];

            for i = 1:length(mVals)
                A = Utility.generateRandomIntMatrix(mVals(i), nVals(i));
                B = Utility.generateRandomIntMatrix(kVals(i), nVals(i));
                W = Utility.generateRandomIntMatrix(mVals(i), kVals(i));

                sum = 0;
                for col = 1:nVals(i)
                    sum = sum + A(:, col)' * W * B(:, col);
                end

                testCase.verifyEqual(sum, trace(W * B * A'));
            end 
        end



        function testIdentity4(testCase)
            %% sum(an^T * W * v) = 1^T * A^T * W * v
            % This one is a special case of one of the above

            % Test matrices of different shapes and sizes
            mVals = [2, 5, 5, 2, 10];
            nVals = [2, 5, 3, 7, 10];
            kVals = [2, 4, 4, 9, 12];

            for i = 1:length(mVals)
                A = Utility.generateRandomIntMatrix(mVals(i), nVals(i));
                W = Utility.generateRandomIntMatrix(mVals(i), kVals(i));
                v = Utility.generateRandomIntMatrix(kVals(i), 1);

                o = ones(nVals(i), 1); % Nx1 vector of '1'

                sum = 0;
                for col = 1:nVals(i)
                    sum = sum + A(:, col)' * W * v;
                end

                testCase.verifyEqual(sum, o' * A' * W * v);
            end 
        end



        function testIdentity5(testCase)
            %% sum(an * an^T) = AA^T
            % Test matrices of different shapes and sizes
            mVals = [5, 5, 2, 10];
            nVals = [5, 3, 7, 10];
            for i = 1:length(mVals)
                A = Utility.generateRandomIntMatrix(mVals(i), nVals(i));

                sum = 0;
                for col = 1:nVals(i)
                    sum = sum + A(:, col) * A(:, col)';
                end
                testCase.verifyEqual(sum, A*A');
            end 
        end

        function testIdentity5_1(testCase)
            %% trace(A^T * A) = sum of squared elements of A
            mVals = [5, 5, 2, 10];
            nVals = [5, 3, 7, 10];
            for i = 1:length(mVals)
                A = Utility.generateRandomIntMatrix(mVals(i), nVals(i));

                testCase.verifyEqual(trace(A' * A), dot(A(:), A(:)));
            end 
        end


        
        
        %% ln q(W) terms when we expand the quadratic form for GFA model
        function testIdentity6(testCase)
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

        function testIdentity7(testCase)
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

        function testIdentity8(testCase)
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



        %% Expectation in q(alpha) for GFA model
        function testIdentity9(testCase)
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

        function testIdentity10(testCase)
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
            temp = diag(X * X' - 2 * W * Z * X' + W * Z * Z' * W');
            for d = 1:D
                sum2 = sum2 + T(d, d) * temp(d);
            end

            testCase.verifyEqual(sum1, sum2);
        end


        %% Sum in E[ln p(X | Z, W, mu, tau)]
        function testIdentity11(testCase)
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
                + 2 * mu' * W * Z * ones(N, 1) + trace(W' * W * Z * Z') + N * mu' * mu;

            testCase.verifyEqual(sum, expr);
        end

        

        %% Code vectorization
        % Vectorization for qZUpdate in BPCA model
        function testIdentity12(testCase)
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

        % Vectorization for qMuUpdate in BPCA model
        function testIdentity13(testCase)
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

        % Vectorization for qWUpdate in BPCA model
        function testIdentity14(testCase)
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

        % Vectorization for qTauUpdate in BPCA model - first step
        function testIdentity15(testCase)
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
    end
end
