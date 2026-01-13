classdef MatrixIdentityTest < matlab.unittest.TestCase
    methods (Test)
        %% [A.27] sum(an^T * an) = Tr(A^T * A)
        function testIdentity1(testCase)
            % Test matrices of different shapes and sizes
            mVals = [5, 5, 2, 10];
            nVals = [5, 3, 7, 10];
            for i = 1:length(mVals)
                A = RandomMatrices.intMatrix(mVals(i), nVals(i));

                sum = 0;
                for col = 1:nVals(i)
                    sum = sum + A(:, col)' * A(:, col);
                end
                testCase.verifyEqual(sum, trace(A' * A));
            end 
        end


        %% [A.29] sum(an^T * v) = vT * A * 1 (vector of 1s)
        function testIdentity2(testCase)
            % Test matrices of different shapes and sizes
            mVals = [2, 5, 5, 2, 10];
            nVals = [2, 5, 3, 7, 10];
            for i = 1:length(mVals)
                A = RandomMatrices.intMatrix(mVals(i), nVals(i));

                v = RandomMatrices.intMatrix(mVals(i), 1); % Mx1
                o = ones(nVals(i), 1); % Nx1 vector of '1'

                sum = 0;
                for col = 1:nVals(i)
                    sum = sum + A(:, col)' * v;
                end

                testCase.verifyEqual(sum, v' * A * o);
            end 
        end


        %% [A.30] sum(an^T * W * bn) = Tr(WBA^T)
        function testIdentity3(testCase)
            % Test matrices of different shapes and sizes
            mVals = [2, 5, 5, 2, 10];
            nVals = [2, 5, 3, 7, 10];
            kVals = [2, 4, 4, 9, 12];

            for i = 1:length(mVals)
                A = RandomMatrices.intMatrix(mVals(i), nVals(i));
                B = RandomMatrices.intMatrix(kVals(i), nVals(i));
                W = RandomMatrices.intMatrix(mVals(i), kVals(i));

                sum = 0;
                for col = 1:nVals(i)
                    sum = sum + A(:, col)' * W * B(:, col);
                end

                testCase.verifyEqual(sum, trace(W * B * A'));
            end 
        end


        %% [A.29 where v = W * v] sum(an^T * W * v) = 1^T * A^T * W * v
        % This one is a special case of one of the above
        function testIdentity4(testCase)
            % Test matrices of different shapes and sizes
            mVals = [2, 5, 5, 2, 10];
            nVals = [2, 5, 3, 7, 10];
            kVals = [2, 4, 4, 9, 12];

            for i = 1:length(mVals)
                A = RandomMatrices.intMatrix(mVals(i), nVals(i));
                W = RandomMatrices.intMatrix(mVals(i), kVals(i));
                v = RandomMatrices.intMatrix(kVals(i), 1);

                o = ones(nVals(i), 1); % Nx1 vector of '1'

                sum = 0;
                for col = 1:nVals(i)
                    sum = sum + A(:, col)' * W * v;
                end

                testCase.verifyEqual(sum, o' * A' * W * v);
            end 
        end


        %% [A.28 where B = A] sum(an * an^T) = AA^T
        function testIdentity5(testCase)
            % Test matrices of different shapes and sizes
            mVals = [5, 5, 2, 10];
            nVals = [5, 3, 7, 10];
            for i = 1:length(mVals)
                A = RandomMatrices.intMatrix(mVals(i), nVals(i));

                sum = 0;
                for col = 1:nVals(i)
                    sum = sum + A(:, col) * A(:, col)';
                end
                testCase.verifyEqual(sum, A*A');
            end 
        end

       
        %% [A.31] trace(A^T * A) = sum of squared elements of A
        function testIdentity6(testCase)
            mVals = [5, 5, 2, 10];
            nVals = [5, 3, 7, 10];
            for i = 1:length(mVals)
                A = RandomMatrices.intMatrix(mVals(i), nVals(i));

                testCase.verifyEqual(trace(A' * A), dot(A(:), A(:)));
            end 
        end
        

        %% [A.32] trace(A * B) = dot(A'(:), B(:))
        function testIdentity7(testCase)
            mVals = [5, 5, 2, 10];
            nVals = [5, 3, 7, 10];
            for i = 1:length(mVals)
                A = RandomMatrices.intMatrix(mVals(i), nVals(i));
                B = RandomMatrices.intMatrix(nVals(i), mVals(i));

                A_tr = A';
                testCase.verifyEqual(trace(A * B), dot(A_tr(:), B(:)));
            end 
        end


        %% [A.33] trace(A * B * C) = dot(AB'(:), C(:))
        % this is the same as trace(CAB) and trace(BCA), so we can choose 
        % which two matrices to multiply
        function testIdentity8(testCase)
            mVals = [5, 5, 2, 10];
            nVals = [5, 3, 7, 10];
            kVals = [5, 5, 12, 15];
            for i = 1:length(mVals)
                A = RandomMatrices.intMatrix(mVals(i), nVals(i));
                B = RandomMatrices.intMatrix(nVals(i), kVals(i));
                C = RandomMatrices.intMatrix(kVals(i), mVals(i));

                AB_tr = (A * B)';
                testCase.verifyEqual(trace(A * B * C), dot(AB_tr(:), C(:)));
            end 
        end
    
        %% WAW' = sum(sum(a_kl * w_k * w_l'))
        function testIdentity9(testCase)
            dVals = [5, 5, 2, 10];
            kVals = [5, 5, 12, 15];
            for i = 1:length(dVals)
                W = RandomMatrices.intMatrix(dVals(i), kVals(i));
                A = RandomMatrices.intMatrix(kVals(i), kVals(i));

                res = W * A * W';
                sum = 0;
                for k = 1:kVals(i)
                    for l = 1:kVals(i)
                        sum = sum + A(k, l) * W(:, k) * W(:, l)';
                    end
                end
                   
                testCase.verifyEqual(sum, res);
            end 
        end

        % FAILS!!!
        % %% sum(sum(a_kl * mu(k, :)' * mu(l, :))) = A x (mu * mu') where x is element size
        %% mu(k, :)' * mu(l, :) are elements of mu * mu'
        function testIdentity10(testCase)
            dVals = [5, 5, 2, 10];
            kVals = [5, 5, 12, 15];
            for i = 1:length(dVals)
                MU = RandomMatrices.intMatrix(kVals(i), dVals(i));
                
                res = MU * MU';

                for k = 1:kVals(i)
                    for l = 1:kVals(i)
                        testCase.verifyEqual(MU(k, :)' * MU(l, :), res(k, l));
                    end
                end
            end 
        end


        % For categorical
        function testIdentity11(testCase)
            dVals = [5, 5, 2, 10];
            kVals = [5, 5, 12, 15];
            for i = 1:length(dVals)
                W = RandomMatrices.intMatrix(dVals(i), kVals(i));
                A = RandomMatrices.intMatrix(dVals(i), dVals(i));
                
                W_t = W';
                res = W_t * A * W;

                sum = zeros(kVals(i), kVals(i));
                for d = 1:dVals(i)
                    for d2 = 1:dVals(i)
                        dCol = W_t(:, d);
                        dRow = W(d2, :);
                        sum = sum + A(d, d2) * dCol * dRow;
                    end
                end
                testCase.verifyEqual(res, sum);
            end 
        end

        function testIdentity12(testCase)
            dVals = [5, 5, 2, 10];
            kVals = [5, 5, 12, 15];
            nVals = [15, 25, 12, 15];
            for i = 1:length(dVals)
                Z = RandomMatrices.intMatrix(kVals(i), nVals(i));
                sum = zeros(kVals(i), kVals(i));

                c = 3;
                for n = 1:nVals(i)
                    sum = sum + Z(:, n) * c * Z(:, n)';
                end

                testCase.verifyEqual(sum, c * (Z * Z'));
            end 
        end

        function testIdentity13(testCase)
            QdVals = [5, 5, 2, 10];
            kVals = [5, 5, 12, 15];
            nVals = [15, 25, 12, 15];

            for i = 1:length(QdVals)
                W = RandomMatrices.intMatrix(QdVals(i), kVals(i));
                A = RandomMatrices.intMatrix(QdVals(i), QdVals(i));

                Z = RandomMatrices.intMatrix(kVals(i), nVals(i));

                sum1 = 0;

                for n = 1:nVals(i)
                    sum1 = sum1 + Z(:, n)' * W' * A * W * Z(:, n);
                end

                sum2 = 0;
                for q = 1:QdVals(i)
                    for q2 = 1:QdVals(i)
                        sum2 = sum2 + W(q, :) * A(q, q2) * (Z * Z') * W(q2, :)';
                    end
                end

                sum3 = 0;
                for q = 1:QdVals(i)
                    innerSum = zeros(kVals(i), 1);
                    for q2 = 1:QdVals(i)
                        innerSum = innerSum + A(q, q2) * W(q2, :)';
                    end

                    testCase.verifyEqual(innerSum, W' * A(q, :)');

                    sum3 = sum3 + W(q, :) * (Z * Z') * innerSum;
                end    

                testCase.verifyEqual(sum1, sum3);

                sum4 = 0;
                for q = 1:QdVals(i)
                    sum4 = sum4 + W(q, :) * (Z * Z') * W' * A(q, :)';
                end    

                testCase.verifyEqual(sum1, sum4);
            end 
        end


        function testIdentitytotal(testCase)
            QdVals = [5, 5, 2, 10];
            kVals = [5, 5, 12, 15];
            nVals = [15, 25, 12, 15];

            for i = 1:length(QdVals)
                W = RandomMatrices.intMatrix(QdVals(i), kVals(i));
                A = RandomMatrices.intMatrix(QdVals(i), QdVals(i));
                Z = RandomMatrices.intMatrix(kVals(i), nVals(i));
                Rd = RandomMatrices.intMatrix(QdVals(i), nVals(i));
                alpha = RandomMatrices.intMatrix(kVals(i), 1);

                sum1 = 0;

                for n = 1:nVals(i)
                    zn = Z(:, n);
                    sum1 = sum1 + zn'* W' * Rd(:, n) ...
                        - 1/2 * zn' * W' * A * W * zn;
                end
                for k = 1:kVals(i)
                    sum1 = sum1 - 1/2 * alpha(k) * W(:, k)' * W(:, k);
                end

              
                sum2 = 0;
                for q = 1:QdVals(i)
                    innerSum = zeros(kVals(i), 1);
                    for q2 = 1:QdVals(i)
                        if q2 ~= q
                            innerSum = innerSum + A(q, q2) * W(q2, :)';
                        end
                    end
                    sum2 = sum2 - 1/2 * W(q, :) * (Z * Z' * A(q, q) + diag(alpha))  * W(q, :)' + ...
                        W(q, :) * (Z * Rd(q, :)' - 1/2 * (Z * Z') * innerSum);
                end

                testCase.verifyEqual(sum1, sum2);
            end 
        end

        function testIdentityColsTerm1(testCase)
            QdVals = [5, 5, 2, 10];
            kVals = [5, 5, 12, 15];
            nVals = [15, 25, 12, 15];

            for i = 1:length(QdVals)
                W = RandomMatrices.intMatrix(QdVals(i), kVals(i));
                A = RandomMatrices.intMatrix(QdVals(i), QdVals(i));
                Z = RandomMatrices.intMatrix(kVals(i), nVals(i));
                Rd = RandomMatrices.intMatrix(QdVals(i), nVals(i));
                alpha = RandomMatrices.intMatrix(kVals(i), 1);

                sum1 = 0;

                for n = 1:nVals(i)
                    zn = Z(:, n);
                    sum1 = sum1 + zn'* W' * Rd(:, n);
                end
                
                sum2 = 0;
                for k = 1:kVals(i)
                    sum2 = sum2 + W(:, k)' * Z(k, :)' * Rd(d, :);
                end

                testCase.verifyEqual(sum1, sum2);
            end 
        end

        function testIdentityColsTerm2(testCase)
            QdVals = [5, 5, 2, 10];
            kVals = [5, 5, 12, 15];
            nVals = [15, 25, 12, 15];

            for i = 1:length(QdVals)
                W = RandomMatrices.intMatrix(QdVals(i), kVals(i));
                A = RandomMatrices.intMatrix(QdVals(i), QdVals(i));
                Z = RandomMatrices.intMatrix(kVals(i), nVals(i));
                Rd = RandomMatrices.intMatrix(QdVals(i), nVals(i));
                alpha = RandomMatrices.intMatrix(kVals(i), 1);

                sum1 = 0;

                for n = 1:nVals(i)
                    zn = Z(:, n);
                    sum1 = sum1 + zn'* W' * A * W * zn;
                end
                
                sum2 = 0;
                for k = 1:kVals(i)
                    for n = 1:nVals(i)
                        for j = 1:kVals(i)
                            sum2 = sum2 + Z(j, n) * W(:, j)' * A * W(:, k) * Z(k, n);
                        end
                    end
                end

                testCase.verifyEqual(sum1, sum2);
            end 
        end


        % Page 38 iz sveske!
        function testIdentityColsTerm21(testCase)
            QdVals = [5, 5, 2, 10];
            kVals = [5, 5, 12, 15];
            nVals = [15, 25, 12, 15];

            for i = 1:length(QdVals)
                W = RandomMatrices.intMatrix(QdVals(i), kVals(i));
                alpha = RandomMatrices.intMatrix(kVals(i), 1);
                
                sum1 = 0;
                for k = 1:kVals(i)
                    w_dk = W(:, k);
                    sum1 = sum1 + w_dk' * alpha(k) * w_dk;
                end

                vec = W(:);
                % Passes the test
                % sum2 = sum( sum( W.^2 .* alpha' ) );
                
                % This WORKS!
                alpha_diag = kron(alpha, ones(QdVals(i), 1));
                % sum2 = vec' * diag(alpha_diag) * vec;


                % THIS ALSO WORKS!!!
                sum2 = vec' * kron(diag(alpha), eye(QdVals(i))) * vec;

                testCase.verifyEqual(sum1, sum2);
            end 
        end
    end
end
