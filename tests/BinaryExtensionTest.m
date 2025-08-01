classdef BinaryExtensionTest < matlab.unittest.TestCase
    % [NOTE] If multiple parameters are defined, the test method will be run for 
    % all combinations (a Cartesian product).
    properties (TestParameter)
        D = {10, 20};
        K = {8};
        N = {50, 100};
    end


    
    % [NOTE] The vectorized form a' * diag(h) * a is used in the derivation steps, 
    % as it is easier to manipulate algebraically using standard matrix 
    % multiplication rules. However, in the implementation, the equivalent 
    % expression sum(h .* a.^2) is used instead, as it is more efficient and 
    % avoids the allocation of the full D-by-D diagonal matrix.
    methods (Test, TestTags = {'qZUpdate', 'transformation'})
        function test_transformation_01(testCase, D)
            %test_transformation_01
            % Verifies that the summation sum_d h_d * a_d^2 matches 
            % the vectorized summation sum(h .* a.^2).
    
            h = Utility.generateRandomIntMatrix(D, 1);
            a = Utility.generateRandomIntMatrix(D, 1);
    
            scalarSum = 0;
            for d = 1:D
                scalarSum = scalarSum + h(d) * a(d)^2;
            end
            vectorizedSum = sum(h .* a.^2);
    
            testCase.verifyEqual(scalarSum, vectorizedSum);
        end
    end


    methods (Test, TestTags = {'qZUpdate', 'vectorization'})
        function test_qZUpdate_01(testCase, D)
            %test_qZUpdate_01
            % Verifies that the manual summation sum_d h_d * a_d^2 matches 
            % the vectorized form a' * diag(h) * a.
 
            h = Utility.generateRandomIntMatrix(D, 1);
            a = Utility.generateRandomIntMatrix(D, 1);
    
            scalarSum = 0;
            for d = 1:D
                scalarSum = scalarSum + h(d) * a(d)^2;
            end
    
            vectorizedForm = a' * diag(h) * a;
            testCase.verifyEqual(scalarSum, vectorizedForm);
        end

        function test_qZUpdate_02(testCase, D, N)
            %test_qZUpdate_02
            % Verifies that H .* mu is equivalent to computing diag(H(:,n)) * mu 
            % column-wise for all n.
    
            mu = Utility.generateRandomIntMatrix(D, 1);
            H = Utility.generateRandomIntMatrix(D, N);
    
            % Vectorized computation: broadcasting mu across columns
            vecRes = H .* mu;
    
            for n = 1:N
                singleColRes = diag(H(:, n)) * mu;
                testCase.verifyEqual(singleColRes, vecRes(:, n));
            end
        end
    end


    methods (Test, TestTags = {'qMuUpdate', 'vectorization'})
        function test_qMuUpdate_01(testCase, D, N)
            %test_qMuUpdate_01
            % Verifies that diag(H * ones(N, 1)) is equivalent to summing diag(H(:,n))
            % over all columns n.
        
            H = Utility.generateRandomIntMatrix(D, N);
        
            diagSumVectorized = diag(H * ones(N, 1));
        
            diagSumManual = zeros(D, D);
            for n = 1:N
                diagSumManual = diagSumManual + diag(H(:, n));
            end
        
            testCase.verifyEqual(diagSumVectorized, diagSumManual);
        end

        function test_qMuUpdate_02(testCase, D, K, N)
            %test_qMuUpdate_02
            % Verifies that the vectorized expression H .* (W * Z) * ones(N, 1) 
            % is equivalent to the manual summation over n of diag(H(:, n)) * W * Z(:, n).
        
            H = Utility.generateRandomIntMatrix(D, N);
            W = Utility.generateRandomIntMatrix(D, K);
            Z = Utility.generateRandomIntMatrix(K, N);
        
            vecSum = H .* (W * Z) * ones(N, 1);
        
            manualSum = zeros(D, 1);
            for n = 1:N
                manualSum = manualSum + diag(H(:, n)) * W * Z(:, n);
            end
        
            testCase.verifyEqual(manualSum, vecSum);
        end


        function test_qMuUpdate_03(testCase, D, N)
            %test_qMuUpdate_03
            % Verifies that when H is a constant matrix with all entries 1/4 (Bohning bound),
            % diag(H * ones(N, 1)) equals (N/4) times the identity matrix.
        
            H = 1/4 * ones(D, N);
        
            expected = (N / 4) * eye(D, D);
            actual = diag(H * ones(N, 1));
        
            testCase.verifyEqual(actual, expected);
        end
    end


    methods (Test, TestTags = {'qWUpdate', 'vectorization'})
        function test_qWUpdate_01(testCase, D, K, N)
            %test_qWUpdate_01
            % Verifies transformation of the term with z_n^T:
            % sum z_n^T * W^T * r_n over n is equivalent to summing
            % W(d,:) * Z * r_d^T over d, where r_n is a column of R and r_d is a row of R.
        
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

        function test_qWUpdate_02(testCase, D, K, N)
            %test_qWUpdate_02
            % Verifies transformation of the term with z_n^T * z_n:
            % sum z_n^T * W^T * diag(h_n) * W * z_n over n is 
            % equivalent to summing W(d,:) * Z * diag(h_d) * Z' * W(d,:)' over d, 
            % where h_n is a column of H and h_d is a row of H.
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

        function test_qWUpdate_03(testCase, D, K, N)
            %test_qWUpdate_03
            % Verifies that the d-th column of sigma * Z * R' equals sigma * Z * r_d',
            % where r_d is the d-th row of R.
                
            sigma = Utility.generateRandomIntMatrix(K, K);
            Z = Utility.generateRandomIntMatrix(K, N);
            R = Utility.generateRandomIntMatrix(D, N);
        
            vecRes = sigma * Z * R';
        
            for d = 1:D
                singleColRes = sigma * Z * R(d, :)';
                testCase.verifyEqual(singleColRes, vecRes(:, d));
            end
        end
    end


    methods (Test, TestTags = {'ELBO', 'vectorization'})
        function test_ELBO_01(testCase, D, K, N)
            %test_ELBO_01
            % Verifies vectorization of the phi_1 term in ELBO: 
            % sum over n of (z_n' * W' + mu') * (x_n + t_n) 
            % equals the trace of (Z' * W' + mu') * (X + T), 
            % where z_n, x_n, and t_n are the n-th columns of Z, X, and T, respectively.
            W = Utility.generateRandomIntMatrix(D, K);
            Z = Utility.generateRandomIntMatrix(K, N);
            T = Utility.generateRandomIntMatrix(D, N);
            X = Utility.generateRandomIntMatrix(D, N);
            mu = Utility.generateRandomIntMatrix(D, 1);

            scalarSum = 0;
            for n = 1:N
                scalarSum = scalarSum + (Z(:, n)' * W'  + mu') * (X(:, n) + T(:, n));
            end

            % Vectorized form
            vecRes = trace((Z' * W' + mu') * (X + T));

            testCase.verifyEqual(scalarSum, vecRes);
        end

        function test_ELBO_02(testCase, D, K, N)
            %test_ELBO_02
            % Verifies vectorization of the phi_2 linear term in ELBO:
            % sum over n of z_n' * W' * diag(h_n) * mu
            % equals ones(1, N) * ((Z' * W') .* H') * mu,
            % where z_n and h_n are the n-th columns of Z and H.
        
            W = Utility.generateRandomIntMatrix(D, K);
            Z = Utility.generateRandomIntMatrix(K, N);
            H = Utility.generateRandomIntMatrix(D, N);
            mu = Utility.generateRandomIntMatrix(D, 1);
        
            scalarSum = 0;
            for n = 1:N
                scalarSum = scalarSum + Z(:, n)' * W' * diag(H(:, n)) * mu;
            end
        
            % Vectorized form
            vecRes = ones(1, N) * ((Z' * W') .* H') * mu;
        
            testCase.verifyEqual(scalarSum, vecRes);
        end

        function test_ELBO_03(testCase, D, K, N)
            %test_ELBO_03
            % Verifies vectorization of the phi_2 square term in ELBO:
            % sum over n of z_n' * W' * diag(h_n) * W * z_n
            % equals trace of (W * Z)' * (H .* (W * Z)),
            % where z_n and h_n are the n-th columns of Z and H.
        
            W = Utility.generateRandomIntMatrix(D, K);
            Z = Utility.generateRandomIntMatrix(K, N);
            H = Utility.generateRandomIntMatrix(D, N);
        
            scalarSum = 0;
            for n = 1:N
                scalarSum = scalarSum + Z(:, n)' * W' * diag(H(:, n)) * W * Z(:, n);
            end
        
            % Vectorized form
            vecRes = trace((W * Z)' * (H .* (W * Z)));
        
            testCase.verifyEqual(scalarSum, vecRes);
        end

        function test_ELBO_04(testCase, D, N)
            %test_ELBO_04
            % Verifies vectorization of the constant term in ELBO:
            % sum over n,d of -c_{d,n} + xi_{d,n} * g_{d,n} - 1/2 * xi_{d,n}^2 * h_{d,n}
            % equals sum of vectorized expression -C + Xi .* G - 1/2 * Xi.^2 .* H.
        
            Xi = Utility.generateRandomIntMatrix(D, N);
            C = Utility.generateRandomIntMatrix(D, N);
            G = Utility.generateRandomIntMatrix(D, N);
            H = Utility.generateRandomIntMatrix(D, N);
        
            scalarSum = 0;
            for n = 1:N
                for d = 1:D
                    scalarSum = scalarSum - C(d, n) + Xi(d, n) * G(d, n) - 1/2 * Xi(d, n)^2 * H(d, n);
                end
            end
        
            vecRes = sum(sum(-C + Xi .* G - 1/2 * Xi.^2 .* H));
        
            testCase.verifyEqual(scalarSum, vecRes);
        end
    end
end
