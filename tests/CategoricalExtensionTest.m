classdef CategoricalExtensionTest < matlab.unittest.TestCase
    % [NOTE] If multiple parameters are defined, the test method will be run for 
    % all combinations (a Cartesian product).
    properties (TestParameter)
        D = {10, 20};
        K = {8};
        N = {50, 100};
        Qd = {5, 7, 10};
    end
    
    methods (Test, TestTags = {'ELBO', 'vectorization'})
        % [NOTE] term 3 in the 'const' term
        function test_transformation_01(testCase, Qd, N)
            xi_d = RandomMatrices.intMatrix(Qd, N);
            A = RandomMatrices.intMatrix(Qd, Qd);
    
            scalarSum = 0;
            for n = 1:N
                xi_dn = xi_d(:, n);
                scalarSum = scalarSum + xi_dn' * A * xi_dn;
            end
            vectorizedSum = trace(xi_d' * A * xi_d);
    
            testCase.verifyEqual(scalarSum, vectorizedSum);
        end

        % [NOTE] term 2 in the 'const' term
        function test_transformation_02(testCase, Qd, N)
            xi_d = RandomMatrices.intMatrix(Qd, N);
            G = RandomMatrices.intMatrix(Qd, N);
    
            scalarSum = 0;
            for n = 1:N
                xi_dn = xi_d(:, n);
                g_dn = G(:, n);
                scalarSum = scalarSum + xi_dn' * g_dn;
            end
            vectorizedSum = trace(xi_d' * G);
    
            testCase.verifyEqual(scalarSum, vectorizedSum);
        end

        % [NOTE] term 1 in the 'non-const' term
        function test_transformation_03(testCase,K, Qd, N)
            X_d = RandomMatrices.intMatrix(Qd, N);
            T_d = RandomMatrices.intMatrix(Qd, N);
            mu = RandomMatrices.intMatrix(Qd, 1);
            Z = RandomMatrices.intMatrix(K, N);
            W = RandomMatrices.intMatrix(Qd, K);
    
            scalarSum = 0;
            for n = 1:N
                x_dn = X_d(:, n);
                t_dn = T_d(:, n);
                z_n = Z(:, n);
                scalarSum = scalarSum + (z_n' * W' + mu') * (x_dn + t_dn);
            end
            % vectorizedSum = trace((Z' * W' + mu') * (X_d + T_d));
            vectorizedSum = trace((W * Z + mu)' * (X_d + T_d));
    
            testCase.verifyEqual(scalarSum, vectorizedSum);
        end

        function test_transformation_04(testCase,K, Qd, N)
            A = RandomMatrices.intMatrix(Qd, Qd);
            Z = RandomMatrices.intMatrix(K, N);
            W = RandomMatrices.intMatrix(Qd, K);
    
            scalarSum = 0;
            for n = 1:N
                z_n = Z(:, n);
                scalarSum = scalarSum + z_n' * W' * A * W * z_n;
            end
            vectorizedSum = trace((Z * Z') * W' * A * W); 
    
            testCase.verifyEqual(scalarSum, vectorizedSum);
        end

        % When I thought that the expectation <W'AW> can be computed using
        % vec(W) representation
         function test_transformation_05(testCase,K, Qd)
            A = RandomMatrices.intMatrix(Qd, Qd);
            W = RandomMatrices.intMatrix(Qd, K);
    
            scalarSum = 0;
            for k = 1:K
                w_k = W(:, k);
                scalarSum = scalarSum +  w_k' * A * w_k;
            end
            vectorizedSum = trace(W' * A * W);
    
            testCase.verifyEqual(scalarSum, vectorizedSum);
         end

          function test_transformation_06(testCase,K, Qd, N)
            A = RandomMatrices.intMatrix(Qd, Qd);
            W = RandomMatrices.intMatrix(Qd, K);
            mu = RandomMatrices.intMatrix(Qd, 1);
            Z = RandomMatrices.intMatrix(K, N);
    
            scalarSum = 0;
            for n = 1:N
                z_n = Z(:, n);
                scalarSum = scalarSum + z_n' * W' * A * mu; 
            end
            vectorizedSum = ones(1, N) * Z' * W' * A * mu;
    
            testCase.verifyEqual(scalarSum, vectorizedSum);
        end
    end
end
