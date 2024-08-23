classdef MatrixIdentityTest < matlab.unittest.TestCase
    methods (Test)
        %% sum(an^T * an) = Tr(A^T * A)
        function testIdentity1(testCase)
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


        %% sum(an^T * v) = vT * A * 1 (vector of 1s)
        function testIdentity2(testCase)
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


        %% sum(an^T * W * bn) = Tr(WBA^T)
        function testIdentity3(testCase)
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


        %% sum(an^T * W * v) = 1^T * A^T * W * v
        % This one is a special case of one of the above
        function testIdentity4(testCase)
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


        %% sum(an * an^T) = AA^T
        function testIdentity5(testCase)
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

       
        %% trace(A^T * A) = sum of squared elements of A
        function testIdentity6(testCase)
            mVals = [5, 5, 2, 10];
            nVals = [5, 3, 7, 10];
            for i = 1:length(mVals)
                A = Utility.generateRandomIntMatrix(mVals(i), nVals(i));

                testCase.verifyEqual(trace(A' * A), dot(A(:), A(:)));
            end 
        end
        

        %% trace(A * B) = dot(A'(:), B(:))
        function testIdentity7(testCase)
            mVals = [5, 5, 2, 10];
            nVals = [5, 3, 7, 10];
            for i = 1:length(mVals)
                A = Utility.generateRandomIntMatrix(mVals(i), nVals(i));
                B = Utility.generateRandomIntMatrix(nVals(i), mVals(i));

                A_tr = A';
                testCase.verifyEqual(trace(A * B), dot(A_tr(:), B(:)));
            end 
        end


        %% trace(A * B * C) = dot(AB'(:), C(:))
            % this is the same as trace(CAB) and trace(BCA), so we can choose 
            % which two matrices to multiply
        function testIdentity8(testCase)
            
            mVals = [5, 5, 2, 10];
            nVals = [5, 3, 7, 10];
            kVals = [5, 5, 12, 15];
            for i = 1:length(mVals)
                A = Utility.generateRandomIntMatrix(mVals(i), nVals(i));
                B = Utility.generateRandomIntMatrix(nVals(i), kVals(i));
                C = Utility.generateRandomIntMatrix(kVals(i), mVals(i));

                AB_tr = (A * B)';
                testCase.verifyEqual(trace(A * B * C), dot(AB_tr(:), C(:)));
            end 
        end
    end
end
