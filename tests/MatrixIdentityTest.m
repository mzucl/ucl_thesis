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
    end
end
