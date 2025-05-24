classdef Copy_of_MatrixTransformationsTest < matlab.unittest.TestCase
    methods (Test)
        % TODO: write a comment
        function testIdentity1(testCase)
            % Setup
            K = 10;
            N = 50;
            chunkSize = 10;

            Z = Utility.generateRandomIntMatrix(K, N);
            
            numChunks = ceil(N / chunkSize);
            Zchunks = cell(1, numChunks);
            
            for i = 1:numChunks
                startIdx = (i-1) * chunkSize + 1;
                endIdx = min(i * chunkSize, N);
                Zchunks{i} = Z(:, startIdx:endIdx);
            end

            sum = 0;
            for i = 1:numChunks
                sum = sum + Zchunks{i} * Zchunks{i}';
            end

            testCase.verifyEqual(Z * Z', sum);
        end
        
        % TODO: write a comment: 
        function testIdentity2(testCase)
            % Setup
            D = 20;
            K = 10;
            N = 50;
            chunkSize = 10;

            X = Utility.generateRandomIntMatrix(D, N);
            Z = Utility.generateRandomIntMatrix(K, N);
            mu = Utility.generateRandomIntMatrix(D, 1);

            numChunks = ceil(N / chunkSize);
            Zchunks = cell(1, numChunks);
            Xchunks = cell(1, numChunks);
            
            for i = 1:numChunks
                startIdx = (i-1) * chunkSize + 1;
                endIdx = min(i * chunkSize, N);
                Zchunks{i} = Z(:, startIdx:endIdx);
                Xchunks{i} = X(:, startIdx:endIdx);
            end
            
            sum = 0;
            for i = 1:numChunks
                sum = sum + Zchunks{i} * (Xchunks{i}' - mu');
            end

            testCase.verifyEqual(Z * (X' - mu'), sum);
        end

        % W * [Z1 | Z2 | ... Zs] = [WZ1 | WZ2...
        function testIdentity3(testCase)
            % Setup
            D = 20;
            K = 10;
            N = 50;
            chunkSize = 10;

            Z = Utility.generateRandomIntMatrix(K, N);
            W = Utility.generateRandomIntMatrix(D, K);

            numChunks = ceil(N / chunkSize);
            WZchunks = cell(1, numChunks);
            
            for i = 1:numChunks
                startIdx = (i-1) * chunkSize + 1;
                endIdx = min(i * chunkSize, N);
                WZchunks{i} = W * Z(:, startIdx:endIdx);
            end
            
            testCase.verifyEqual(W * Z, [WZchunks{:}]);
        end

        % Tr(X' * X) = sum (Tr(X_i' * X_i))
        function testIdentity4(testCase)
            % Setup
            D = 20;
            N = 50;
            chunkSize = 10;

            X = Utility.generateRandomIntMatrix(D, N);

            numChunks = ceil(N / chunkSize);
            Xchunks = cell(1, numChunks);
            
            for i = 1:numChunks
                startIdx = (i-1) * chunkSize + 1;
                endIdx = min(i * chunkSize, N);
                Xchunks{i} = X(:, startIdx:endIdx);
            end

            sum = 0;
            for i = 1:numChunks
                sum = sum + trace(Xchunks{i}' * Xchunks{i});
            end
            
            testCase.verifyEqual(trace(X' * X), sum);
        end

        function testIdentity5(testCase)
            % Setup
            K = 10;
            D = 20;
            N = 50;
            chunkSize = 10;

            X = Utility.generateRandomIntMatrix(D, N);
            Z = Utility.generateRandomIntMatrix(K, N);
            W = Utility.generateRandomIntMatrix(D, K);

            numChunks = ceil(N / chunkSize);
            Xchunks = cell(1, numChunks);
            Zchunks = cell(1, numChunks);
            
            for i = 1:numChunks
                startIdx = (i-1) * chunkSize + 1;
                endIdx = min(i * chunkSize, N);
                Xchunks{i} = X(:, startIdx:endIdx);
                Zchunks{i} = Z(:, startIdx:endIdx);
            end

            sum = 0;
            for i = 1:numChunks
                sum = sum + trace(W * Zchunks{i} * Xchunks{i}');
            end
            
            testCase.verifyEqual(trace(W * Z * X'), sum);
        end
    end
end
