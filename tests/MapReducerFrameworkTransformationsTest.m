classdef MapReducerFrameworkTransformationsTest < matlab.unittest.TestCase
    methods (Test)
        %% Z * Z^T = ?
        % This method verifies that the full matrix product Z * Z^T can be reconstructed 
        % by summing Zc * Zc^T over all chunks Zc, where each Zc contains the latent 
        % variables corresponding to a subset of observations.
        function testIdentity1(testCase)
            % Setup
            K = 10;
            N = 50;
            chunkSize = 10;

            Z = Utility.generateRandomIntMatrix(K, N);
            
            numChunks = ceil(N / chunkSize);
            Zchunks = cell(1, numChunks);
            
            % Z = [Zchunks{1} | Zchunks{2} | ... | Zchunks{numChunks}]
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

        % This method performs the same verification as the one above, but allows 
        % for variable chunk sizes.
        function testIdentity1b(testCase)
            % Setup
            K = 10;
            N = 50;
            chunkSize = 2;

            Z = Utility.generateRandomIntMatrix(K, N);
            
            Zchunks = {};
            i = 1;
            startIdx = 1;
            
            while true
                endIdx = min(startIdx + chunkSize - 1, N);
                Zchunks{i} = Z(:, startIdx:endIdx);
            
                if endIdx == N
                    break;
                end
            
                startIdx = endIdx + 1;
                chunkSize = chunkSize * 2;  % Double the chunk size
                i = i + 1;
            end
            numChunks = i;

            sum = 0;
            for i = 1:numChunks
                sum = sum + Zchunks{i} * Zchunks{i}';
            end

            testCase.verifyEqual(Z * Z', sum);
        end
        
        %% Z * (X^T - mu^T) = ?
        % This method verifies that Z * (X^T - mu^T) can be reconstructed 
        % by summing Zc * (Xc^T - mu^T) over all chunks of latent
        % variables Zc and corresponding observations Xc.
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


        % TODO: comment
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

        % TODO: comment
        function testIdentity3b(testCase)
            % Setup
            D = 20;
            K = 10;
            N = 50;
            chunkSize = 10;

            Z = Utility.generateRandomIntMatrix(K, N);
            W = Utility.generateRandomIntMatrix(D, K);

            numChunks = ceil(N / chunkSize);
            
            sum = 0;
            for i = 1:numChunks
                startIdx = (i-1) * chunkSize + 1;
                endIdx = min(i * chunkSize, N);
                WZchunk = W * Z(:, startIdx:endIdx);
                
                Nc = size(WZchunk, 2);
                sum = sum + WZchunk * ones(Nc, 1);
            end
            
            testCase.verifyEqual(W * Z * ones(N, 1), sum);
        end

        %% Tr(X^T * X) = ?
        % This method verifies that Tr(X^T * X) can be reconstructed by 
        % summing Tr(Xc^T * Xc) over all chunks of observations.
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

        %% Tr(W * Z * X^T) = ?
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
