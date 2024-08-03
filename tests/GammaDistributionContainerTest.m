classdef GammaDistributionContainerTest < matlab.unittest.TestCase
    methods (Test)
        %% Skipped this for now, it is a very long code with many cases
        % function testConstructor(testCase)
        % 
        % end
        %%

        function testDependentProperties(testCase)
            % Test 1
            a = 1; b = 2; size = 10;
            obj = GammaDistributionContainer(a, b, size);

            testCase.verifyEqual(obj.Size, size);
            
            for i = 1:size
                testCase.verifyEqual(obj.Expectation{i}, a / b);
            end

            % Test 2
            size = 5;
            bVals = ones(1, size);
            aVals = 1:size;
            
            obj = GammaDistributionContainer(aVals, bVals);
            for i = 1:size
                testCase.verifyEqual(obj.Expectation{i}, aVals(i) / bVals(i));
            end
            
            a = 1; b = 1;
            obj = GammaDistributionContainer(a, b, size);
            for i = 1:size
                testCase.verifyEqual(obj.H{i}, 1);
            end
            testCase.verifyEqual(obj.HC, obj.Size);


        end

        function testSingleDistributionMethods(testCase)
            bVals = ones(1, 5);
            aVals = 1:5;
            
            obj = GammaDistributionContainer(aVals, bVals);

            % 'getDistribution' method
            idx = 3;
            testCase.verifyEqual(obj.getDistribution(idx).a, idx);

            % 'updateDistributionParams' method
            deltaA = 0.1; deltaB = 0.01;
            obj.updateDistributionParams(idx, deltaA, deltaB);

            testCase.verifyEqual(obj.distributions(idx).a, aVals(idx) + deltaA);
            testCase.verifyEqual(obj.distributions(idx).b, bVals(idx) + deltaB);

            % 'updateDistributionParams' method with 'inc' = false
            obj.updateDistributionParams(idx, deltaA, deltaB, false);

            testCase.verifyEqual(obj.distributions(idx).a, deltaA);
            testCase.verifyEqual(obj.distributions(idx).b, deltaB);
            
            obj.updateAllDistributionsParams(aVals, bVals, false);

            for i = 1:obj.Size
                testCase.verifyEqual(obj.distributions(i).a, aVals(i));
                testCase.verifyEqual(obj.distributions(i).b, bVals(i));
            end
        end
    end
end
