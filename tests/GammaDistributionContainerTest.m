% TODO(medium/high) If any of the priors passed in 'priors' array is NaN
% then some of tests fail. If this (some of distributions have prior and
% some doesn't) is needed as a feature fix it!

classdef GammaDistributionContainerTest < matlab.unittest.TestCase
    methods (Test)
        %% Constructor tests
        function testDefaultConstructor(testCase)
            obj = GammaDistributionContainer();
            testCase.verifyEqual(obj.Size, 1); % Single element container

            GammaDistributionTest.verifyObject(testCase, obj.distributions(1), Constants.DEFAULT_GAMMA_A, ...
                Constants.DEFAULT_GAMMA_B, NaN);
        end

        function testConstructorWithOneParameter(testCase)
            % Test 1: 'a' is a scalar
            dim = 5;
            obj = GammaDistributionContainer(dim);

            testCase.verifyEqual(obj.Size, 5);
            for i = 1:dim
                GammaDistributionTest.verifyObject(testCase, obj.distributions(i), Constants.DEFAULT_GAMMA_A, ...
                    Constants.DEFAULT_GAMMA_B, NaN);
            end

            % Test 2: 'a' is single GammaDistribution
            truePrior = GammaDistribution();
            obj = GammaDistributionContainer(truePrior);

            testCase.verifyEqual(obj.Size, 1);
            GammaDistributionTest.verifyObject(testCase, obj.distributions(1), Constants.DEFAULT_GAMMA_A, ...
                Constants.DEFAULT_GAMMA_B, truePrior);


            % Test 3: 'a' is an array of GammaDistribution objects
            truePriors = [GammaDistribution(), GammaDistribution(1, 2), GammaDistribution(4, 3)];
            obj = GammaDistributionContainer(truePriors);
            numDists = length(truePriors);

            testCase.verifyEqual(obj.Size, numDists);
            for i = 1:numDists
                GammaDistributionTest.verifyObject(testCase, obj.distributions(i), truePriors(i).a, ...
                truePriors(i).b, truePriors(i));
            end
        end

        function testConstructorWithTwoParameters(testCase)
            % Test 1: 
            % a: scalar
            % b: scalar
            a = 1; b = 2;
            obj = GammaDistributionContainer(a, b);

            testCase.verifyEqual(obj.Size, 1); % Single element container
            for i = 1:obj.Size
                testCase.verifyTrue(isnan(obj.distributions(i).prior), 'All priors are not NaN');
            end

            % Test 2: 
            % a: scalar
            % b: array
            a = 1; b = [1, 2, 3, 4];
            obj = GammaDistributionContainer(a, b);

            testCase.verifyEqual(obj.Size, length(b));
            for i = 1:obj.Size
                GammaDistributionTest.verifyObject(testCase, obj.distributions(i), a, b(i), NaN);
            end

            % Test 3: 
            % a: array
            % b: scalar
            a = 1:10; b = 2;
            obj = GammaDistributionContainer(a, b);

            testCase.verifyEqual(obj.Size, length(a));
            for i = 1:obj.Size
                GammaDistributionTest.verifyObject(testCase, obj.distributions(i), a(i), b, NaN);
            end

            % Test 4: 
            % a: array
            % b: array
            a = 1:10; b = ones(10, 1);
            obj = GammaDistributionContainer(a, b);

            testCase.verifyEqual(obj.Size, length(a));
            for i = 1:obj.Size
                GammaDistributionTest.verifyObject(testCase, obj.distributions(i), a(i), b(i), NaN);
            end
        end

        function testConstructorWithThreeParameters(testCase)
            aPrior = 4; bPrior = 5;
            priors = GammaDistribution(aPrior, bPrior);

            % Test 1: 
            % a: scalar
            % b: scalar
            a = 1; b = 2;
            
            obj = GammaDistributionContainer(a, b, priors);

            testCase.verifyEqual(obj.Size, 1); % Single element container
            for i = 1:obj.Size
                GammaDistributionTest.verifyObject(testCase, obj.distributions(i), a, b, priors);
                GammaDistributionTest.verifyObject(testCase, obj.distributions(i).prior, aPrior, bPrior);
            end

            % Test 2: 
            % a: scalar
            % b: array
            a = 1; b = [1, 2, 3, 4];
            obj = GammaDistributionContainer(a, b, priors);

            testCase.verifyEqual(obj.Size, length(b));
            for i = 1:obj.Size
                GammaDistributionTest.verifyObject(testCase, obj.distributions(i), a, b(i), priors);
                GammaDistributionTest.verifyObject(testCase, obj.distributions(i).prior, aPrior, bPrior);
            end

            % Test 3: 
            % a: array
            % b: scalar
            a = 1:10; b = 2;
            obj = GammaDistributionContainer(a, b, priors);

            testCase.verifyEqual(obj.Size, length(a));
            for i = 1:obj.Size
                GammaDistributionTest.verifyObject(testCase, obj.distributions(i), a(i), b, priors);
                GammaDistributionTest.verifyObject(testCase, obj.distributions(i).prior, aPrior, bPrior);
            end

            % Test 4: 
            % a: array
            % b: array
            a = 1:10; b = ones(10, 1);
            obj = GammaDistributionContainer(a, b, priors);

            testCase.verifyEqual(obj.Size, length(a));
            for i = 1:obj.Size
                GammaDistributionTest.verifyObject(testCase, obj.distributions(i), a(i), b(i), priors);
                GammaDistributionTest.verifyObject(testCase, obj.distributions(i).prior, aPrior, bPrior);
            end
            
            % Test 5: 
            % a: scalar
            % b: array
            a = 1; b = [1, 2, 3, 4];
            priorVal = GammaDistribution(aPrior, bPrior);
            priors = [priorVal.copy(), priorVal.copy(), priorVal.copy(), priorVal.copy()];
            obj = GammaDistributionContainer(a, b, priors);

            testCase.verifyEqual(obj.Size, length(b));
            for i = 1:obj.Size
                GammaDistributionTest.verifyObject(testCase, obj.distributions(i), a, b(i), priors(i));
            end
        end

        function testConstructorWithFourParameters(testCase)
            % Test 1: 'priors' set to NaN
            a = 1; b = 2; size = 10;
            obj = GammaDistributionContainer(a, b, NaN, size);

            testCase.verifyEqual(obj.Size, size);
            for i = 1:size
                GammaDistributionTest.verifyObject(testCase, obj.distributions(i), a, b, NaN);
            end

            % Test 2: 'priors' is set to a single GammaDistribution
            aPrior = 5;
            bPrior = 2;
            priors = GammaDistribution(aPrior, bPrior);
            obj = GammaDistributionContainer(a, b, priors, 10);

            testCase.verifyEqual(obj.Size, size);
            for i = 1:obj.Size
                testCase.verifyTrue(obj.distributions(i).prior == priors);
            end

            % Test 3: 'priors' is set to an array of GammaDistribution
            % objects
            aPrior = 5;
            bPrior = 2;
            priors = [GammaDistribution(aPrior, bPrior), GammaDistribution(aPrior, bPrior)];
            obj = GammaDistributionContainer(a, b, priors, length(priors));

            testCase.verifyEqual(obj.Size, length(priors));
            for i = 1:obj.Size
                testCase.verifyTrue(obj.distributions(i).prior == priors(i));
            end
        end



        %% Private properties
        function testSetters(testCase)
            a = 1; b = 2; size = 10;
            obj = GammaDistributionContainer(a, b, NaN, size);

            testCase.verifyTrue(all(isequal(obj.getExpCInit(), obj.ExpectationC)));
            
            % randomArray = minValue + (maxValue - minValue) * rand(1, n);
            newExpC = 0.1 + (5 - 0.1) * rand(1, obj.Size);
            obj.setExpCInit(newExpC);

            testCase.verifyTrue(isequal(obj.getExpCInit(), newExpC));
        end



        %% Dependent properties
        function testDependentProperties(testCase)
            % Test 1.1: Expectation, ExpectationC
            a = 1; b = 2; size = 10;
            obj = GammaDistributionContainer(a, b, NaN, size);

            testCase.verifyEqual(obj.Size, size);
            for i = 1:size
                testCase.verifyEqual(obj.Expectation{i}, a / b);
            end
            testCase.verifyEqual(obj.ExpectationC, a/b * ones(obj.Size, 1));

            % Test 1.2: Expectation, ExpectationC
            size = 5;
            bVals = ones(1, size);
            aVals = 1:size;
            
            obj = GammaDistributionContainer(aVals, bVals);
            for i = 1:size
                testCase.verifyEqual(obj.Expectation{i}, aVals(i) / bVals(i));
            end
            
            % Test 2: H, HC
            a = 1; b = 1;
            obj = GammaDistributionContainer(a, b, NaN, size);
            for i = 1:size
                testCase.verifyEqual(obj.H{i}, 1);
            end
            testCase.verifyEqual(obj.HC, obj.Size);
        end



        %% Single distribution and update methods
        function testSingleDistributionMethods(testCase)
            bVals = ones(1, 5);
            aVals = 1:5;
            
            obj = GammaDistributionContainer(aVals, bVals);

            % 'getDistribution' method
            idx = 3;
            testCase.verifyEqual(obj.getDistribution(idx).a, idx);

            % 'updateDistributionParams' method
            deltaA = 0.1; deltaB = 0.01;
            obj.updateDistributionParams(idx, deltaA, deltaB, true);

            testCase.verifyEqual(obj.distributions(idx).a, aVals(idx) + deltaA);
            testCase.verifyEqual(obj.distributions(idx).b, bVals(idx) + deltaB);

            % 'updateDistributionParams' method with 'inc' = false
            obj.updateDistributionParams(idx, deltaA, deltaB);

            testCase.verifyEqual(obj.distributions(idx).a, deltaA);
            testCase.verifyEqual(obj.distributions(idx).b, deltaB);
            
            obj.updateAllDistributionsParams(aVals, bVals);

            for i = 1:obj.Size
                testCase.verifyEqual(obj.distributions(i).a, aVals(i));
                testCase.verifyEqual(obj.distributions(i).b, bVals(i));
            end
        end

        function testUpdateAllDistributionsA(testCase)
            % Test 1: 
            % deltaA: scalar
            % inc: false
            numOfDistr = 5;
            aVals = 1:numOfDistr;
            bVals = ones(1, numOfDistr);
            deltaA = 0.1;

            obj = GammaDistributionContainer(aVals, bVals);
            obj.updateAllDistributionsA(deltaA);

            % Add test for priors - update shouldn't affect priors
            for i = 1:obj.Size
                testCase.verifyTrue(isnan(obj.distributions(i).prior), 'All priors are not NaN');
                testCase.verifyEqual(obj.distributions(i).a, deltaA);
                testCase.verifyEqual(obj.distributions(i).b, bVals(i)); % Not changed
            end

            % Test 2: 
            % deltaA: array
            % inc: true
            obj = GammaDistributionContainer(aVals, bVals);
            deltaA = 0.25 * 1:numOfDistr;
            obj.updateAllDistributionsA(deltaA, true);

            % Add test for priors - update shouldn't affect priors
            for i = 1:obj.Size
                testCase.verifyTrue(isnan(obj.distributions(i).prior), 'All priors are not NaN');
                testCase.verifyEqual(obj.distributions(i).a, aVals(i) + deltaA(i));
                testCase.verifyEqual(obj.distributions(i).b, bVals(i)); % Not changed
            end
        end

        function testUpdateAllDistributionsParams(testCase)
            % Test 1: 
            % deltaA: scalar
            % deltaB: scalar
            % inc: true

            numOfDistr = 5;
            aVals = 1:numOfDistr;
            bVals = ones(1, numOfDistr);
          
            obj = GammaDistributionContainer(aVals, bVals);
    
            deltaA = 0.1; deltaB = 0.01;
            obj.updateAllDistributionsParams(deltaA, deltaB, true);

            % Add test for priors - update shouldn't affect priors
            for i = 1:obj.Size
                testCase.verifyTrue(isnan(obj.distributions(i).prior), 'All priors are not NaN');
                testCase.verifyEqual(obj.distributions(i).a, aVals(i) + deltaA);
                testCase.verifyEqual(obj.distributions(i).b, bVals(i) + deltaB);
            end

            % Test 2: 
            % deltaA: scalar
            % deltaB: scalar
            % inc: false
            obj = GammaDistributionContainer(aVals, bVals);
            obj.updateAllDistributionsParams(deltaA, deltaB);

            % Add test for priors - update shouldn't affect priors
            for i = 1:obj.Size
                testCase.verifyTrue(isnan(obj.distributions(i).prior), 'All priors are not NaN');
                testCase.verifyEqual(obj.distributions(i).a, deltaA);
                testCase.verifyEqual(obj.distributions(i).b, deltaB);
            end

            % Test 3: 
            % deltaA: array
            % deltaB: scalar
            % inc: true
            obj = GammaDistributionContainer(aVals, bVals);

            deltaA = 0.25 * 1:numOfDistr;
            obj.updateAllDistributionsParams(deltaA, deltaB, true);

            % Add test for priors - update shouldn't affect priors
            for i = 1:obj.Size
                testCase.verifyTrue(isnan(obj.distributions(i).prior), 'All priors are not NaN');
                testCase.verifyEqual(obj.distributions(i).a, aVals(i) + deltaA(i));
                testCase.verifyEqual(obj.distributions(i).b, bVals(i) + deltaB);
            end

            % Test 4: 
            % deltaA: array
            % deltaB: scalar
            % inc: false
            obj = GammaDistributionContainer(aVals, bVals);

            deltaA = 0.25 * 1:numOfDistr;
            obj.updateAllDistributionsParams(deltaA, deltaB);

            % Add test for priors - update shouldn't affect priors
            for i = 1:obj.Size
                testCase.verifyTrue(isnan(obj.distributions(i).prior), 'All priors are not NaN');
                testCase.verifyEqual(obj.distributions(i).a, deltaA(i));
                testCase.verifyEqual(obj.distributions(i).b, deltaB);
            end


            % Test 5: 
            % deltaA: scalar
            % deltaB: array
            % inc: true
            obj = GammaDistributionContainer(aVals, bVals);
            
            deltaA = 5;
            deltaB = 0.25 * 1:numOfDistr;
            obj.updateAllDistributionsParams(deltaA, deltaB, true);

            % Add test for priors - update shouldn't affect priors
            for i = 1:obj.Size
                testCase.verifyTrue(isnan(obj.distributions(i).prior), 'All priors are not NaN');
                testCase.verifyEqual(obj.distributions(i).a, aVals(i) + deltaA);
                testCase.verifyEqual(obj.distributions(i).b, bVals(i) + deltaB(i));
            end

            % Test 6: 
            % deltaA: scalar
            % deltaB: array
            % inc: false
            obj = GammaDistributionContainer(aVals, bVals);
            
            deltaA = 5;
            deltaB = 0.25 * 1:numOfDistr;
            obj.updateAllDistributionsParams(deltaA, deltaB);

            % Add test for priors - update shouldn't affect priors
            for i = 1:obj.Size
                testCase.verifyTrue(isnan(obj.distributions(i).prior), 'All priors are not NaN');
                testCase.verifyEqual(obj.distributions(i).a, deltaA);
                testCase.verifyEqual(obj.distributions(i).b, deltaB(i));
            end

            % Test 7: 
            % deltaA: array
            % deltaB: array
            % inc: true
            obj = GammaDistributionContainer(aVals, bVals);
            
            deltaA = 0.1 * 1:numOfDistr;
            deltaB = 0.25 * 1:numOfDistr;
            obj.updateAllDistributionsParams(deltaA, deltaB, true);

            % Add test for priors - update shouldn't affect priors
            for i = 1:obj.Size
                testCase.verifyTrue(isnan(obj.distributions(i).prior), 'All priors are not NaN');
                testCase.verifyEqual(obj.distributions(i).a, aVals(i) + deltaA(i));
                testCase.verifyEqual(obj.distributions(i).b, bVals(i) + deltaB(i));
            end
            
            % Test 8: 
            % deltaA: array
            % deltaB: array
            % inc: false
            obj = GammaDistributionContainer(aVals, bVals);
            
            deltaA = 0.1 * 1:numOfDistr;
            deltaB = 0.25 * 1:numOfDistr;
            obj.updateAllDistributionsParams(deltaA, deltaB);

            % Add test for priors - update shouldn't affect priors
            for i = 1:obj.Size
                testCase.verifyTrue(isnan(obj.distributions(i).prior), 'All priors are not NaN');
                testCase.verifyEqual(obj.distributions(i).a, deltaA(i));
                testCase.verifyEqual(obj.distributions(i).b, deltaB(i));
            end
        end
    end
end
