classdef GammaDistributionTest < matlab.unittest.TestCase
    methods (Static)
        function verifyObject(testCase, obj, a, b, prior)
            testCase.verifyEqual(obj.a, a);
            testCase.verifyEqual(obj.b, b);
            % Verify obj.prior if prior is passed in as a parameter
            if nargin > 4
                testCase.verifyEqual(obj.prior, prior);
            end
        end
        
    end

    methods (Test)
        function testConstructor(testCase)
            % Default constructor
            obj = GammaDistribution();

            GammaDistributionTest.verifyObject(testCase, obj, Constants.DEFAULT_GAMMA_A, Constants.DEFAULT_GAMMA_B);

            % One parameter constructor - the only parameter is the prior
            truePrior = GammaDistribution();
            obj = GammaDistribution(truePrior);

            % Validate prior and obj
            GammaDistributionTest.verifyObject(testCase, obj.prior, Constants.DEFAULT_GAMMA_A, Constants.DEFAULT_GAMMA_B, NaN);
            GammaDistributionTest.verifyObject(testCase, obj, Constants.DEFAULT_GAMMA_A, Constants.DEFAULT_GAMMA_B, truePrior);

            a = 5;
            % One parameter constructor
            obj = GammaDistribution(a);

            GammaDistributionTest.verifyObject(testCase, obj, a, a);

            b = 7;
            % Two parameters constructor
            obj = GammaDistribution(a, b);

            GammaDistributionTest.verifyObject(testCase, obj, a, b);

            aPrior = 1; bPrior = 2;
            % Three parameters constructor
            prior = GammaDistribution(aPrior, bPrior);
            obj = GammaDistribution(a, b, prior);

            GammaDistributionTest.verifyObject(testCase, obj, a, b, prior);
            GammaDistributionTest.verifyObject(testCase, obj.prior, aPrior, bPrior);

        end

        function testDependentProperties(testCase)
            a = 1; b = 2;
            obj = GammaDistribution(a, b);

            testCase.verifyEqual(obj.Expectation, 0.5);
            testCase.verifyEqual(obj.Variance, 0.25);

            a = 1; b = 1;
            obj.updateParameters(a, b, false);
            testCase.verifyEqual(obj.H, 1);

            a = 1; b = exp(1);
            obj.updateParameters(a, b, false);
            testCase.verifyEqual(obj.H, 0);

            testCase.verifyEqual(obj.ExpectationLn, psi(obj.a) - log(obj.b));

            aPrior = 1; bPrior = 1;
            a = 14; b = 7;
            prior = GammaDistribution(aPrior, bPrior);
            obj = GammaDistribution(a, b, prior);

            testCase.verifyEqual(obj.ExpectationLnP, -2);
        end
        
        % TODO (low): This can be split into multiple test methods
        function testUpdateMethods(testCase)
            a = 1; b = 2;
            obj = GammaDistribution(a, b);

            %% updateA
            deltaA = 10;
            for i=1:10
                obj.updateA(deltaA, true);
                GammaDistributionTest.verifyObject(testCase, obj, a + deltaA * i, b);
            end
    
            obj.updateA(deltaA);

            GammaDistributionTest.verifyObject(testCase, obj, deltaA, b);
            
            % Get current parameter values
            a = obj.a; b = obj.b;

            %% updateB
            deltaB = 100;

            for i=1:10
                obj.updateB(deltaB, true);
                GammaDistributionTest.verifyObject(testCase, obj, a, b + deltaB * i);
            end
    
            obj.updateB(deltaB);

            GammaDistributionTest.verifyObject(testCase, obj, a, deltaB);

            a = obj.a; b = obj.b;
            %% updateParameters
            obj.updateParameters(deltaA, deltaB, true);

            GammaDistributionTest.verifyObject(testCase, obj, a + deltaA, b + deltaB);

            obj.updateParameters(deltaA, deltaB);
            
            GammaDistributionTest.verifyObject(testCase, obj, deltaA, deltaB);
        end
    end
end
