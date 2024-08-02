classdef GammaDistributionTest < matlab.unittest.TestCase
    methods (Static, Access = private)
        function verifyObject(testCase, obj, a, b)
            testCase.verifyEqual(obj.a, a);
            testCase.verifyEqual(obj.b, b);
        end
        
    end

    methods (Test)
        function testConstructor(testCase)
            % Default constructor
            obj = GammaDistribution();

            GammaDistributionTest.verifyObject(testCase, obj, Constants.DEFAULT_GAMMA_A, Constants.DEFAULT_GAMMA_B);

            a = 5;
            % One parameter constructor
            obj = GammaDistribution(a);

            GammaDistributionTest.verifyObject(testCase, obj, a, a);

            b = 7;
            % Two parameters constructor
            obj = GammaDistribution(a, b);

            GammaDistributionTest.verifyObject(testCase, obj, a, b);
        end

        function testDependentProperties(testCase)
            a = 1; b = 2;
            obj = GammaDistribution(a, b);

            testCase.verifyEqual(obj.Expectation, 0.5);
            testCase.verifyEqual(obj.Variance, 0.25);
        end
        
        % TODO (low): This can be split into multiple test methods
        function testUpdateMethods(testCase)
            a = 1; b = 2;
            obj = GammaDistribution(a, b);

            %% updateA
            deltaA = 10;
            for i=1:10
                obj.updateA(deltaA);
                GammaDistributionTest.verifyObject(testCase, obj, a + deltaA * i, b);
            end
    
            obj.updateA(deltaA, false);

            GammaDistributionTest.verifyObject(testCase, obj, deltaA, b);
            
            % Get current parameter values
            a = obj.a; b = obj.b;

            %% updateB
            deltaB = 100;

            for i=1:10
                obj.updateB(deltaB);
                GammaDistributionTest.verifyObject(testCase, obj, a, b + deltaB * i);
            end
    
            obj.updateB(deltaB, false);

            GammaDistributionTest.verifyObject(testCase, obj, a, deltaB);

            a = obj.a; b = obj.b;
            %% updateParameters
            obj.updateParameters(deltaA, deltaB);

            GammaDistributionTest.verifyObject(testCase, obj, a + deltaA, b + deltaB);

            obj.updateParameters(deltaA, deltaB, false);
            
            GammaDistributionTest.verifyObject(testCase, obj, deltaA, deltaB);
        end
    end
end
