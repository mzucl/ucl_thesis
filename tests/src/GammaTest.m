classdef GammaTest < matlab.unittest.TestCase
    methods (Static)
        function verifyObject(testCase, obj, a, b, prior)
            testCase.verifyEqual(obj.a, a);
            testCase.verifyEqual(obj.b, b);

            % Verify obj.prior if prior is passed in as a parameter
            if nargin > 4
                areEqual = Utility.isNaN(obj.prior) && Utility.isNaN(prior) || ...
                    obj.prior == prior;

                testCase.verifyTrue(areEqual);
            end
        end
    end


    
    methods (Test)
        %% Deep copy and operators overloading
        function testEq(testCase)
            a1 = 1; b1 = 2; 
            a2 = 3; b2 = 4;

            % Test 1: Same objects; no priors;
            obj1 = Gamma(a1, b1);
            obj2 = Gamma(a1, b1);
            testCase.verifyTrue(obj1 == obj2);

            % Test 2: Same objects; same priors;
            obj1 = Gamma(a1, b1, Gamma());
            obj2 = Gamma(a1, b1, Gamma());
            testCase.verifyTrue(obj1 == obj2);

            % Test 3.1: Same objects; different priors (one in NaN);
            obj1 = Gamma(a1, b1, Gamma());
            obj2 = Gamma(a1, b1);
            testCase.verifyTrue(obj1 ~= obj2);

            % Test 3.2: Same objects; different priors (one in NaN);
            obj1 = Gamma(a1, b1);
            obj2 = Gamma(a1, b1, Gamma());
            testCase.verifyTrue(obj1 ~= obj2);

            % Test 3.3: Same objects; different priors;
            obj1 = Gamma(a1, b1, Gamma());
            obj2 = Gamma(a1, b1, Gamma(10, 20));
            testCase.verifyTrue(obj1 ~= obj2);

            % Test 4: Different objects; no priors;
            obj1 = Gamma(a1, b1);
            obj2 = Gamma(a2, b2);
            testCase.verifyTrue(obj1 ~= obj2);

            % Test 5: Different objects; same priors;
            obj1 = Gamma(a1, b1, Gamma());
            obj2 = Gamma(a2, b2, Gamma());
            testCase.verifyTrue(obj1 ~= obj2);

            % Test 6.1: Different objects; different priors (one is NaN);
            obj1 = Gamma(a1, b1, Gamma());
            obj2 = Gamma(a2, b2);
            testCase.verifyTrue(obj1 ~= obj2);

            % Test 6.2: Different objects; different priors;
            obj1 = Gamma(a1, b1, Gamma(10, 20));
            obj2 = Gamma(a2, b2, Gamma(30, 40));
            testCase.verifyTrue(obj1 ~= obj2);
        end

        function testDeepCopy(testCase)
            a = 1; b = 2;
            aNew = 5; bNew = 6;
            obj = Gamma(a, b);

            % Test 1: Deep copy
            deepCopy = obj.copy();
            testCase.verifyTrue(obj == deepCopy);
            %   Update deepCopy
            deepCopy.updateParameters(aNew, bNew);
            testCase.verifyTrue(obj ~= deepCopy);

            % Test 2: Shallow copy
            shallowCopy = obj;
            testCase.verifyTrue(obj == shallowCopy);
            %   Update shallow copy
            shallowCopy.updateParameters(aNew, bNew);
            testCase.verifyTrue(obj == shallowCopy);
        end

        % Deep copy with prior set
        function testDeepCopy2(testCase)
            a = 1; b = 2;
            aNew = 5; bNew = 6;

            obj = Gamma(a, b, Gamma());
            deepCopy = obj.copy();

            testCase.verifyTrue(obj == deepCopy);
          
            obj.updateParameters(aNew, bNew);

            testCase.verifyTrue(obj ~= deepCopy);

            % They still should have the same prior 
            %   (this will test if the prior is copied correctly)
            testCase.verifyTrue(obj.prior == deepCopy.prior);
        end



        %% Constructor
        function testConstructor(testCase)
            % Test 1: Default constructor
            obj = Gamma();
            GammaTest.verifyObject(testCase, obj, Constants.DEFAULT_GAMMA_A, Constants.DEFAULT_GAMMA_B);

            % Test 2: One parameter constructor - parameter is the prior
            truePrior = Gamma();
            obj = Gamma(truePrior);
            GammaTest.verifyObject(testCase, obj, Constants.DEFAULT_GAMMA_A, Constants.DEFAULT_GAMMA_B, truePrior);
            GammaTest.verifyObject(testCase, obj.prior, Constants.DEFAULT_GAMMA_A, Constants.DEFAULT_GAMMA_B, NaN);
          
            % Test 3: One parameter constructor - parameter is numeric,
            % thus value for a
            a = 5;
            obj = Gamma(a);
            GammaTest.verifyObject(testCase, obj, a, a, NaN);

            % Test 4: Two parameters constructor
            b = 7;
            obj = Gamma(a, b);
            GammaTest.verifyObject(testCase, obj, a, b, NaN);

            % Test 5: Three parameters constructor
            aPrior = 1; bPrior = 2;
            prior = Gamma(aPrior, bPrior);
            obj = Gamma(a, b, prior);

            GammaTest.verifyObject(testCase, obj, a, b, prior);
            GammaTest.verifyObject(testCase, obj.prior, aPrior, bPrior, NaN);
        end



        %% Dependent properties
        function testDependentProperties(testCase)
            a = 1; b = 2;
            obj = Gamma(a, b);
            testCase.verifyEqual(obj.E, 0.5);
            testCase.verifyEqual(obj.Var, 0.25);
            testCase.verifyEqual(obj.E_LnP, NaN);

            a = 1; b = 1;
            obj = Gamma(a, b);
            testCase.verifyEqual(obj.E, 1);
            testCase.verifyEqual(obj.Var, 1);
            testCase.verifyEqual(obj.H, 1);
            testCase.verifyEqual(obj.E_LnP, NaN);

            a = 1; b = exp(1);
            obj = Gamma(a, b);
            testCase.verifyEqual(obj.E, exp(-1));
            testCase.verifyEqual(obj.Var, exp(-2));
            testCase.verifyEqual(obj.H, 0);
            testCase.verifyTrue(abs(obj.E_LnX - (psi(obj.a) - log(obj.b))) < 1e-12);
            testCase.verifyEqual(obj.E_LnP, NaN);

            a = 14; b = 7;
            aPrior = 1; bPrior = 1;
            prior = Gamma(aPrior, bPrior);
            obj = Gamma(a, b, prior);

            testCase.verifyEqual(obj.E_LnP, -2);
        end
        


        %% Private properties
        function testSetters(testCase)
            a = 2; b = 1;
            obj = Gamma(a, b);
            testCase.verifyTrue(obj.getExpInit() == obj.E);

            expInit = 34;
            obj.setExpInit(expInit);
            testCase.verifyTrue(obj.getExpInit() == expInit);
        end



        %% Update methods
        function testUpdateMethods(testCase)
            a = 1; b = 2;
            obj = Gamma(a, b);

            % updateA
            newA = 10;
    
            obj.updateA(newA);
            GammaTest.verifyObject(testCase, obj, newA, b);
            
            a = obj.a;

            % updateB
            newB = 100;
            obj.updateB(newB);
            GammaTest.verifyObject(testCase, obj, a, newB);

            % updateParameters
            obj.updateParameters(newA, newB);
            GammaTest.verifyObject(testCase, obj, newA, newB);
        end
    end
end
