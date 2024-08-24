classdef GammaContainerTest < matlab.unittest.TestCase
    methods (Test)
        %% Constructor tests
        function testConstructorWithOneParameter(testCase)
            obj = GammaContainer("SD");

            testCase.verifyEqual(obj.Size, 1);

            testCase.verifyEqual(obj.a, Constants.DEFAULT_GAMMA_A);
            testCase.verifyEqual(obj.b, Constants.DEFAULT_GAMMA_B);
            testCase.verifyTrue(obj.prior == Gamma());
        end

        function testConstructorWithTwoParameters(testCase)
            size = 5;
            obj = GammaContainer("SD", size);

            testCase.verifyEqual(obj.Size, size);

            testCase.verifyEqual(obj.a, Constants.DEFAULT_GAMMA_A);
            testCase.verifyEqual(obj.b, repmat(Constants.DEFAULT_GAMMA_B, size, 1));
            testCase.verifyTrue(obj.prior == Gamma());
        end

        function testConstructorWithThreeParameters(testCase)
            size = 5;
            a = 4;
            obj = GammaContainer("SD", size, a);

            testCase.verifyEqual(obj.Size, size);

            testCase.verifyEqual(obj.a, a);
            testCase.verifyEqual(obj.b, repmat(Constants.DEFAULT_GAMMA_B, size, 1));
            testCase.verifyTrue(obj.prior == Gamma());
        end

        function testConstructorWithFourParameters(testCase)
            size = 5;
            a = 4;

            % Test 1: 'b' is scalar
            b = 12;
            obj = GammaContainer("SD", size, a, b);

            testCase.verifyEqual(obj.Size, size);

            testCase.verifyEqual(obj.a, a);
            testCase.verifyEqual(obj.b, repmat(b, size, 1));
            testCase.verifyTrue(obj.prior == Gamma());

            % Test 2: 'b' is column vector
            b = [1; 2; 3; 4; 5];
            obj = GammaContainer("SD", size, a, b);

            testCase.verifyEqual(obj.Size, size);

            testCase.verifyEqual(obj.a, a);
            testCase.verifyEqual(obj.b, b);
            testCase.verifyTrue(obj.prior == Gamma());
        end

        function testConstructorWithFiveParameters(testCase)
            size = 5;
            a = 4;
            aPrior = 3; bPrior = 6;
            prior = Gamma(aPrior, bPrior);

            % Test 1: 'b' is scalar
            b = 12;
            obj = GammaContainer("SD", size, a, b, prior);

            testCase.verifyEqual(obj.Size, size);

            testCase.verifyEqual(obj.a, a);
            testCase.verifyEqual(obj.b, repmat(b, size, 1));
            testCase.verifyTrue(obj.prior == prior);

            % Test 2: 'b' is array
            b = (1:size)';
            obj = GammaContainer("SD", size, a, b, prior);

            testCase.verifyEqual(obj.Size, size);

            testCase.verifyEqual(obj.a, a);
            testCase.verifyEqual(obj.b, b);
            testCase.verifyTrue(obj.prior == prior);
        end





        %% Dependent properties
        function testDependentProperties(testCase)
            % Setup
            size = 4;
            a = 2;
            b = [2; 4; 8; 10];
            obj = GammaContainer("SD", size, a, b);

            % Size, E, E_Diag, H, E_LnP
            testCase.verifyEqual(obj.Size, size);
            testCase.verifyEqual(obj.E, [1; 0.5; 0.25; 0.2]);
            testCase.verifyEqual(obj.E_Diag, diag([1; 0.5; 0.25; 0.2]));
            
            H = 0; E_LnP = 0;
            for i = 1:size
                ds = Gamma(a, b(i), Gamma());
                H = H + ds.H;
                E_LnP = E_LnP + ds.E_LnP;
            end
            testCase.verifyEqual(obj.H, H);
            testCase.verifyTrue(obj.E_LnP - E_LnP < 1e-12);

            % E_LnP (special test)
            a = 14; b = 7;
            aPrior = 1; bPrior = 1;
            prior = Gamma(aPrior, bPrior);
            obj = GammaContainer("SD", size, a, b, prior);
            
            testCase.verifyEqual(obj.E_LnP, -2 * size);

            % E_LnX
            obj.updateAllDistributionsB([exp(1); exp(2); exp(6); exp(10)]);
            testCase.verifyTrue(obj.E_LnX - (4 * psi(a) - 19) < 1e-12);
        end


 


        %% Update methods
        function testUpdateAllDistributionsA(testCase)
            % Setup
            size = 4;
            a = 2;
            b = [2; 4; 8; 10];
            obj = GammaContainer("SD", size, a, b);

            testCase.verifyEqual(obj.E, [1; 0.5; 0.25; 0.2]);
            testCase.verifyTrue(obj.prior == Gamma());

            % Update
            aNew = 16;
            obj.updateAllDistributionsA(aNew);
            testCase.verifyEqual(obj.a, aNew);
            testCase.verifyEqual(obj.E, [8; 4; 2; 1.6]);

            % Update shouldn't affect prior
            testCase.verifyTrue(obj.prior == Gamma());
        end


        function testUpdateAllDistributionsB(testCase)
            % Setup
            size = 4;
            a = 2;
            b = [2; 4; 8; 10];
            obj = GammaContainer("SD", size, a, b);

            testCase.verifyEqual(obj.E, [1; 0.5; 0.25; 0.2]);
            testCase.verifyTrue(obj.prior == Gamma());

            % Test 1: scalar 'b'
            bNew = 16;
            obj.updateAllDistributionsB(bNew);
            testCase.verifyEqual(obj.Size, size);
            testCase.verifyEqual(obj.b, repmat(bNew, size, 1));

            % Update shouldn't affect prior
            testCase.verifyTrue(obj.prior == Gamma());

            % Test 2: array 'b'
            bNew = [1; 2; 4; 5];
            obj.updateAllDistributionsB(bNew);
            testCase.verifyEqual(obj.Size, size);
            testCase.verifyEqual(obj.b, bNew);
            testCase.verifyEqual(obj.E, [2; 1; 0.5; 0.4]);

            % Update shouldn't affect prior
            testCase.verifyTrue(obj.prior == Gamma());
        end


        function testRemoveDistributions(testCase) 
            % Setup
            size = 4;
            a = 2;
            b = [2; 4; 8; 10];
            obj = GammaContainer("SD", size, a, b);

            obj.removeDistributions([1, 2]);

            testCase.verifyEqual(obj.Size, 2);
            testCase.verifyEqual(obj.E, [0.25; 0.2]);

            % Update shouldn't affect prior
            testCase.verifyTrue(obj.prior == Gamma());
        end





        %% Private properties
        function testSetters(testCase)
            a = 1; b = 2; size = 10;
            obj = GammaContainer("SD", size, a, b);

            % Init value is actual expectation
            testCase.verifyTrue(isequal(obj.getExpCInit(), obj.E));
            
            % randomArray = minValue + (maxValue - minValue) * rand(1, n);
            newExpC = (0.1 + (5 - 0.1) * rand(1, obj.Size))';
            obj.setExpCInit(newExpC);

            testCase.verifyTrue(isequal(obj.getExpCInit(), newExpC));
        end
    end
end
