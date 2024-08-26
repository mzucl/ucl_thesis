classdef GaussianTest < matlab.unittest.TestCase
    %% Static methods
    methods (Static, Access = public)
        function verifyObject(testCase, obj, mu, cov, priorPrec, dim)
            testCase.verifyEqual(obj.mu, mu);
            testCase.verifyEqual(obj.cov, cov);

            if nargin > 5
                testCase.verifyEqual(obj.priorPrec, priorPrec);
                if nargin == 6
                    testCase.verifyEqual(obj.dim, dim);
                end
            end
        end
        
    end





    methods (Test)
        %% Deep copy and operators overloading
        function testDeepCopy(testCase)
            mu = 1; cov = 1;
            obj = Gaussian(1, mu, cov);

            % Test 1: Deep copy
            deepCopy = obj.copy();
          
            muNew = 5; covNew = 6;
            deepCopy.mu = muNew;
            deepCopy.cov = covNew;

            % Only 'deepCopy' is updated
            GaussianTest.verifyObject(testCase, obj, mu, cov);
            GaussianTest.verifyObject(testCase, deepCopy, muNew, covNew);

            % Test 2: Shallow copy
            shallowCopy = obj;
            shallowCopy.updateMu(muNew);
            shallowCopy.updateCovariance(covNew);
            
            % Both 'obj' and 'shallowCopy' are updated
            GaussianTest.verifyObject(testCase, obj, muNew, covNew);
            GaussianTest.verifyObject(testCase, shallowCopy, muNew, covNew);
        end

        function testDeepCopy2(testCase)
            mu = 1; cov = 2;
            muNew = 5; covNew = 6;

            obj = Gaussian(1, mu, cov, 10); % priorPrec = 10
            deepCopy = obj.copy();

            testCase.verifyTrue(obj == deepCopy);
          
            obj.updateMu(muNew);
            obj.updateCovariance(covNew);

            testCase.verifyTrue(obj ~= deepCopy);

            testCase.verifyTrue(obj.priorPrec == deepCopy.priorPrec);
        end
        
        function testEq(testCase)
            mu1 = 1; cov1 = 2; 
            mu2 = 3; cov2 = 4;

            % Test 1: Same objects; no priors;
            obj1 = Gaussian(1, mu1, cov1);
            obj2 = Gaussian(1, mu1, cov1);
            testCase.verifyTrue(obj1 == obj2);

            % Test 2: Same objects; same priors;
            obj1 = Gaussian(1, mu1, cov1, 12);
            obj2 = Gaussian(1, mu1, cov1, 12);
            testCase.verifyTrue(obj1 == obj2);

            % Test 3.1: Same objects; different priors (one in not set);
            obj1 = Gaussian(1, mu1, cov1, 10);
            obj2 = Gaussian(1, mu1, cov1);
            testCase.verifyTrue(obj1 ~= obj2);

            % Test 3.2: Same objects; different priors (one in not set);
            obj1 = Gaussian(1, mu1, cov1);
            obj2 = Gaussian(1, mu1, cov1, 12);
            testCase.verifyTrue(obj1 ~= obj2);

            % Test 3.3: Same objects; different priors;
            obj1 = Gaussian(1, mu1, cov1, 12);
            obj2 = Gaussian(1, mu1, cov1, 14);
            testCase.verifyTrue(obj1 ~= obj2);

            % Test 4: Different objects; no priors;
            obj1 = Gaussian(1, mu1, cov1);
            obj2 = Gaussian(1, mu2, cov2);
            testCase.verifyTrue(obj1 ~= obj2);

            % Test 5: Different objects; same priors;
            obj1 = Gaussian(1, mu1, cov1, 1);
            obj2 = Gaussian(1, mu2, cov2, 1);
            testCase.verifyTrue(obj1 ~= obj2);

            % Test 6.1: Different objects; different priors (one is not set);
            obj1 = Gaussian(1, mu1, cov1, 12);
            obj2 = Gaussian(1, mu2, cov2);
            testCase.verifyTrue(obj1 ~= obj2);

            % Test 6.2: Different objects; different priors;
            obj1 = Gaussian(1, mu1, cov1, 12);
            obj2 = Gaussian(1, mu2, cov2, 14);
            testCase.verifyTrue(obj1 ~= obj2);
        end





        %% Constructors
        function testDefaultConstructor(testCase)
            obj = Gaussian();

            GaussianTest.verifyObject(testCase, obj, ...
                Constants.DEFAULT_GAUSS_MU, 1/Constants.DEFAULT_GAUSS_PRECISION, ...
                Constants.DEFAULT_GAUSS_PRECISION, Constants.DEFAULT_GAUSS_DIM);
        end

        function testOneParameterConstructor(testCase)
            % OPTION 1: Parameter is an instance of Gaussian
            param = Gaussian();
            obj = Gaussian(param);

            GaussianTest.verifyObject(testCase, obj, param.mu, param.cov, param.priorPrec, param.dim);

            % OPTION 2: Parameter is 'dim' (scalar)
            dim = 5;
            obj = Gaussian(dim);
            GaussianTest.verifyObject(testCase, obj, ...
                repmat(Constants.DEFAULT_GAUSS_MU, dim, 1), eye(5), ...
                Constants.DEFAULT_GAUSS_PRECISION, dim);
        end

        function testTwoParameterConstructor(testCase)
            dim = 5;

            % OPTION 1: 'mu' is a scalar
            mu = 7;
            obj = Gaussian(dim, 7);    
            GaussianTest.verifyObject(testCase, obj, repmat(mu, dim, 1), eye(dim), ...
                Constants.DEFAULT_GAUSS_PRECISION, dim);

            % OPTION 2: 'mu' is a column vector
            mu = [1; 2; 3; 4; 5];
            obj = Gaussian(dim, mu);

            GaussianTest.verifyObject(testCase, obj, mu, eye(dim), ...
                Constants.DEFAULT_GAUSS_PRECISION, dim);
        end

        function testThreeParameterConstructor(testCase)
            % Test 1
            % mu: array
            % cov: scalar
            dim = 3;
            mu = [1; 2; 3];
            cov = 4;
            obj = Gaussian(dim, mu, cov);

            GaussianTest.verifyObject(testCase, obj, mu, cov * eye(dim), ...
                Constants.DEFAULT_GAUSS_PRECISION, dim);

            % Test 2
            % mu: array
            % cov: array
            dim = 5;
            mu = [1; 2; 3; 4; 5];
            cov = [1, 1, 2, 2, 3];
            obj = Gaussian(dim, mu, cov);

            GaussianTest.verifyObject(testCase, obj, mu, diag(cov), ...
                Constants.DEFAULT_GAUSS_PRECISION, dim);

            % Test 3
            % mu: array
            % cov: matrix
            dim = 3;
            mu = [1; 2; 3];
            cov = Utility.generateRandomSPDMatrix(dim);
            obj = Gaussian(dim, mu, cov);

            GaussianTest.verifyObject(testCase, obj, mu, cov, ...
                Constants.DEFAULT_GAUSS_PRECISION, dim);

            % Test 4
            % mu: scalar
            % cov: scalar
            dim = 10;
            mu = 5;
            cov = 1;
            obj = Gaussian(dim, mu, cov);

            GaussianTest.verifyObject(testCase, obj, repmat(mu, dim, 1), cov * eye(dim), ...
                Constants.DEFAULT_GAUSS_PRECISION, dim);

            % Test 5
            % mu: scalar
            % cov: array
            dim = 3;
            mu = 5;
            cov = [1; 2; 3];
            obj = Gaussian(dim, mu, cov);

            GaussianTest.verifyObject(testCase, obj, repmat(mu, dim, 1), diag(cov), ...
                Constants.DEFAULT_GAUSS_PRECISION, dim);

            % Test 6
            % mu: scalar
            % cov: matrix
            dim = 5;
            mu = 5;
            cov = Utility.generateRandomSPDMatrix(dim);
            obj = Gaussian(dim, mu, cov);

            GaussianTest.verifyObject(testCase, obj, repmat(mu, dim, 1), cov, ...
                Constants.DEFAULT_GAUSS_PRECISION, dim);
        end

        function testFourParameterConstructor(testCase)
            % Test 1
            % mu: array
            % cov: scalar
            dim = 3;
            mu = [1; 2; 3];
            cov = 4;
            priorPrec = 12;
            obj = Gaussian(dim, mu, cov, priorPrec);

            GaussianTest.verifyObject(testCase, obj, mu, cov * eye(dim), priorPrec, dim);

            % Test 2
            % mu: array
            % cov: array
            dim = 5;
            mu = [1; 2; 3; 4; 5];
            cov = [1, 1, 2, 2, 3];
            obj = Gaussian(dim, mu, cov, priorPrec);

            GaussianTest.verifyObject(testCase, obj, mu, diag(cov), priorPrec, dim);

            % Test 3
            % mu: array
            % cov: matrix
            dim = 3;
            mu = [1; 2; 3];
            cov = Utility.generateRandomSPDMatrix(dim);
            obj = Gaussian(dim, mu, cov, priorPrec);

            GaussianTest.verifyObject(testCase, obj, mu, cov, priorPrec, dim);

            % Test 4
            % mu: scalar
            % cov: scalar
            dim = 10;
            mu = 5;
            cov = 1;
            obj = Gaussian(dim, mu, cov, priorPrec);

            GaussianTest.verifyObject(testCase, obj, repmat(mu, dim, 1), cov * eye(dim), priorPrec, dim);

            % Test 5
            % mu: scalar
            % cov: array
            dim = 3;
            mu = 5;
            cov = [1; 2; 3];
            obj = Gaussian(dim, mu, cov, priorPrec);

            GaussianTest.verifyObject(testCase, obj, repmat(mu, dim, 1), diag(cov), priorPrec, dim);

            % Test 6
            % mu: scalar
            % cov: matrix
            dim = 5;
            mu = 5;
            cov = Utility.generateRandomSPDMatrix(dim);
            obj = Gaussian(dim, mu, cov, priorPrec);

            GaussianTest.verifyObject(testCase, obj, repmat(mu, dim, 1), cov, priorPrec, dim);
        end

        

        
        %% Dependent properties
        function testDependentProperties(testCase)
            % Test 1.1: Expectations and Variance
            dim = 10;
            obj = Gaussian(dim, 0, 1);
            
            testCase.verifyEqual(obj.E, zeros(dim, 1));
            testCase.verifyEqual(obj.Var, eye(dim));
            testCase.verifyEqual(obj.E_Xt, zeros(1, dim));
            testCase.verifyEqual(obj.E_XtX, dim);
            testCase.verifyEqual(obj.E_XXt, eye(dim));

            % Test 1.2: Expectations and Variance
            dim = 10;
            cov = Utility.generateRandomSPDMatrix(dim);
            mu = ones(dim, 1);
            obj = Gaussian(dim, mu, cov);
            
            testCase.verifyEqual(obj.E, mu);
            testCase.verifyEqual(obj.Var, cov);
            testCase.verifyEqual(obj.E_Xt, mu');
            testCase.verifyEqual(obj.E_XtX, mu' * mu + trace(cov));
            testCase.verifyEqual(obj.E_XXt, mu * mu' + cov);

            % Test 2: H
            dim = 2;
            obj = Gaussian(dim, 0, 2*pi);
            testCase.verifyEqual(obj.H, 1 + 2 * log(2*pi));

            % Test 3: E_LnP
            obj = Gaussian(dim, 0, 1, 2 * pi);
            testCase.verifyEqual(obj.E_LnP, -pi * dim);
        end
        




        %% Private properties
        function testSetters(testCase)
            obj = Gaussian();
            testCase.verifyTrue(obj.getExpInit() == obj.E);

            expInit = 34;
            obj.setExpInit(expInit);
            testCase.verifyTrue(obj.getExpInit() == expInit);
        end




        
        %% Update methods
        function testUpdateMu(testCase)
            dim = 10;
            obj = Gaussian(dim, 0, 1);

            newMu = 5 * ones(dim, 1); % We want to keep the same dimension

            obj.updateMu(newMu);

            GaussianTest.verifyObject(testCase, obj, newMu, eye(dim), Constants.DEFAULT_GAUSS_PRECISION, dim);
        end

        function testUpdateCovariance(testCase)
            dim = 10;
            obj = Gaussian(dim, 0);

            testCase.verifyEqual(obj.cov, eye(dim));

            newCov = Utility.generateRandomSPDMatrix(dim);
            newMu = 5 * ones(dim, 1);

            obj.updateParameters(newMu, newCov);

            GaussianTest.verifyObject(testCase, obj, newMu, newCov, Constants.DEFAULT_GAUSS_PRECISION, dim);
        end

        function testUpdateParameters(testCase)
            dim = 10;
            obj = Gaussian(dim, 0);

            testCase.verifyEqual(obj.cov, eye(dim));

            newCov = Utility.generateRandomSPDMatrix(dim);

            obj.updateCovariance(newCov);

            GaussianTest.verifyObject(testCase, obj, zeros(dim, 1), newCov, Constants.DEFAULT_GAUSS_PRECISION, dim);
        end

        function testRemoveDimensions(testCase)
            % Setup
            dim = 4;
            mu = (1:dim)';
            cov = [
                4, 2, 1, 3;
                2, 5, 2, 1;
                1, 2, 3, 0;
                3, 1, 0, 6
            ];

            obj = Gaussian(dim, mu, cov);
            exp = obj.E;

            obj.removeDimensions([1, 3]);

            testCase.verifyTrue(obj.dim == 2);
            testCase.verifyEqual(obj.E, [exp(2); exp(4)]);
            testCase.verifyEqual(obj.cov, [5, 1; 1, 6]);
            testCase.verifyEqual(obj.mu, [2; 4]);
            testCase.verifyEqual(obj.priorPrec, Constants.DEFAULT_GAUSS_PRECISION);
        end
    end
end
