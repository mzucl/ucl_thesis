classdef GaussianDistributionTest < matlab.unittest.TestCase
    %% Static methods
    methods (Static, Access = public)
        function verifyObject(testCase, obj, mu, cov, prior, dim)
            testCase.verifyEqual(obj.mu, mu);
            testCase.verifyEqual(obj.cov, cov);
            if nargin > 5
                areEqual = Utility.isNaN(obj.prior) && Utility.isNaN(prior) || ...
                    obj.prior == prior;
                testCase.verifyTrue(areEqual);
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
            obj = GaussianDistribution(mu, cov);

            % Test 1: Deep copy
            deepCopy = obj.copy();
          
            muNew = 5; covNew = 6;
            deepCopy.mu = muNew;
            deepCopy.cov = covNew;

            % Only 'deepCopy' is updated
            GaussianDistributionTest.verifyObject(testCase, obj, mu, cov);
            GaussianDistributionTest.verifyObject(testCase, deepCopy, muNew, covNew);

            % Test 2: Shallow copy
            shallowCopy = obj;
            shallowCopy.updateParameters(muNew, covNew);
            
            % Both 'obj' and 'shallowCopy' are updated
            GaussianDistributionTest.verifyObject(testCase, obj, muNew, covNew);
            GaussianDistributionTest.verifyObject(testCase, shallowCopy, muNew, covNew);
        end

        function testDeepCopy2(testCase)
            mu = 1; cov = 2;
            muNew = 5; covNew = 6;

            obj = GaussianDistribution(mu, cov, GaussianDistribution());
            deepCopy = obj.copy();

            testCase.verifyTrue(obj == deepCopy);
          
            obj.updateParameters(muNew, covNew);

            testCase.verifyTrue(obj ~= deepCopy);

            % They still should have the same prior 
            %   (this will test if the prior is copied correctly)
            testCase.verifyTrue(obj.prior == deepCopy.prior);
        end
        
        function testEq(testCase)
            mu1 = 1; cov1 = 2; 
            mu2 = 3; cov2 = 4;

            % Test 1: Same objects; no priors;
            obj1 = GaussianDistribution(mu1, cov1);
            obj2 = GaussianDistribution(mu1, cov1);
            testCase.verifyTrue(obj1 == obj2);

            % Test 2: Same objects; same priors;
            obj1 = GaussianDistribution(mu1, cov1, GaussianDistribution());
            obj2 = GaussianDistribution(mu1, cov1, GaussianDistribution());
            testCase.verifyTrue(obj1 == obj2);

            % Test 3.1: Same objects; different priors (one in NaN);
            obj1 = GaussianDistribution(mu1, cov1, GaussianDistribution());
            obj2 = GaussianDistribution(mu1, cov1);
            testCase.verifyTrue(obj1 ~= obj2);

            % Test 3.2: Same objects; different priors (one in NaN);
            obj1 = GaussianDistribution(mu1, cov1);
            obj2 = GaussianDistribution(mu1, cov1, GaussianDistribution());
            testCase.verifyTrue(obj1 ~= obj2);

            % Test 3.3: Same objects; different priors;
            obj1 = GaussianDistribution(mu1, cov1, GaussianDistribution());
            obj2 = GaussianDistribution(mu1, cov1, GaussianDistribution(10, 20));
            testCase.verifyTrue(obj1 ~= obj2);

            % Test 4: Different objects; no priors;
            obj1 = GaussianDistribution(mu1, cov1);
            obj2 = GaussianDistribution(mu2, cov2);
            testCase.verifyTrue(obj1 ~= obj2);

            % Test 5: Different objects; same priors;
            obj1 = GaussianDistribution(mu1, cov1, GaussianDistribution());
            obj2 = GaussianDistribution(mu2, cov2, GaussianDistribution());
            testCase.verifyTrue(obj1 ~= obj2);

            % Test 6.1: Different objects; different priors (one is NaN);
            obj1 = GaussianDistribution(mu1, cov1, GaussianDistribution());
            obj2 = GaussianDistribution(mu2, cov2);
            testCase.verifyTrue(obj1 ~= obj2);

            % Test 6.2: Different objects; different priors;
            obj1 = GaussianDistribution(mu1, cov1, GaussianDistribution(10, 20));
            obj2 = GaussianDistribution(mu2, cov2, GaussianDistribution(30, 40));
            testCase.verifyTrue(obj1 ~= obj2);
        end



        %% Constructors
        function testDefaultConstructor(testCase)
            obj = GaussianDistribution();

            GaussianDistributionTest.verifyObject(testCase, obj, ...
                Constants.DEFAULT_GAUSS_MU, 1/Constants.DEFAULT_GAUSS_PRECISION, ...
                NaN, Constants.DEFAULT_GAUSS_DIM);
        end

        function testOneParameterConstructor(testCase)
            % Test 1: 'mu' is a column vector
            mu = [1; 2; 3];
            obj = GaussianDistribution(mu);

            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu, eye(length(mu)), NaN, length(mu));

            % Test 2: 'mu' is a row vector
            mu = [1, 2, 3];
            obj = GaussianDistribution(mu);

            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu', eye(length(mu)), NaN, length(mu));

            % Test 3: When parameter is a prior distribution
            prior = GaussianDistribution();
            obj = GaussianDistribution(prior);
            
            GaussianDistributionTest.verifyObject(testCase, obj, ...
                prior.mu, prior.cov, prior, prior.dim);

        end

        function testTwoParameterConstructor(testCase)
            % Test 1
            % mu: array
            % cov: scalar
            mu = [1; 2; 3];
            cov = 4;
            obj = GaussianDistribution(mu, cov);

            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu, cov * eye(length(mu)), NaN, length(mu));


            % Test 2
            % mu: array
            % cov: array
            mu = [1; 2; 3; 4; 5];
            cov = [1, 1, 2, 2, 3];
            obj = GaussianDistribution(mu, cov);
            
            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu, diag(cov), NaN, length(mu));

            % Test 3
            % mu: array
            % cov: matrix
            mu = [1; 2; 3];
            cov = Utility.generateRandomSPDMatrix(length(mu));
            obj = GaussianDistribution(mu, cov);
            
            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu, cov, NaN, length(mu));


            % Test 4
            % mu: scalar
            % cov: scalar
            mu = 5;
            cov = 1;
            obj = GaussianDistribution(mu, cov);
            
            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu, cov, NaN, 1);


            % Test 5
            % mu: scalar
            % cov: array
            mu = 5;
            cov = [1; 2; 3];
            obj = GaussianDistribution(mu, cov);
            
            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu * ones(length(cov), 1), diag(cov), NaN, length(cov));


            % Test 6
            % mu: scalar
            % cov: matrix
            mu = 5;
            cov = Utility.generateRandomSPDMatrix(5);
            obj = GaussianDistribution(mu, cov);
            
            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu * ones(length(cov), 1), cov, NaN, size(cov, 1));
        end

        function testThreeParameterConstructor(testCase)
            % Test 1
            % mu: array
            % cov: scalar
            mu = [1; 2; 3];
            cov = 4;

            % 'true value' for the prior
            prior = GaussianDistribution(0, 1/2 * eye(length(mu)));

            obj = GaussianDistribution(mu, cov, prior);

            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu, cov * eye(length(mu)), prior, length(mu));

            % Test 6
            % mu: scalar
            % cov: matrix
            mu = 5;
            cov = Utility.generateRandomSPDMatrix(5);

            % 'true value' for the dim and prior
            trueDim = size(cov, 1);
            prior = GaussianDistribution(0, 1/2 * eye(trueDim));

            obj = GaussianDistribution(mu, cov, prior);
            
            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu * ones(trueDim, 1), cov, prior, trueDim);

            % NaN for 'prior'
            obj = GaussianDistribution(mu, cov, NaN);
            
            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu * ones(trueDim, 1), cov, NaN, trueDim);
        end

        function testFourParameterConstructor(testCase)
            % Test 1: full test
            mu = 5;
            cov = 1;
            dim = 10;

            priorMu = 1:dim;
            prior = GaussianDistribution(priorMu); % mu = [1, 2, ..., dim]

            obj = GaussianDistribution(mu, cov, prior, dim);

            % Verify object
            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu * ones(dim, 1), cov * eye(dim), prior, dim);

            % Test 2: prior is NaN
            obj = GaussianDistribution(mu, cov, NaN, dim);
            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu * ones(dim, 1), cov * eye(dim), NaN, dim);
        end
        

        
        %% Dependent properties
        function testDependentProperties(testCase)
            % Test 1.1: Expectations and Variance
            dim = 10;
            obj = GaussianDistribution(0, 1, NaN, dim);
            
            testCase.verifyEqual(obj.Expectation, zeros(dim, 1));
            testCase.verifyEqual(obj.Variance, eye ...
                (dim));
            testCase.verifyEqual(obj.ExpectationXt, zeros(1, dim));
            testCase.verifyEqual(obj.ExpectationXtX, dim);
            testCase.verifyEqual(obj.ExpectationXXt, eye(dim));

            % Test 1.2: Expectations and Variance
            dim = 10;
            cov = Utility.generateRandomSPDMatrix(dim);
            mu = ones(dim, 1);
            obj = GaussianDistribution(mu, cov);
            
            testCase.verifyEqual(obj.Expectation, mu);
            testCase.verifyEqual(obj.Variance, cov);
            testCase.verifyEqual(obj.ExpectationXt, mu');
            testCase.verifyEqual(obj.ExpectationXtX, mu' * mu + trace(cov));
            testCase.verifyEqual(obj.ExpectationXXt, mu * mu' + cov);

            % Test 2: H
            dim = 2;
            obj = GaussianDistribution(0, 2*pi, NaN, dim);
            testCase.verifyEqual(obj.H, 1 + 2 * log(2*pi));
        end
        
        function testPriorPrecision(testCase)
            % Test 1: covariance is identity matrix
            prior = GaussianDistribution();
            obj = GaussianDistribution(2, 4, prior);
            testCase.verifyEqual(obj.PriorPrecision, Constants.DEFAULT_GAUSS_PRECISION);

            % Test 2: covariance is full matrix
            %   this test can fail if generateRandomSPDMatrix returns
            %   a diagonal matrix, which is possible!
            prior = GaussianDistribution([1, 2], Utility.generateRandomSPDMatrix(2));
            obj = GaussianDistribution([0, 0], 1, prior);
            testCase.verifyEqual(obj.PriorPrecision, NaN);

            % Test 3: spherical covariance
            prior = GaussianDistribution([1, 2], diag(4 * ones(2, 1)));
            obj = GaussianDistribution([0, 0], 1, prior);
            testCase.verifyEqual(obj.PriorPrecision, 1/4);

            % Test 4: diagonal covariance
            prior = GaussianDistribution([1, 2], diag([2, 4]));
            obj = GaussianDistribution([0, 0], 1, prior);
            testCase.verifyEqual(obj.PriorPrecision, [1/2; 1/4]);
        end

        

        %% Private properties
        function testSetters(testCase)
            obj = GammaDistribution();
            testCase.verifyTrue(obj.getExpInit() == obj.Expectation);

            expInit = 34;
            obj.setExpInit(expInit);
            testCase.verifyTrue(obj.getExpInit() == expInit);
        end


        
        %% Update methods
        function testUpdateParameters(testCase)
            % Test 1: No prior
            % ----------------------------------------------------
            dim = 10;
            obj = GaussianDistribution(0, 1, NaN, dim);

            % New values for 'mu' and 'cov'
            newMu = 5;
            newCov = Utility.generateRandomSPDMatrix(dim); % We want to keep the same dimension
            obj.updateParameters(newMu, newCov);

            GaussianDistributionTest.verifyObject(testCase, obj, ...
                newMu * ones(dim, 1), newCov, NaN, dim);

            % Test 2: With prior
            % ----------------------------------------------------
            dim = 10;
            prior = GaussianDistribution(1:dim, Utility.generateRandomSPDMatrix(dim));
            obj = GaussianDistribution(0, 1, prior, dim);

            % New values for 'mu' and 'cov'
            newMu = 13;
            newCov = Utility.generateRandomSPDMatrix(dim); % We want to keep the same dimension

            % Before update
            testCase.verifyTrue(obj.prior == prior);
            % After update (prior is verified in the verifyObject call)
            obj.updateParameters(newMu, newCov);
            GaussianDistributionTest.verifyObject(testCase, obj, ...
                newMu * ones(dim, 1), newCov, prior, dim);
        end

        function testUpdateParameters2(testCase)
            % Test if prior is affected by the update method
            dim = 10;
            priorPrec = 2;
            % 'true value' for the prior
            prior = GaussianDistribution(0, 1/priorPrec * eye(dim));

            obj = GaussianDistribution(0, 1, prior, dim);

            GaussianDistributionTest.verifyObject(testCase, obj, ...
                zeros(dim, 1), eye(dim), prior, dim);

            % Update
            newMu = 5;
            newCov = Utility.generateRandomSPDMatrix(dim); % We want to keep the same dimension
            obj.updateParameters(newMu, newCov);

            GaussianDistributionTest.verifyObject(testCase, obj, ...
                newMu * ones(dim, 1), newCov, prior, dim);
        end

        function testUpdateCovariance(testCase)
            dim = 10;
            priorPrec = 2;
            % 'true value' for the prior
            prior = GaussianDistribution(0, 1/priorPrec * eye(dim));
            obj = GaussianDistribution(0, 1, prior, dim);

            newCov = Utility.generateRandomSPDMatrix(dim); % We want to keep the same dimension

            obj.updateCovariance(newCov);

            GaussianDistributionTest.verifyObject(testCase, obj, ...
                zeros(dim, 1), newCov, prior, dim);
        end

        function testUpdateMu(testCase)
            dim = 10;
            obj = GaussianDistribution(0, 1, NaN, dim);

            newMu = 5 * ones(dim); % We want to keep the same dimension

            obj.updateMu(newMu);

            GaussianDistributionTest.verifyObject(testCase, obj, ...
                newMu, eye(dim), NaN, dim);
        end
    end
end
