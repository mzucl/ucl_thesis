classdef GaussianDistributionTest < matlab.unittest.TestCase
    methods (Static, Access = public)
        function verifyObject(testCase, obj, mu, cov, prior, dim)
            testCase.verifyEqual(obj.mu, mu);
            testCase.verifyEqual(obj.cov, cov);
            if nargin > 5
                testCase.verifyEqual(obj.prior, prior);
                if nargin == 6
                    testCase.verifyEqual(obj.dim, dim);
                end
            end
        end
        
    end

    methods (Test)
        function testDeepCopy(testCase)
            mu = 1; cov = 1;
            obj = GaussianDistribution(mu, cov);

            % Test 1: Deep copy
            objCopy = obj.copy();
          
            muNew = 5; covNew = 6;
            objCopy.mu = muNew;
            objCopy.cov = covNew;

            % Only 'objCopy' is updated
            GaussianDistributionTest.verifyObject(testCase, obj, mu, cov);
            GaussianDistributionTest.verifyObject(testCase, objCopy, muNew, covNew);

            % Test 2: Shallow copy
            objCopySh = obj;
            objCopySh.updateParameters(muNew, covNew);
            
            % Both 'obj' and 'objCopySh' are updated
            GaussianDistributionTest.verifyObject(testCase, obj, muNew, covNew);
            GaussianDistributionTest.verifyObject(testCase, objCopy, muNew, covNew);
        end

        function testDeepCopy2(testCase)
            a = 1; b = 2;
            obj = GammaDistribution(a, b, GammaDistribution());

            % Test 1: Deep copy
            objCopy = obj.copy();
          
            aNew = 5; bNew = 6;
            objCopy.a = aNew;
            objCopy.b = bNew;

            % Only 'objCopy' is updated
            GammaDistributionTest.verifyObject(testCase, obj, a, b);
            GammaDistributionTest.verifyObject(testCase, objCopy, aNew, bNew);
        end

        %% Constructors
        function testDefaultConstructor(testCase)
            obj = GaussianDistribution();

            GaussianDistributionTest.verifyObject(testCase, obj, ...
                Constants.DEFAULT_GAUSS_MU, 1/Constants.DEFAULT_GAUSS_PRECISION, NaN, Constants.DEFAULT_GAUSS_DIM);
        end

        function testOneParameterConstructor(testCase)
            % Test 1
            mu = [1; 2; 3]; % Column vector
            obj = GaussianDistribution(mu);

            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu, eye(length(mu)), NaN, length(mu));

            % Test 2
            mu = [1, 2, 3]; % Row vector
            obj = GaussianDistribution(mu);

            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu', eye(length(mu)), NaN, length(mu));

            % Test 3: When parameter is a prior distribution
            prior = GaussianDistribution();
            obj = GaussianDistribution(prior);
            
            % Test both 'prior' and 'obj'
            GaussianDistributionTest.verifyObject(testCase, obj.prior, ...
                prior.mu, prior.cov, NaN, prior.dim);

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
            mu = 5;
            cov = 1;
            dim = 10;

            priorMu = 1:dim;
            prior = GaussianDistribution(priorMu); % mu = [1, 2, ..., dim]

            obj = GaussianDistribution(mu, cov, prior, dim);

            % Verify prior
            GaussianDistributionTest.verifyObject(testCase, obj.prior, ...
                priorMu', eye(length(priorMu)), NaN, dim);

            % Verify object
            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu * ones(dim, 1), cov * eye(dim), prior, dim);

            obj = GaussianDistribution(mu, cov, NaN, dim); % Four parameter where priorPrec is NaN
            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu * ones(dim, 1), cov * eye(dim), NaN, dim);
        end
        

        
        %% Dependent properties
        function testDependentProperties(testCase)
            % Test 1
            dim = 10;
            obj = GaussianDistribution(0, 1, NaN, dim);
            
            testCase.verifyEqual(obj.Expectation, zeros(dim, 1));
            testCase.verifyEqual(obj.Variance, eye ...
                (dim));
            testCase.verifyEqual(obj.ExpectationXt, zeros(1, dim));

            testCase.verifyEqual(obj.ExpectationXtX, dim);
            testCase.verifyEqual(obj.ExpectationXXt, eye(dim));

            % Test 2
            dim = 10;
            cov = Utility.generateRandomSPDMatrix(dim);
            mu = ones(dim, 1);
            obj = GaussianDistribution(mu, cov);
            
            testCase.verifyEqual(obj.Expectation, mu);
            testCase.verifyEqual(obj.Variance, cov);
            testCase.verifyEqual(obj.ExpectationXt, mu');

            testCase.verifyEqual(obj.ExpectationXtX, mu' * mu + trace(cov));
            testCase.verifyEqual(obj.ExpectationXXt, mu * mu' + cov);

            dim = 2;
            obj = GaussianDistribution(0, 2*pi, NaN, dim);
            testCase.verifyEqual(obj.H, 1 + 2 * log(2*pi));
        end
        
        function testDependentProperties2(testCase)
            % Test 1
            prior = GaussianDistribution();
            obj = GaussianDistribution(2, 4, prior);
            testCase.verifyEqual(obj.PriorPrecision, Constants.DEFAULT_GAUSS_PRECISION);

            % Test 2: this test can fail if generateRandomSPDMatrix returns
            % a diagonal matrix, which is possible!
            prior = GaussianDistribution([1, 2], Utility.generateRandomSPDMatrix(2));
            obj = GaussianDistribution([0, 0], 1, prior);
            testCase.verifyEqual(obj.PriorPrecision, NaN);

            % Test 3
            prior = GaussianDistribution([1, 2], diag(4 * ones(2, 1)));
            obj = GaussianDistribution([0, 0], 1, prior);
            testCase.verifyEqual(obj.PriorPrecision, 1/4);

            % Test 4
            prior = GaussianDistribution([1, 2], diag([2, 4]));
            obj = GaussianDistribution([0, 0], 1, prior);
            testCase.verifyEqual(obj.PriorPrecision, [1/2; 1/4]);
        end

        %% Update methods
        function testUpdateParameters(testCase)
            dim = 10;
            obj = GaussianDistribution(0, 1, NaN, dim);

            mu = 5;
            cov = Utility.generateRandomSPDMatrix(dim); % We want to keep the same dimension

            obj.updateParameters(mu, cov);

            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu * ones(dim, 1), cov, NaN, dim);
        end

        function testUpdateParameters2(testCase)
            % Test if prior is affected by the update method
            dim = 10;
            priorPrec = 2;

            % 'true value' for the prior
            prior = GaussianDistribution(0, 1/priorPrec * eye(dim));

            obj = GaussianDistribution(0, 1, prior, dim);

            % Verify both prior and the obj
            GaussianDistributionTest.verifyObject(testCase, obj.prior, ...
                zeros(dim, 1), 1/priorPrec * eye(dim), NaN, dim);

            GaussianDistributionTest.verifyObject(testCase, obj, ...
                zeros(dim, 1), eye(dim), prior, dim);

           
            % Update
            mu = 5;
            cov = Utility.generateRandomSPDMatrix(dim); % We want to keep the same dimension

            obj.updateParameters(mu, cov);

            % Verify both prior and the obj
            GaussianDistributionTest.verifyObject(testCase, obj.prior, ...
                zeros(dim, 1), 1/priorPrec * eye(dim), NaN, dim);

            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu * ones(dim, 1), cov, prior, dim);
            
        end

        function testUpdateCovariance(testCase)
            dim = 10;
            priorPrec = 2;
            % 'true value' for the prior
            prior = GaussianDistribution(0, 1/priorPrec * eye(dim));

            obj = GaussianDistribution(0, 1, prior, dim);

            cov = Utility.generateRandomSPDMatrix(dim); % We want to keep the same dimension

            obj.updateCovariance(cov);

            GaussianDistributionTest.verifyObject(testCase, obj, ...
                zeros(dim, 1), cov, prior, dim);
        end

        function testUpdateMu(testCase)
            dim = 10;
            obj = GaussianDistribution(0, 1, NaN, dim);

            mu = 5 * ones(dim); % We want to keep the same dimension

            obj.updateMu(mu);

            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu, eye(dim), NaN, dim); % Don't validate 'cov', it hasn't been updated
        end
    end
end
