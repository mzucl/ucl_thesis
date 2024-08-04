classdef GaussianDistributionTest < matlab.unittest.TestCase
    methods (Static, Access = public)
        function verifyObject(testCase, obj, mu, cov, dim)
            if ~isnan(mu)
                testCase.verifyEqual(obj.mu, mu);
            end
            if ~isnan(cov)
                testCase.verifyEqual(obj.cov, cov);
            end
            testCase.verifyEqual(obj.dim, dim);
        end
        
    end

    methods (Test)
        %% Constructors
        function testDefaultConstructor(testCase)
            obj = GaussianDistribution();

            GaussianDistributionTest.verifyObject(testCase, obj, ...
                Constants.DEFAULT_GAUSS_MU, Constants.DEFAULT_GAUSS_COV, Constants.DEFAULT_GAUSS_DIM);
        end

        function testOneParameterConstructor(testCase)
            mu = [1; 2; 3];
            obj = GaussianDistribution(mu);

            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu, eye(length(mu)), length(mu));
        end

        function testTwoParameterConstructor(testCase)
            % Test 1
            % mu: array
            % cov: scalar
            mu = [1; 2; 3];
            cov = 4;
            obj = GaussianDistribution(mu, cov);

            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu, cov * eye(length(mu)), length(mu));


            % Test 2
            % mu: array
            % cov: array
            mu = [1; 2; 3; 4; 5];
            cov = [1, 1, 2, 2, 3];
            obj = GaussianDistribution(mu, cov);
            
            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu, diag(cov), length(mu));

            % Test 3
            % mu: array
            % cov: matrix
            mu = [1; 2; 3];
            cov = Utility.generateRandomSPDMatrix(length(mu));
            obj = GaussianDistribution(mu, cov);
            
            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu, cov, length(mu));


            % Test 4
            % mu: scalar
            % cov: scalar
            mu = 5;
            cov = 1;
            obj = GaussianDistribution(mu, cov);
            
            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu, cov, 1);


            % Test 5
            % mu: scalar
            % cov: array
            mu = 5;
            cov = [1; 2; 3];
            obj = GaussianDistribution(mu, cov);
            
            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu * ones(length(cov), 1), diag(cov), length(cov));


            % Test 6
            % mu: scalar
            % cov: matrix
            mu = 5;
            cov = Utility.generateRandomSPDMatrix(5);
            obj = GaussianDistribution(mu, cov);
            
            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu * ones(length(cov), 1), cov, size(cov, 1));
        end

        function testThreeParameterConstructor(testCase)
            mu = 5;
            cov = 1;
            dim = 10;

            obj = GaussianDistribution(mu, cov, dim);

            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu * ones(dim, 1), cov * eye(dim), dim);
        end
        

        
        %% Dependent properties
        function testDependentProperties(testCase)
            % Test 1
            dim = 10;
            obj = GaussianDistribution(0, 1, dim);
            
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
            obj = GaussianDistribution(0, 2*pi, dim);
            testCase.verifyEqual(obj.H, 1 + 2 * log(2*pi));
        end
        


        %% Update methods
        function testUpdateParameters(testCase)
            dim = 10;
            obj = GaussianDistribution(0, 1, dim);

            mu = 5;
            cov = Utility.generateRandomSPDMatrix(dim); % We want to keep the same dimension

            obj.updateParameters(mu, cov);

            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu * ones(dim, 1), cov, dim);
        end

        function testUpdateCovariance(testCase)
            dim = 10;
            obj = GaussianDistribution(0, 1, dim);

            cov = Utility.generateRandomSPDMatrix(dim); % We want to keep the same dimension

            obj.updateCovariance(cov);

            GaussianDistributionTest.verifyObject(testCase, obj, ...
                NaN, cov, dim); % Don't validate 'mu', it hasn't been updated
        end

        function testUpdateMu(testCase)
            dim = 10;
            obj = GaussianDistribution(0, 1, dim);

            mu = 5 * ones(dim); % We want to keep the same dimension

            obj.updateMu(mu);

            GaussianDistributionTest.verifyObject(testCase, obj, ...
                mu, NaN, dim); % Don't validate 'cov', it hasn't been updated
        end
    end
end
