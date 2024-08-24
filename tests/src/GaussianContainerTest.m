classdef GaussianContainerTest < matlab.unittest.TestCase
    methods (Test)
        %% Constructors
        function testFourParameterConstructor(testCase) % Default values
            % Test 1: type = "DS"
            type = "DS";
            dim = 5;
            size_ = 2;
            cols = true;

            obj = GaussianContainer(type, size_, cols, dim);

            testCase.verifyEqual(obj.mu, zeros(dim, size_));
            testCase.verifyEqual(obj.cov, eye(dim));

            % Test 2: type = "DD"
            type = "DD";
            dim = 5;
            size_ = 10;
            cols = true;

            obj = GaussianContainer(type, size_, cols, dim);

            testCase.verifyEqual(obj.mu, zeros(dim, size_));
            for i = 1:size_
                testCase.verifyEqual(obj.cov(:, :, i), eye(dim));
            end
        end

        function testFiveParameterConstructor(testCase) % mu
            type = "DS";
            dim = 5;
            size_ = 2;
            cols = true;

            % Test 1: 'mu' is a scalar
            mu = 5;
            obj = GaussianContainer(type, size_, cols, dim, mu);

            testCase.verifyEqual(obj.mu, mu * ones(dim, size_));

            % Test 2: 'mu' is a vector
            mu = (1:dim)';
            obj = GaussianContainer(type, size_, cols, dim, mu);

            for c = 1:size(mu, 2) % for each column in mu
                testCase.verifyEqual(obj.mu(:, c), mu);
            end
            
            % Test 3: 'mu' is a matrix
            mu = Utility.generateRandomIntMatrix(dim, size_);
            obj = GaussianContainer(type, size_, cols, dim, mu);

            testCase.verifyEqual(obj.mu, mu);
        end

        function testSixParameterConstructor(testCase) % mu, cov
            type = "DS";
            dim = 5;
            size_ = 2;
            cols = true;

            % Test 1
            % 'mu': scalar
            % 'cov': scalar
            mu = 5;
            cov = 6;
            obj = GaussianContainer(type, size_, cols, dim, mu, cov);

            testCase.verifyEqual(obj.mu, mu * ones(dim, size_));
            testCase.verifyEqual(obj.cov, cov * eye(dim));

            % Test 2
            % 'mu': scalar
            % 'cov': array
            mu = 5;
            cov = (1:dim)';
            obj = GaussianContainer(type, size_, cols, dim, mu, cov);

            testCase.verifyEqual(obj.mu, mu * ones(dim, size_));
            testCase.verifyEqual(obj.cov, diag(cov));

            % Test 3
            % 'mu': scalar
            % 'cov': matrix
            mu = 5;
            cov = Utility.generateRandomSPDMatrix(dim);
            obj = GaussianContainer(type, size_, cols, dim, mu, cov);

            testCase.verifyEqual(obj.mu, mu * ones(dim, size_));
            testCase.verifyEqual(obj.cov, cov);

            % Test 4
            % 'mu': array
            % 'cov': scalar
            mu = (1:dim)';
            cov = 5;
            obj = GaussianContainer(type, size_, cols, dim, mu, cov);

            for c = 1:size(mu, 2) % for each column in mu
                testCase.verifyEqual(obj.mu(:, c), mu);
            end
            testCase.verifyEqual(obj.cov, cov * eye(dim));

            % Test 5
            % 'mu': array
            % 'cov': array
            mu = (1:dim)';
            cov = (1:dim)';
            obj = GaussianContainer(type, size_, cols, dim, mu, cov);

            for c = 1:size(mu, 2) % for each column in mu
                testCase.verifyEqual(obj.mu(:, c), mu);
            end
            testCase.verifyEqual(obj.cov, diag(cov));

            % Test 6
            % 'mu': array
            % 'cov': matrix
            mu = (1:dim)';
            cov = Utility.generateRandomSPDMatrix(dim);
            obj = GaussianContainer(type, size_, cols, dim, mu, cov);

            for c = 1:size(mu, 2) % for each column in mu
                testCase.verifyEqual(obj.mu(:, c), mu);
            end
            testCase.verifyEqual(obj.cov, cov);

            % Test 7 
            % 'mu': matrix
            % 'cov': scalar
            mu = Utility.generateRandomIntMatrix(dim, size_);
            cov = 5;
            obj = GaussianContainer(type, size_, cols, dim, mu, cov);

            testCase.verifyEqual(obj.mu, mu);
            testCase.verifyEqual(obj.cov, cov * eye(dim));

            % Test 8
            % 'mu': matrix
            % 'cov': array
            mu = Utility.generateRandomIntMatrix(dim, size_);
            cov = (1:dim)';
            obj = GaussianContainer(type, size_, cols, dim, mu, cov);

            testCase.verifyEqual(obj.mu, mu);
            testCase.verifyEqual(obj.cov, diag(cov));

            % Test 9
            % 'mu': matrix
            % 'cov': matrix
            mu = Utility.generateRandomIntMatrix(dim, size_);
            cov = Utility.generateRandomSPDMatrix(dim);
            obj = GaussianContainer(type, size_, cols, dim, mu, cov);

            testCase.verifyEqual(obj.mu, mu);
            testCase.verifyEqual(obj.cov, cov);
        end

        function testSixParameterConstructor_2(testCase) % mu, cov
            type = "DD";
            dim = 5;
            size_ = 2;
            cols = true;

            % Test 1
            % 'mu': scalar
            % 'cov': scalar
            mu = 5;
            cov = 6;
            obj = GaussianContainer(type, size_, cols, dim, mu, cov);
          
            testCase.verifyEqual(size(obj.cov), [dim, dim, size_])
            testCase.verifyEqual(obj.mu, mu * ones(dim, size_));

            for i = 1:size_
                testCase.verifyEqual(obj.cov(:, :, i), cov * eye(dim));
            end

            % Test 2
            % 'mu': scalar
            % 'cov': array
            mu = 5;
            cov = (1:dim)';
            obj = GaussianContainer(type, size_, cols, dim, mu, cov);

            testCase.verifyEqual(obj.mu, mu * ones(dim, size_));
            for i = 1:size_
                testCase.verifyEqual(obj.cov(:, :, i), diag(cov));
            end
            
            % Test 3
            % 'mu': scalar
            % 'cov': matrix
            mu = 5;
            cov = Utility.generateRandomSPDMatrix(dim);
            obj = GaussianContainer(type, size_, cols, dim, mu, cov);

            testCase.verifyEqual(obj.mu, mu * ones(dim, size_));

            for i = 1:size_
                testCase.verifyEqual(obj.cov(:, :, i), cov);
            end

            % Test 4
            % 'mu': array
            % 'cov': scalar
            mu = (1:dim)';
            cov = 5;
            obj = GaussianContainer(type, size_, cols, dim, mu, cov);

            for c = 1:size(mu, 2) % for each column in mu
                testCase.verifyEqual(obj.mu(:, c), mu);
            end

            for i = 1:size_
                testCase.verifyEqual(obj.cov(:, :, i), cov * eye(dim));
            end

            % Test 5
            % 'mu': array
            % 'cov': array
            mu = (1:dim)';
            cov = (1:dim)';
            obj = GaussianContainer(type, size_, cols, dim, mu, cov);

            for c = 1:size(mu, 2) % for each column in mu
                testCase.verifyEqual(obj.mu(:, c), mu);
            end

            for i = 1:size_
                testCase.verifyEqual(obj.cov(:, :, i), diag(cov));
            end
            

            % Test 6
            % 'mu': array
            % 'cov': matrix
            mu = (1:dim)';
            cov = Utility.generateRandomSPDMatrix(dim);
            obj = GaussianContainer(type, size_, cols, dim, mu, cov);

            for c = 1:size(mu, 2) % for each column in mu
                testCase.verifyEqual(obj.mu(:, c), mu);
            end

            for i = 1:size_
                testCase.verifyEqual(obj.cov(:, :, i), cov);
            end

            % Test 7 
            % 'mu': matrix
            % 'cov': scalar
            mu = Utility.generateRandomIntMatrix(dim, size_);
            cov = 5;
            obj = GaussianContainer(type, size_, cols, dim, mu, cov);

            testCase.verifyEqual(obj.mu, mu);
            for i = 1:size_
                testCase.verifyEqual(obj.cov(:, :, i), cov * eye(dim));
            end

            % Test 8
            % 'mu': matrix
            % 'cov': array
            mu = Utility.generateRandomIntMatrix(dim, size_);
            cov = (1:dim)';
            obj = GaussianContainer(type, size_, cols, dim, mu, cov);

            testCase.verifyEqual(obj.mu, mu);
            
            for i = 1:size_
                testCase.verifyEqual(obj.cov(:, :, i), diag(cov));
            end

            % Test 9
            % 'mu': matrix
            % 'cov': matrix
            mu = Utility.generateRandomIntMatrix(dim, size_);
            cov = Utility.generateRandomSPDMatrix(dim);
            obj = GaussianContainer(type, size_, cols, dim, mu, cov);

            testCase.verifyEqual(obj.mu, mu);
            
            for i = 1:size_
                testCase.verifyEqual(obj.cov(:, :, i), cov);
            end
        end

         

        % 
        % 
        % %% Dependent properties (dependent of the format)
        % function testDependentPropertiesColumnFormat(testCase)
        %     dim = 2;
        %     cols = true;
        %     numDistributions = 3;
        %     prior = GaussianDistribution(zeros(dim, 1), eye(dim));
        % 
        %     obj = GaussianContainer(numDistributions, prior, cols);
        % 
        %     % Size
        %     testCase.verifyEqual(obj.Size, numDistributions);
        % 
        %     % Setup
        %     %   standard normal
        %     %   mu = 1, cov: diag(2)
        %     %   mu = [1; 2], cov: newCov
        %     obj.updateDistributionParams(2, ones(dim, 1), 2);
        %     newCov = [10, 1; 1, 2];
        %     obj.updateDistributionParams(3, [1; 2], newCov);
        % 
        %     % E, E_XXt, E_XtX
        %     for i = 1:numDistributions
        %         testCase.verifyEqual(obj.E{i}, obj.ds(i).E);
        %         testCase.verifyEqual(obj.E_XXt{i}, obj.ds(i).E_XXt);
        %         testCase.verifyEqual(obj.E_XtX{i}, obj.ds(i).E_XtX);
        %     end
        % 
        %     % EC, E_Ct, E_CtC
        %     % EC in this format must have dimension (dim x numDistributions)
        %     testCase.verifyEqual(size(obj.EC, 1), dim);
        %     testCase.verifyEqual(size(obj.EC, 2), numDistributions);
        % 
        %     expectedVal = [obj.ds(1).E, obj.ds(2).E ...
        %             obj.ds(3).E];
        % 
        %     testCase.verifyEqual(obj.EC, expectedVal);
        %     testCase.verifyEqual(obj.E_Ct, expectedVal');
        %     testCase.verifyEqual(obj.E_CtC, [[2, 0, 0]; [0, 6, 3]; [0, 3, 17]]);
        % 
        %     testCase.verifyTrue(Utility.isSymmetricMatrix(obj.E_CtC));
        %     testCase.verifyTrue(Utility.isSymmetricMatrix(obj.E_CCt));
        % 
        %     colSqNorm = obj.getExpectationOfColumnsNormSq();
        %     testCase.verifyEqual(colSqNorm, [2; 6; 17]);
        % 
        %     % Test 2: Same test when none of the ds are st. normal
        %     obj.updateDistributionParams(1, [5; 7], [3, 4; 4, 7]);
        %     testCase.verifyEqual(obj.E_CtC, [[84, 12, 19]; [12, 6, 3]; [19, 3, 17]]);
        % 
        %     colSqNorm = obj.getExpectationOfColumnsNormSq();
        %     testCase.verifyEqual(colSqNorm, [84; 6; 17]);
        % 
        %     testCase.verifyTrue(Utility.isSymmetricMatrix(obj.E_CtC));
        %     testCase.verifyTrue(Utility.isSymmetricMatrix(obj.E_CCt));
        % end
        % 
        % function testExpectationOfColumnsNormSq(testCase)
        %     dim = 2;
        %     cols = false;
        %     numDistributions = 3;
        % 
        %     prior = GaussianDistribution(zeros(dim, 1), eye(dim));
        % 
        %     obj = GaussianContainer(numDistributions, prior, cols);
        % 
        %     % Size
        %     testCase.verifyEqual(obj.Size, numDistributions);
        % 
        %     cov = Utility.generateRandomSPDMatrix(2);
        %     % Setup
        %     %   mu = 0, cov
        %     %   mu = 1, cov
        %     %   mu = [1; 2], cov
        %     obj.updateDistributionParams(1, zeros(dim, 1), cov);
        %     obj.updateDistributionParams(2, ones(dim, 1), cov);
        %     obj.updateDistributionParams(3, [1; 2], cov);
        % 
        %     % Test: non-vectorized version
        %     colSqNorm = obj.getExpectationOfColumnsNormSq();
        % 
        %     % Test: vectorized version
        %     MU = [obj.ds(1).mu, obj.ds(2).mu, obj.ds(3).mu];
        %     colSqNorm_Vect = diag(MU * MU') + obj.Size * diag(obj.ds(1).cov);
        % 
        %     testCase.verifyTrue(norm(colSqNorm - colSqNorm_Vect) < 1e-12);
        % end
        % 
        % % cols = true
        % function testDependentPropertiesColumnFormat2(testCase)
        %     dim = 3;
        %     cols = true;
        %     numDistributions = 4;
        %     prior = GaussianDistribution(zeros(dim, 1), eye(dim));
        % 
        %     obj = GaussianContainer(numDistributions, prior, cols);
        % 
        %     % Size
        %     testCase.verifyEqual(obj.Size, numDistributions);
        % 
        %     % Setup
        %     %   mu = [0; 0; 0], cov = eye(3)
        %     %   mu = [1; 2; 3], cov: [1, 0.5, 0; 0.5, 3, 1; 0, 1, 3]
        %     %   mu = [2; 0; 4], cov: [4, -0.5, 0; -0.5, 4, -1; 0, -1, 5]
        %     %   mu = [1; 1; 1], cov: [1, 0, 1.5; 0, 3, 0.5; 1.5, 0.5, 5]
        %     obj.updateDistributionParams(1, zeros(dim, 1), eye(dim));
        %     obj.updateDistributionParams(2, [1; 2; 3], [1, 0.5, 0; 0.5, 3, 1; 0, 1, 3]);
        %     obj.updateDistributionParams(3, [2; 0; 4], [4, -0.5, 0; -0.5, 4, -1; 0, -1, 5]);
        %     obj.updateDistributionParams(4, [1; 1; 1], [1, 0, 1.5; 0, 3, 0.5; 1.5, 0.5, 5]);
        % 
        %     % EC, E_Ct, E_CtC
        %     % EC in this format must have dimension (dim x numDistributions)
        %     testCase.verifyEqual(size(obj.EC, 1), dim);
        %     testCase.verifyEqual(size(obj.EC, 2), numDistributions);
        % 
        %     expectedEC = [0, 1, 2, 1; 0, 2, 0, 1; 0, 3, 4, 1];
        % 
        %     testCase.verifyEqual(obj.EC, expectedEC);
        %     testCase.verifyEqual(obj.E_Ct, expectedEC');
        %     testCase.verifyEqual(obj.E_CtC, [3, 0, 0, 0; 0, 21, 14, 6; 0, 14, 33, 6; 0, 6, 6, 12]);
        % 
        %     colSqNorm = obj.getExpectationOfColumnsNormSq();
        %     testCase.verifyEqual(colSqNorm, [3; 21; 33; 12]);
        % 
        %     testCase.verifyEqual(obj.E_CCt, obj.ds(1).E_XXt + obj.ds(2).E_XXt + ...
        %         obj.ds(3).E_XXt + obj.ds(4).E_XXt)
        % 
        %     testCase.verifyTrue(Utility.isSymmetricMatrix(obj.E_CtC));
        %     testCase.verifyTrue(Utility.isSymmetricMatrix(obj.E_CCt));
        % end
        % 
        % function testDependentPropertiesRowFormat(testCase)
        %     dim = 2;
        %     cols = false;
        %     numDistributions = 3;
        % 
        %     prior = GaussianDistribution(zeros(dim, 1), eye(dim));
        % 
        %     obj = GaussianContainer(numDistributions, prior, cols);
        % 
        %     % Size
        %     testCase.verifyEqual(obj.Size, numDistributions);
        % 
        %     % Setup
        %     %   standard normal
        %     %   mu = 1, cov: diag(2)
        %     %   mu = [1; 2], cov: newCov
        %     obj.updateDistributionParams(2, ones(dim, 1), 2);
        %     newCov = [10, 1; 1, 2];
        %     obj.updateDistributionParams(3, [1; 2], newCov);
        % 
        %     % E, E_XXt, E_XtX
        %     for i = 1:numDistributions
        %         testCase.verifyEqual(obj.E{i}, obj.ds(i).E);
        %         testCase.verifyEqual(obj.E_XXt{i}, obj.ds(i).E_XXt);
        %         testCase.verifyEqual(obj.E_XtX{i}, obj.ds(i).E_XtX);
        %     end
        % 
        %     testCase.verifyEqual(size(obj.EC, 1), numDistributions);
        %     testCase.verifyEqual(size(obj.EC, 2), dim);
        % 
        %     % EC, E_Ct, E_CtC
        %     expectedVal = [obj.ds(1).E'; obj.ds(2).E'; ...
        %             obj.ds(3).E'];
        % 
        %     testCase.verifyEqual(obj.EC, expectedVal);
        %     testCase.verifyEqual(obj.E_Ct, expectedVal');
        %     testCase.verifyEqual(obj.E_CtC, [[15, 4]; [4, 10]]);
        % 
        %     colSqNorm = obj.getExpectationOfColumnsNormSq();
        %     testCase.verifyEqual(colSqNorm, [15; 10]);
        % 
        %     testCase.verifyTrue(Utility.isSymmetricMatrix(obj.E_CtC));
        %     testCase.verifyTrue(Utility.isSymmetricMatrix(obj.E_CCt));
        % 
        %     % Test 2: Same test when none of the ds are st. normal
        %     obj.updateDistributionParams(1, [5; 7], [3, 4; 4, 7]);
        %     testCase.verifyEqual(obj.E_CtC, [[42, 43]; [43, 65]]);
        % 
        %     colSqNorm = obj.getExpectationOfColumnsNormSq();
        %     testCase.verifyEqual(colSqNorm, [42; 65]);
        %     testCase.verifyTrue(Utility.isSymmetricMatrix(obj.E_CtC));
        %     testCase.verifyTrue(Utility.isSymmetricMatrix(obj.E_CCt));
        % end
        % 
        % % cols = false
        % function testDependentPropertiesRowFormat2(testCase)
        %     dim = 3;
        %     cols = false;
        %     numDistributions = 4;
        %     prior = GaussianDistribution(zeros(dim, 1), eye(dim));
        % 
        %     obj = GaussianContainer(numDistributions, prior, cols);
        % 
        %     % Size
        %     testCase.verifyEqual(obj.Size, numDistributions);
        % 
        %     % Setup
        %     %   mu = [0; 0; 0], cov = eye(3)
        %     %   mu = [1; 2; 3], cov: [1, 0.5, 0; 0.5, 3, 1; 0, 1, 3]
        %     %   mu = [2; 0; 4], cov: [4, -0.5, 0; -0.5, 4, -1; 0, -1, 5]
        %     %   mu = [1; 1; 1], cov: [1, 0, 1.5; 0, 3, 0.5; 1.5, 0.5, 5]
        %     obj.updateDistributionParams(1, zeros(dim, 1), eye(dim));
        %     obj.updateDistributionParams(2, [1; 2; 3], [1, 0.5, 0; 0.5, 3, 1; 0, 1, 3]);
        %     obj.updateDistributionParams(3, [2; 0; 4], [4, -0.5, 0; -0.5, 4, -1; 0, -1, 5]);
        %     obj.updateDistributionParams(4, [1; 1; 1], [1, 0, 1.5; 0, 3, 0.5; 1.5, 0.5, 5]);
        % 
        %     % EC, E_Ct, E_CtC
        %     % EC in this format must have dimension (numDistributions x dim)
        %     testCase.verifyEqual(size(obj.EC, 1), numDistributions);
        %     testCase.verifyEqual(size(obj.EC, 2), dim);
        % 
        %     expectedEC = [0, 0, 0; 1, 2, 3; 2, 0, 4; 1, 1, 1];
        % 
        %     testCase.verifyEqual(obj.EC, expectedEC);
        %     testCase.verifyEqual(obj.E_Ct, expectedEC');
        %     testCase.verifyEqual(obj.E_CtC, obj.ds(1).E_XXt + obj.ds(2).E_XXt + ...
        %         obj.ds(3).E_XXt + obj.ds(4).E_XXt)
        % 
        %     testCase.verifyEqual(obj.E_CCt, [3, 0, 0, 0; 0, 21, 14, 6; 0, 14, 33, 6; 0, 6, 6, 12]);
        % 
        %     colSqNorm = obj.getExpectationOfColumnsNormSq();
        %     testCase.verifyEqual(colSqNorm, [13; 16; 40]);
        % 
        %     testCase.verifyTrue(Utility.isSymmetricMatrix(obj.E_CtC));
        %     testCase.verifyTrue(Utility.isSymmetricMatrix(obj.E_CCt));
        % end
        % 
        % % TODO (high): Implement this properly, for now it works for cols
        % % format only!
        % function testTraceDependentProperties(testCase)
        %     dim = 2;
        %     cols = true;
        %     numDistributions = 3;
        %     prior = GaussianDistribution(zeros(dim, 1), eye(dim));
        % 
        %     obj = GaussianContainer(numDistributions, prior, cols);
        % 
        %     % Setup
        %     %   standard normal
        %     %   mu = 1, cov: diag(2)
        %     %   mu = [1; 2], cov: newCov
        %     obj.updateDistributionParams(2, ones(dim, 1), 2);
        %     newCov = [10, 1; 1, 2];
        %     obj.updateDistributionParams(3, [1; 2], newCov);
        % 
        %     testCase.verifyEqual(obj.Tr_CtC, trace(obj.E_CtC));
        % end
        % 
        % 
        % 
        % %% Dependent properties (independent of the format)
        % function testEntropyProperties(testCase)
        %     dim = 2;
        %     numDistributions = 2;
        % 
        %     % Covariance matrix for each of the components
        %     compCov = 2 * pi * eye(dim);
        %     prior = GaussianDistribution(ones(dim, 1), compCov);
        %     obj = GaussianContainer(numDistributions, prior, true);
        % 
        %     for i = 1:obj.Size
        %         testCase.verifyEqual(obj.H{i}, 1 + 2 * log(2*pi));
        %     end
        % 
        %     testCase.verifyEqual(obj.HC, numDistributions * (1 + 2 * log(2*pi)));
        % end
        % 
        % function testPriorPrecisionProperty(testCase)
        %     % Test 1: spherical covariance
        %     dim = 10;
        %     numDistributions = 2;
        % 
        %     prior = GaussianDistribution(zeros(dim, 1), 5 * eye(dim));
        %     obj = GaussianContainer(numDistributions, prior, true);
        % 
        %     testCase.verifyEqual(obj.PPrec{1}, 1/5);
        % 
        %     % Test 2: diagonal covariance
        %     dim = 10;
        %     numDistributions = 2;
        % 
        %     precisions = [2, 2, 2, 2, 2, 4, 4, 4, 4, 4]';
        %     prior = GaussianDistribution(zeros(dim, 1), diag(precisions));
        %     obj = GaussianContainer(numDistributions, prior, true);
        % 
        %     for i = 1:numDistributions
        %         testCase.verifyEqual(obj.PPrec{i}, 1./precisions);
        %     end
        % 
        %     % Test 3: full covariance
        %     dim = 10;
        %     numDistributions = 2;
        % 
        %     prior = GaussianDistribution(zeros(dim, 1), Utility.generateRandomSPDMatrix(dim));
        %     obj = GaussianContainer(numDistributions, prior, true);
        % 
        %     for i = 1:numDistributions
        %         testCase.verifyEqual(obj.PPrec{i}, NaN);
        %     end
        % end
        % 
        % 
        % 
        %  %% Private properties
        % function testSetters(testCase)
        %     % Test 1: cols = true
        %     dim = 10;
        %     numDistributions = 2;
        % 
        %     prior = GaussianDistribution(zeros(dim, 1), 5 * eye(dim));
        %     obj = GaussianContainer(numDistributions, prior, true);
        % 
        %     testCase.verifyTrue(all(isequal(obj.getExpCInit(), obj.EC)));
        % 
        %     % Set init expectation
        %     newExpC = Utility.generateRandomIntMatrix(dim, numDistributions);
        %     obj.setExpCInit(newExpC);
        % 
        %     testCase.verifyTrue(isequal(obj.getExpCInit(), newExpC));
        % 
        %     % Test 2: cols = false
        %     obj = GaussianContainer(numDistributions, prior, false);
        % 
        %     testCase.verifyTrue(all(isequal(obj.getExpCInit(), obj.EC)));
        % 
        %     % Set init expectation
        %     newExpC = Utility.generateRandomIntMatrix(numDistributions, dim);
        %     obj.setExpCInit(newExpC);
        % 
        %     testCase.verifyTrue(isequal(obj.getExpCInit(), newExpC));
        % end
        % 
        % 
        % 
        % %% Update and retrieve methods
        % function testGetDistribution(testCase)
        %     dim = 10;
        %     numDistributions = 2;
        % 
        %     prior = GaussianDistribution(zeros(dim, 1));
        %     obj = GaussianContainer(numDistributions, prior, true);
        % 
        %     idx = 1;
        %     % Update 'mu' of the distribution on this index
        %     newMu = 1:dim;
        %     obj.updateDistributionMu(idx, newMu);
        % 
        %     dist = obj.getDistribution(idx);
        % 
        %     testCase.verifyEqual(dist.mu, newMu');
        %     testCase.verifyEqual(dist.cov, prior.cov); % 'cov' hasn't change
        % end
        % 
        % function testUpdateDistribution(testCase)
        %     dim = 10;
        %     numDistributions = 2;
        %     prior = GaussianDistribution(zeros(dim, 1));
        %     obj = GaussianContainer(numDistributions, prior, true);
        % 
        %     % New distribution
        %     covNew = Utility.generateRandomSPDMatrix(dim);
        %     newDist = GaussianDistribution(ones(dim, 1), covNew); % prior: NaN
        % 
        %     idx = 2;
        %     obj.updateDistribution(idx, newDist);
        % 
        %     for i = 1:obj.Size
        %         if i ~= idx
        %             GaussianDistributionTest.verifyObject(testCase, obj.ds(i), ...
        %                 zeros(dim, 1), eye(dim), prior, dim);
        %         else
        %             GaussianDistributionTest.verifyObject(testCase, obj.ds(i), ...
        %                 newDist.mu, newDist.cov, NaN, newDist.dim); % After update prior should be NaN
        %         end
        %     end
        % 
        % end
        % 
        % function testUpdateDistributionParams(testCase)
        %     dim = 10;
        %     numDistributions = 2;
        % 
        %     prior = GaussianDistribution(zeros(dim, 1));
        %     obj = GaussianContainer(numDistributions, prior, true);
        % 
        %     covNew = Utility.generateRandomSPDMatrix(dim);
        %     idx = 1;
        %     obj.updateDistributionParams(idx, ones(dim, 1), covNew);
        % 
        % 
        %     for i = 1:obj.Size
        %         if i ~= idx
        %             GaussianDistributionTest.verifyObject(testCase, obj.ds(i), ...
        %                 zeros(dim, 1), eye(dim), prior, dim);
        %         else
        %             GaussianDistributionTest.verifyObject(testCase, obj.ds(i), ...
        %                 ones(dim, 1), covNew, prior, dim);
        %         end
        %     end
        % 
        % end
        % 
        % function testUpdateDistributionMu(testCase)
        %     dim = 10;
        %     numDistributions = 2;
        %     prior = GaussianDistribution(zeros(dim,1));
        %     obj = GaussianContainer(numDistributions, prior, true);
        % 
        %     idx = 1;
        %     newMu = 1:dim;
        %     obj.updateDistributionMu(idx, newMu);
        % 
        %     for i = 1:obj.Size
        %         if i ~= idx
        %             GaussianDistributionTest.verifyObject(testCase, obj.ds(i), ...
        %                 zeros(dim, 1), eye(dim), prior, dim);
        %         else
        %             GaussianDistributionTest.verifyObject(testCase, obj.ds(i), ...
        %                 (newMu)', eye(dim), prior, dim);
        %         end
        %     end
        % end
        % 
        % function testUpdateAllDistributionsMu(testCase)
        %     dim = 3;
        %     numDistributions = 2;
        %     prior = GaussianDistribution(zeros(dim,1));
        %     obj = GaussianContainer(numDistributions, prior, true);
        % 
        %     newMU = [[1; 2; 3;], [4; 5; 6]];
        %     obj.updateAllDistributionsMu(newMU);
        % 
        %     for i = 1:obj.Size
        %         GaussianDistributionTest.verifyObject(testCase, obj.ds(i), ...
        %                 newMU(:, i), eye(dim), prior, dim);
        %     end
        % end
        % 
        % function testUpdateDistributionCovariance(testCase)
        %     dim = 10;
        %     numDistributions = 2;
        %     prior = GaussianDistribution(zeros(dim,1));
        %     obj = GaussianContainer(numDistributions, prior, true);
        % 
        %     idx = 1;
        %     newCov = Utility.generateRandomSPDMatrix(dim);
        %     obj.updateDistributionCovariance(idx, newCov);
        % 
        %     for i = 1:obj.Size
        %         if i ~= idx
        %             GaussianDistributionTest.verifyObject(testCase, obj.ds(i), ...
        %                 zeros(dim, 1), eye(dim), prior, dim);
        %         else
        %             GaussianDistributionTest.verifyObject(testCase, obj.ds(i), ...
        %                 zeros(dim, 1), newCov, prior, dim);
        %         end
        %     end
        % end
        % 
        % function testUpdateAllDistributionsCovariance(testCase)
        %     dim = 10;
        %     numDistributions = 2;
        % 
        %     prior = GaussianDistribution(zeros(dim, 1));
        %     obj = GaussianContainer(numDistributions, prior, true);
        % 
        %     newCov = Utility.generateRandomSPDMatrix(dim);
        %     obj.updateAllDistributionsCovariance(newCov);
        % 
        %     % Test that for each distribution value of 'mu' hasn't changed,
        %     % but the value of 'cov' is set to 'newCov'
        %     for i = 1:obj.Size
        %         GaussianDistributionTest.verifyObject(testCase, obj.ds(i), ...
        %             obj.ds(i).mu, newCov, prior, dim);
        %     end
        % end
        % 
        % function testRemoveDimensions(testCase)
        %     dim = 5;
        %     numDistributions = 2;
        %     prior = GaussianDistribution(zeros(dim,1));
        %     obj = GaussianContainer(numDistributions, prior, true);
        % 
        %     indices = [1, 3, 5];
        % 
        %     obj.removeDimensions(indices);
        % 
        %     testCase.verifyEqual(obj.Size, numDistributions);
        %     for i = 1:numDistributions
        %         testCase.verifyEqual(obj.ds(i).dim, 2);
        %     end
        % end
    end
end
