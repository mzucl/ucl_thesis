classdef GaussianDistributionContainerTest < matlab.unittest.TestCase
    methods (Test)
        %% Constructors
        function testThreeParameterConstructor(testCase)
            dim = 10;
            numDistributions = 2;

            prior = GaussianDistribution(zeros(dim, 1));
            obj = GaussianDistributionContainer(numDistributions, prior, true);

            testCase.verifyEqual(obj.Size, numDistributions);
            testCase.verifyEqual(obj.cols, true);

            for i = 1:obj.Size
                GaussianDistributionTest.verifyObject(testCase, obj.ds(i), ...
                    prior.mu, prior.cov, prior, prior.dim);
            end
        end



        %% Dependent properties (dependent of the format)
        function testDependentPropertiesColumnFormat(testCase)
            dim = 2;
            cols = true;
            numDistributions = 3;
            prior = GaussianDistribution(zeros(dim, 1), eye(dim));

            obj = GaussianDistributionContainer(numDistributions, prior, cols);

            % Size
            testCase.verifyEqual(obj.Size, numDistributions);

            % Setup
            %   standard normal
            %   mu = 1, cov: diag(2)
            %   mu = [1; 2], cov: newCov
            obj.updateDistributionParams(2, ones(dim, 1), 2);
            newCov = [10, 1; 1, 2];
            obj.updateDistributionParams(3, [1; 2], newCov);

            % E, E_XXt, E_XtX
            for i = 1:numDistributions
                testCase.verifyEqual(obj.E{i}, obj.ds(i).E);
                testCase.verifyEqual(obj.E_XXt{i}, obj.ds(i).E_XXt);
                testCase.verifyEqual(obj.E_XtX{i}, obj.ds(i).E_XtX);
            end

            % EC, E_Ct, E_CtC
            % EC in this format must have dimension (dim x numDistributions)
            testCase.verifyEqual(size(obj.EC, 1), dim);
            testCase.verifyEqual(size(obj.EC, 2), numDistributions);

            expectedVal = [obj.ds(1).E, obj.ds(2).E ...
                    obj.ds(3).E];

            testCase.verifyEqual(obj.EC, expectedVal);
            testCase.verifyEqual(obj.E_Ct, expectedVal');
            testCase.verifyEqual(obj.E_CtC, [[2, 0, 0]; [0, 6, 3]; [0, 3, 17]]);

            % [!!!]
            colSqNorm = obj.getExpectationOfColumnsNormSq();
            testCase.verifyEqual(colSqNorm, [2; 6; 17]);

            % Test 2: Same test when none of the ds are st. normal
            obj.updateDistributionParams(1, [5; 7], [3, 4; 4, 7]);
            testCase.verifyEqual(obj.E_CtC, [[84, 12, 19]; [12, 6, 3]; [19, 3, 17]]);

            colSqNorm = obj.getExpectationOfColumnsNormSq();
            testCase.verifyEqual(colSqNorm, [84; 6; 17]);
        end

        function testDependentPropertiesRowFormat(testCase)
            dim = 2;
            cols = false;
            numDistributions = 3;

            prior = GaussianDistribution(zeros(dim, 1), eye(dim));

            obj = GaussianDistributionContainer(numDistributions, prior, cols);

            % Size
            testCase.verifyEqual(obj.Size, numDistributions);

            % Setup
            %   standard normal
            %   mu = 1, cov: diag(2)
            %   mu = [1; 2], cov: newCov
            obj.updateDistributionParams(2, ones(dim, 1), 2);
            newCov = [10, 1; 1, 2];
            obj.updateDistributionParams(3, [1; 2], newCov);

            % E, E_XXt, E_XtX
            for i = 1:numDistributions
                testCase.verifyEqual(obj.E{i}, obj.ds(i).E);
                testCase.verifyEqual(obj.E_XXt{i}, obj.ds(i).E_XXt);
                testCase.verifyEqual(obj.E_XtX{i}, obj.ds(i).E_XtX);
            end

            testCase.verifyEqual(size(obj.EC, 1), numDistributions);
            testCase.verifyEqual(size(obj.EC, 2), dim);

            % EC, E_Ct, E_CtC
            expectedVal = [obj.ds(1).E'; obj.ds(2).E'; ...
                    obj.ds(3).E'];

            testCase.verifyEqual(obj.EC, expectedVal);
            testCase.verifyEqual(obj.E_Ct, expectedVal');
            testCase.verifyEqual(obj.E_CtC, [[15, 4]; [4, 10]]);

            % [!!!]
            colSqNorm = obj.getExpectationOfColumnsNormSq();
            testCase.verifyEqual(colSqNorm, [15; 10]);

            % Test 2: Same test when none of the ds are st. normal
            obj.updateDistributionParams(1, [5; 7], [3, 4; 4, 7]);
            testCase.verifyEqual(obj.E_CtC, [[42, 43]; [43, 65]]);

            colSqNorm = obj.getExpectationOfColumnsNormSq();
            testCase.verifyEqual(colSqNorm, [42; 65]);
        end



        %% Dependent properties (independent of the format)
        function testEntropyProperties(testCase)
            dim = 2;
            numDistributions = 2;

            % Covariance matrix for each of the components
            compCov = 2 * pi * eye(dim);
            prior = GaussianDistribution(ones(dim, 1), compCov);
            obj = GaussianDistributionContainer(numDistributions, prior, true);

            for i = 1:obj.Size
                testCase.verifyEqual(obj.H{i}, 1 + 2 * log(2*pi));
            end

            testCase.verifyEqual(obj.HC, numDistributions * (1 + 2 * log(2*pi)));
        end
        
        function testPriorPrecisionProperty(testCase)
            % Test 1: spherical covariance
            dim = 10;
            numDistributions = 2;

            prior = GaussianDistribution(zeros(dim, 1), 5 * eye(dim));
            obj = GaussianDistributionContainer(numDistributions, prior, true);
            
            testCase.verifyEqual(obj.PPrec{1}, 1/5);

            % Test 2: diagonal covariance
            dim = 10;
            numDistributions = 2;

            precisions = [2, 2, 2, 2, 2, 4, 4, 4, 4, 4]';
            prior = GaussianDistribution(zeros(dim, 1), diag(precisions));
            obj = GaussianDistributionContainer(numDistributions, prior, true);
            
            for i = 1:numDistributions
                testCase.verifyEqual(obj.PPrec{i}, 1./precisions);
            end

            % Test 3: full covariance
            dim = 10;
            numDistributions = 2;

            prior = GaussianDistribution(zeros(dim, 1), Utility.generateRandomSPDMatrix(dim));
            obj = GaussianDistributionContainer(numDistributions, prior, true);
            
            for i = 1:numDistributions
                testCase.verifyEqual(obj.PPrec{i}, NaN);
            end
        end
        


         %% Private properties
        function testSetters(testCase)
            % Test 1: cols = true
            dim = 10;
            numDistributions = 2;

            prior = GaussianDistribution(zeros(dim, 1), 5 * eye(dim));
            obj = GaussianDistributionContainer(numDistributions, prior, true);

            testCase.verifyTrue(all(isequal(obj.getExpCInit(), obj.EC)));

            % Set init expectation
            newExpC = Utility.generateRandomIntMatrix(dim, numDistributions);
            obj.setExpCInit(newExpC);

            testCase.verifyTrue(isequal(obj.getExpCInit(), newExpC));

            % Test 2: cols = false
            obj = GaussianDistributionContainer(numDistributions, prior, false);

            testCase.verifyTrue(all(isequal(obj.getExpCInit(), obj.EC)));

            % Set init expectation
            newExpC = Utility.generateRandomIntMatrix(numDistributions, dim);
            obj.setExpCInit(newExpC);

            testCase.verifyTrue(isequal(obj.getExpCInit(), newExpC));
        end



        %% Update and retrieve methods
        function testGetDistribution(testCase)
            dim = 10;
            numDistributions = 2;

            prior = GaussianDistribution(zeros(dim, 1));
            obj = GaussianDistributionContainer(numDistributions, prior, true);

            idx = 1;
            % Update 'mu' of the distribution on this index
            newMu = 1:dim;
            obj.updateDistributionMu(idx, newMu);

            dist = obj.getDistribution(idx);

            testCase.verifyEqual(dist.mu, newMu');
            testCase.verifyEqual(dist.cov, prior.cov); % 'cov' hasn't change
        end

        function testUpdateDistribution(testCase)
            dim = 10;
            numDistributions = 2;
            prior = GaussianDistribution(zeros(dim, 1));
            obj = GaussianDistributionContainer(numDistributions, prior, true);

            % New distribution
            covNew = Utility.generateRandomSPDMatrix(dim);
            newDist = GaussianDistribution(ones(dim, 1), covNew); % prior: NaN

            idx = 2;
            obj.updateDistribution(idx, newDist);

            for i = 1:obj.Size
                if i ~= idx
                    GaussianDistributionTest.verifyObject(testCase, obj.ds(i), ...
                        zeros(dim, 1), eye(dim), prior, dim);
                else
                    GaussianDistributionTest.verifyObject(testCase, obj.ds(i), ...
                        newDist.mu, newDist.cov, NaN, newDist.dim); % After update prior should be NaN
                end
            end
            
        end

        function testUpdateDistributionParams(testCase)
            dim = 10;
            numDistributions = 2;

            prior = GaussianDistribution(zeros(dim, 1));
            obj = GaussianDistributionContainer(numDistributions, prior, true);

            covNew = Utility.generateRandomSPDMatrix(dim);
            idx = 1;
            obj.updateDistributionParams(idx, ones(dim, 1), covNew);


            for i = 1:obj.Size
                if i ~= idx
                    GaussianDistributionTest.verifyObject(testCase, obj.ds(i), ...
                        zeros(dim, 1), eye(dim), prior, dim);
                else
                    GaussianDistributionTest.verifyObject(testCase, obj.ds(i), ...
                        ones(dim, 1), covNew, prior, dim);
                end
            end
            
        end
        
        function testUpdateDistributionMu(testCase)
            dim = 10;
            numDistributions = 2;
            prior = GaussianDistribution(zeros(dim,1));
            obj = GaussianDistributionContainer(numDistributions, prior, true);

            idx = 1;
            newMu = 1:dim;
            obj.updateDistributionMu(idx, newMu);

            for i = 1:obj.Size
                if i ~= idx
                    GaussianDistributionTest.verifyObject(testCase, obj.ds(i), ...
                        zeros(dim, 1), eye(dim), prior, dim);
                else
                    GaussianDistributionTest.verifyObject(testCase, obj.ds(i), ...
                        (newMu)', eye(dim), prior, dim);
                end
            end
        end

        function testUpdateAllDistributionsCovariance(testCase)
            dim = 10;
            numDistributions = 2;

            prior = GaussianDistribution(zeros(dim, 1));
            obj = GaussianDistributionContainer(numDistributions, prior, true);

            newCov = Utility.generateRandomSPDMatrix(dim);
            obj.updateAllDistributionsCovariance(newCov);

            % Test that for each distribution value of 'mu' hasn't changed,
            % but the value of 'cov' is set to 'newCov'
            for i = 1:obj.Size
                GaussianDistributionTest.verifyObject(testCase, obj.ds(i), ...
                    obj.ds(i).mu, newCov, prior, dim);
            end
        end
    end
end
