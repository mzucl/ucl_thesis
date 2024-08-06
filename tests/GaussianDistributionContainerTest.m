classdef GaussianDistributionContainerTest < matlab.unittest.TestCase
    % methods (Static, Access = private)
    %     function verifyObject(testCase, obj, mu, cov, dim)
    %         if ~isempty(mu)
    %             testCase.verifyEqual(obj.mu, mu);
    %         end
    %         if ~isempty(cov)
    %             testCase.verifyEqual(obj.cov, cov);
    %         end
    %         testCase.verifyEqual(obj.dim, dim);
    %     end
    % 
    % end

    methods (Test)
        %% Constructors
        function testThreeParameterConstructor(testCase)
            dim = 10;
            numDistributions = 2;
            obj = GaussianDistributionContainer(dim, true, numDistributions);

            testCase.verifyEqual(obj.Size, numDistributions);
            testCase.verifyEqual(obj.cols, true);

            % Test that each distribution is a multivariate standard normal
            for i = 1:obj.Size
                GaussianDistributionTest.verifyObject(testCase, obj.distributions(i), zeros(dim, 1), eye(dim), dim);
            end
        end

        function testFourParameterConstructor(testCase)
            dim = 10;
            numDistributions = 2;
            prec = 1:numDistributions;

            obj = GaussianDistributionContainer(dim, false, numDistributions, prec);

            testCase.verifyEqual(obj.Size, numDistributions);
            testCase.verifyEqual(obj.cols, false);

            % Test that each distribution is a multivariate standard normal
            for i = 1:obj.Size
                GaussianDistributionTest.verifyObject(testCase, obj.distributions(i), zeros(dim, 1), 1./prec(i) * eye(dim), dim);
            end
        end

        function testFourParameterConstructorWithCovariance(testCase)
            dim = 5;
            numDistributions = 2;
            prec = diag(1:dim); % 'prec' is a covariance matrix now

            obj = GaussianDistributionContainer(dim, false, numDistributions, prec);

            testCase.verifyEqual(obj.Size, numDistributions);
            testCase.verifyEqual(obj.cols, false);

            % Test that each distribution is a multivariate normal with
            % diagonal covariance matrix
            for i = 1:obj.Size
                GaussianDistributionTest.verifyObject(testCase, obj.distributions(i), zeros(dim, 1), prec, dim);
            end
        end



        %% Dependent properties
        function testDependentPropertiesColumnFormat(testCase)
            dim = 2;
            cols = true;
            numDistributions = 3;

            obj = GaussianDistributionContainer(dim, cols, numDistributions);

            % Size
            testCase.verifyEqual(obj.Size, numDistributions);

            % Setup
            %   standard normal
            %   mu = 1, cov: diag(2)
            %   mu = [1; 2], cov: newCov

            obj.updateDistributionParams(2, ones(dim, 1), 2);
            newCov = [10, 1; 1, 2];
            obj.updateDistributionParams(3, [1; 2], newCov);

            % assignin('base', 'objT', obj);

            % Expectation, ExpectationXXt, ExpectationXtX
            for i = 1:numDistributions
                testCase.verifyEqual(obj.Expectation{i}, obj.distributions(i).Expectation);
                testCase.verifyEqual(obj.ExpectationXXt{i}, obj.distributions(i).ExpectationXXt);
                testCase.verifyEqual(obj.ExpectationXtX{i}, obj.distributions(i).ExpectationXtX);
            end

            % ExpectationC, ExpectationCt, ExpectationCtC

            expectedVal = [obj.distributions(1).Expectation, obj.distributions(2).Expectation ...
                    obj.distributions(3).Expectation];

            testCase.verifyEqual(obj.ExpectationC, expectedVal);
            testCase.verifyEqual(obj.ExpectationCt, expectedVal');
            testCase.verifyEqual(obj.ExpectationCtC, [[2, 0, 0]; [0, 6, 3]; [0, 3, 17]]);
        end

        function testDependentPropertiesRowFormat(testCase)
            dim = 2;
            cols = false;
            numDistributions = 3;

            obj = GaussianDistributionContainer(dim, cols, numDistributions);

            % Size
            testCase.verifyEqual(obj.Size, numDistributions);

            % Setup
            %   standard normal
            %   mu = 1, cov: diag(2)
            %   mu = [1; 2], cov: newCov

            obj.updateDistributionParams(2, ones(dim, 1), 2);
            newCov = [10, 1; 1, 2];
            obj.updateDistributionParams(3, [1; 2], newCov);

            % assignin('base', 'objT', obj);

            % Expectation, ExpectationXXt, ExpectationXtX
            for i = 1:numDistributions
                testCase.verifyEqual(obj.Expectation{i}, obj.distributions(i).Expectation);
                testCase.verifyEqual(obj.ExpectationXXt{i}, obj.distributions(i).ExpectationXXt);
                testCase.verifyEqual(obj.ExpectationXtX{i}, obj.distributions(i).ExpectationXtX);
            end

            % ExpectationC, ExpectationCt, ExpectationCtC

            expectedVal = [obj.distributions(1).Expectation'; obj.distributions(2).Expectation'; ...
                    obj.distributions(3).Expectation'];
            % else
            %     expectedVal = [obj.distributions(1).Expectation'; obj.distributions(2).Expectation'; ...
            %         obj.distributions(3).Expectation'];
                
            % end
            testCase.verifyEqual(obj.ExpectationC, expectedVal);
            testCase.verifyEqual(obj.ExpectationCt, expectedVal');
            testCase.verifyEqual(obj.ExpectationCtC, [[15, 4]; [4, 10]]);

            col1SqNorm = obj.getExpectationOfColumnNormSq();
            testCase.verifyEqual(col1SqNorm, [107 14]);
        end


        
        % [NOTE] Entropy properties are independent of the format col/row
        function testEntropyProperties(testCase)
            numDistributions = 3;
            dim = 2;
            
            % Covariance matrix for each of the components
            compCov = 2 * pi * eye(dim);
            obj = GaussianDistributionContainer(dim, true, numDistributions, compCov);
            
            for i = 1:obj.Size
                testCase.verifyEqual(obj.H{i}, 1 + 2 * log(2*pi));
            end

            testCase.verifyEqual(obj.HC, numDistributions * (1 + 2 * log(2*pi)));
        end
        
        % [NOTE] PriorPrecision property is independent of the format col/row
        function testPriorPrecisionProperty(testCase)
            testCase.verifyEqual(3, 4); % Fail on purpose
            % numDistributions = 3;
            % dim = 2;
            % 
            % % Covariance matrix for each of the components
            % compCov = 2 * pi * eye(dim);
            % obj = GaussianDistributionContainer(dim, true, numDistributions, compCov);
            % 
            % for i = 1:obj.Size
            %     testCase.verifyEqual(obj.H{i}, 1 + 2 * log(2*pi));
            % end
            % 
            % testCase.verifyEqual(obj.HC, numDistributions * (1 + 2 * log(2*pi)));
        end
        

        %% Update/retrieve methods
        function testGetDistribution(testCase)
            dim = 10;
            numDistributions = 5;
            prec = 1:numDistributions;

            obj = GaussianDistributionContainer(dim, false, numDistributions, prec);

            idx = 4;

            GaussianDistributionTest.verifyObject(testCase, obj.distributions(idx), zeros(dim, 1), 1/prec(idx) * eye(dim), dim);
        end

        function testUpdateDistribution(testCase)
            dim = 10;
            numDistributions = 5;

            obj = GaussianDistributionContainer(dim, false, numDistributions);

            covNew = Utility.generateRandomSPDMatrix(dim);
            newDist = GaussianDistribution(ones(dim, 1), covNew);

            idx = 3;
            obj.updateDistribution(idx, newDist);


            % Test that each distribution, except the updated one is a multivariate standard normal
            for i = 1:obj.Size
                if i ~= idx
                    GaussianDistributionTest.verifyObject(testCase, obj.distributions(i), zeros(dim, 1), eye(dim), dim);
                else
                    GaussianDistributionTest.verifyObject(testCase, obj.distributions(i), newDist.mu, newDist.cov, newDist.dim);
                end
            end
            
        end

        function testUpdateDistributionParams(testCase)
            dim = 10;
            numDistributions = 5;

            obj = GaussianDistributionContainer(dim, false, numDistributions);

            covNew = Utility.generateRandomSPDMatrix(dim);
            idx = 3;
            obj.updateDistributionParams(idx, ones(dim, 1), covNew);


            % Test that each distribution, except the updated one is a multivariate standard normal
            for i = 1:obj.Size
                if i ~= idx
                    GaussianDistributionTest.verifyObject(testCase, obj.distributions(i), zeros(dim, 1), eye(dim), dim);
                else
                    GaussianDistributionTest.verifyObject(testCase, obj.distributions(i), ones(dim, 1), covNew, dim);
                end
            end
            
        end
        
        function testUpdateDistributionMu(testCase)
            dim = 10;
            numDistributions = 5;

            obj = GaussianDistributionContainer(dim, false, numDistributions);

            idx = 3;
            obj.updateDistributionMu(idx, 1:dim);


            % Test that each distribution, except the updated one is a multivariate standard normal
            for i = 1:obj.Size
                if i ~= idx
                    GaussianDistributionTest.verifyObject(testCase, obj.distributions(i), zeros(dim, 1), eye(dim), dim);
                else
                    GaussianDistributionTest.verifyObject(testCase, obj.distributions(i), 1:dim, eye(dim), dim);
                end
            end
            
        end

        function testUpdateAllDistributionsCovariance(testCase)
            dim = 10;
            numDistributions = 5;

            obj = GaussianDistributionContainer(dim, false, numDistributions);

            newCov = Utility.generateRandomSPDMatrix(dim);
            obj.updateAllDistributionsCovariance(newCov);

            % Test that for each distribution value of 'mu' hasn't changed,
            % but the value of 'cov' is set to 'newCov'
            for i = 1:obj.Size
                GaussianDistributionTest.verifyObject(testCase, obj.distributions(i), obj.distributions(i).mu, newCov, dim);
            end
            
        end
    end
end
