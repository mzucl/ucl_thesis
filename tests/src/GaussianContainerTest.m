classdef GaussianContainerTest < matlab.unittest.TestCase
    methods (Test)
        %% Constructors
        % 'mu' and 'cov' don't depend on 'cols'
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
            cols = false;

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
            testCase.verifyEqual(obj.priorPrec, Constants.DEFAULT_GAUSS_PRECISION);
            
            for i = 1:size_
                testCase.verifyEqual(obj.cov(:, :, i), cov);
            end
        end
   
        function testSevenParameterConstructor(testCase)
            % Test 1
            type = "DS";
            dim = 5;
            size_ = 2;
            cols = false;
            priorPrec = 100;
           
            mu = 5;
            cov = 6;
            obj = GaussianContainer(type, size_, cols, dim, mu, cov, priorPrec);

            testCase.verifyEqual(obj.mu, mu * ones(dim, size_));
            testCase.verifyEqual(obj.cov, cov * eye(dim));
            testCase.verifyEqual(obj.priorPrec, priorPrec);

            % --------------------------------------------------------
            % Test 2
            type = "DD";
            dim = 5;
            size_ = 2;
            cols = true;
            priorPrec = [10, 100];

            mu = 5;
            cov = 6;
            obj = GaussianContainer(type, size_, cols, dim, mu, cov, priorPrec);
          
            testCase.verifyEqual(size(obj.cov), [dim, dim, size_])
            testCase.verifyEqual(obj.mu, mu * ones(dim, size_));
            testCase.verifyEqual(obj.priorPrec, priorPrec);

            for i = 1:size_
                testCase.verifyEqual(obj.cov(:, :, i), cov * eye(dim));
            end

        end
        

        


        %% Update methods
        function testUpdateDistributionsMu(testCase)
            % Test 1: cols = true
            type = "DS";
            dim = 5;
            size_ = 2;
            cols = true;

            obj = GaussianContainer(type, size_, cols, dim);

            testCase.verifyEqual(obj.mu, zeros(dim, size_));
            testCase.verifyEqual(obj.E, zeros(dim, size_));
            testCase.verifyEqual(obj.cov, eye(dim));

            newMu = Utility.generateRandomIntMatrix(dim, size_);
            obj.updateDistributionsMu(newMu);

            testCase.verifyEqual(obj.mu, newMu);
            testCase.verifyEqual(obj.cov, eye(dim)); % No change
            testCase.verifyEqual(obj.E, newMu);

            
            % Test 1: cols = false
            type = "DS";
            dim = 5;
            size_ = 2;
            cols = false;

            obj = GaussianContainer(type, size_, cols, dim);

            testCase.verifyEqual(obj.mu, zeros(dim, size_));
            testCase.verifyEqual(obj.E, zeros(size_, dim));
            testCase.verifyEqual(obj.cov, eye(dim));

            newMu = Utility.generateRandomIntMatrix(dim, size_);
            obj.updateDistributionsMu(newMu);

            testCase.verifyEqual(obj.mu, newMu);
            testCase.verifyEqual(obj.cov, eye(dim)); % No change
            testCase.verifyEqual(obj.E_Xt, newMu);
        end

        function testUpdateDistributionsCovariance(testCase)
            % Test 1: type = "DS"
            type = "DS";
            dim = 5;
            size_ = 2;
            cols = true;

            obj = GaussianContainer(type, size_, cols, dim);

            testCase.verifyEqual(obj.mu, zeros(dim, size_));
            testCase.verifyEqual(obj.cov, eye(dim));

            newCov = Utility.generateRandomSPDMatrix(dim);
            % Update
            obj.updateDistributionsCovariance(newCov);

            testCase.verifyEqual(obj.mu, zeros(dim, size_)); % No change
            testCase.verifyEqual(obj.cov, newCov); 

            % Test 2: type = "DD"
            type = "DD";
            dim = 5;
            size_ = 2;
            cols = true;

            obj = GaussianContainer(type, size_, cols, dim);

            testCase.verifyEqual(obj.mu, zeros(dim, size_));
            for i = 1:size_
                testCase.verifyEqual(obj.cov(:, :, i), eye(dim));
            end

            newCov = zeros(dim, dim, size_); % Preallocate
            newCov1 = Utility.generateRandomSPDMatrix(dim);
            newCov2 = Utility.generateRandomSPDMatrix(dim);
            newCov(:, :, 1) = newCov1;
            newCov(:, :, 2) = newCov2;

            % Update
            obj.updateDistributionsCovariance(newCov);

            testCase.verifyEqual(obj.mu, zeros(dim, size_)); % No change
            testCase.verifyEqual(obj.cov, newCov);
        end

        function testUpdateDistributionsParameters(testCase)
            % Test 1: type = "DS"
            type = "DS";
            dim = 5;
            size_ = 2;
            cols = true;

            obj = GaussianContainer(type, size_, cols, dim);

            testCase.verifyEqual(obj.mu, zeros(dim, size_));
            testCase.verifyEqual(obj.cov, eye(dim));

            newMu = Utility.generateRandomIntMatrix(dim, size_);
            newCov = Utility.generateRandomSPDMatrix(dim);

            % Update
            obj.updateDistributionsParameters(newMu, newCov);

            testCase.verifyEqual(obj.mu, newMu);
            testCase.verifyEqual(obj.cov, newCov); 

            % Test 2: type = "DD"
            type = "DD";
            dim = 5;
            size_ = 2;
            cols = true;

            obj = GaussianContainer(type, size_, cols, dim);

            testCase.verifyEqual(obj.mu, zeros(dim, size_));
            for i = 1:size_
                testCase.verifyEqual(obj.cov(:, :, i), eye(dim));
            end

            newMu = Utility.generateRandomIntMatrix(dim, size_);

            newCov = zeros(dim, dim, size_); % Preallocate
            newCov1 = Utility.generateRandomSPDMatrix(dim);
            newCov2 = Utility.generateRandomSPDMatrix(dim);
            newCov(:, :, 1) = newCov1;
            newCov(:, :, 2) = newCov2;

            % Update
            obj.updateDistributionsParameters(newMu, newCov);

            testCase.verifyEqual(obj.mu, newMu);
            testCase.verifyEqual(obj.cov, newCov);
        end

        function testRemoveDimensions(testCase)
            % Test 1: type = "DS"
            type = "DS";
            dim = 4;
            size_ = 2;
            cols = true;

            obj = GaussianContainer(type, size_, cols, dim);

            testCase.verifyEqual(obj.mu, zeros(dim, size_));
            testCase.verifyEqual(obj.cov, eye(dim));

            newMu = [1, 2; 3, 4; 5, 6; 7, 8];
            newCov = [
                4, 2, 1, 3;
                2, 5, 2, 1;
                1, 2, 3, 0;
                3, 1, 0, 6
            ];
            % Update
            obj.updateDistributionsMu(newMu);
            obj.updateDistributionsCovariance(newCov);

            % Remove dimensions
            obj.removeDimensions([1, 3]);
            testCase.verifyTrue(obj.dim == 2);
            testCase.verifyEqual(obj.mu, [3, 4; 7, 8]);
            testCase.verifyEqual(obj.cov, [5, 1; 1, 6]);
            
            % ---------------------------------------------------------

            % Test 2: type = "DD"
            type = "DD";
            dim = 4;
            size_ = 2;
            cols = true;

            obj = GaussianContainer(type, size_, cols, dim);

            newMu = [1, 2; 3, 4; 5, 6; 7, 8];
            newCov = [
                4, 2, 1, 3;
                2, 5, 2, 1;
                1, 2, 3, 0;
                3, 1, 0, 6
            ];
            newCov = repmat(newCov, 1, 1, size_);
            

            % Update
            obj.updateDistributionsMu(newMu);
            obj.updateDistributionsCovariance(newCov);

            testCase.verifyEqual(obj.mu, newMu);
            testCase.verifyEqual(obj.cov, newCov);

            % Remove dimensions
            obj.removeDimensions([1, 3]);
            testCase.verifyTrue(obj.dim == 2);
            testCase.verifyEqual(obj.mu, [3, 4; 7, 8]);
            testCase.verifyEqual(obj.cov, repmat([5, 1; 1, 6], 1, 1, obj.Size));
        end



        

        %% Private properties
        % This is independent of the type because it is all about 'mu' and
        % type affects 'cov'
        function testSetters(testCase)
            % Test 1: cols = true
            type = "DS";
            dim = 5;
            size_ = 2;
            cols = true;
            mu = Utility.generateRandomIntMatrix(dim, size_);
            cov = Utility.generateRandomSPDMatrix(dim);
            obj = GaussianContainer(type, size_, cols, dim, mu, cov);

            testCase.verifyTrue(isequal(obj.getExpInit(), obj.E));

            % Set init expectation
            newExp = Utility.generateRandomIntMatrix(dim, size_);
            obj.setExpInit(newExp);

            testCase.verifyTrue(isequal(obj.getExpInit(), newExp));
            

            % Test 2: cols = false
            type = "DS";
            dim = 5;
            size_ = 2;
            cols = false;
            mu = Utility.generateRandomIntMatrix(dim, size_);
            cov = Utility.generateRandomSPDMatrix(dim);
            obj = GaussianContainer(type, size_, cols, dim, mu, cov);

            testCase.verifyTrue(isequal(obj.getExpInit(), obj.E));

            % Set init expectation
            newExp = Utility.generateRandomIntMatrix(size_, dim);
            obj.setExpInit(newExp);

            testCase.verifyTrue(isequal(obj.getExpInit(), newExp));
        end


         


        %% Dependent properties
        % Independent of cols
        % obj.H
        function testEntropyProperties(testCase)
            % Test 1: type = "DS"; cov is spherical
            type = "DS";
            dim = 4;
            size_ = 2;
            cols = true;

            obj = GaussianContainer(type, size_, cols, dim);

            testCase.verifyEqual(obj.mu, zeros(dim, size_));
            testCase.verifyEqual(obj.cov, eye(dim));

            testCase.verifyEqual(obj.H, obj.Size * obj.dim/2 * (1 + log(2*pi)));

            % -------------------------------------------------------------------

            % Test 2: type = "DD"
            type = "DD";
            dim = 4;
            size_ = 2;
            cols = true;

            obj = GaussianContainer(type, size_, cols, dim);
            testCase.verifyEqual(obj.mu, zeros(dim, size_));
            testCase.verifyEqual(obj.cov, repmat(eye(dim), 1, 1, size_));

            testCase.verifyEqual(obj.H, obj.Size * obj.dim/2 * (1 + log(2*pi)));

            % Test 3: type = "DS"; cov is full matrix
            type = "DS";
            dim = 4;
            size_ = 2;
            cols = true;

            obj = GaussianContainer(type, size_, cols, dim);

            testCase.verifyEqual(obj.mu, zeros(dim, size_));
            testCase.verifyEqual(obj.cov, eye(dim));

            newCov = Utility.generateRandomSPDMatrix(dim);
            obj.updateDistributionsCovariance(newCov);
            testCase.verifyEqual(obj.cov, newCov);

            H = obj.Size/2 * log(det(newCov)) + obj.Size * obj.dim/2 * (1 + log(2*pi));
            testCase.verifyTrue(abs(obj.H - H) < 1e-12);
        end

        % Test 1: 
        % type: "DS"
        % cols = true
        function testFormatDependentProperties_1(testCase)
            % Setup
            %   cov: [10, 1; 1, 2];
            type = "DS";
            size_ = 3;
            dim = 2;
            cols = true;
            mu = [[0; 0], [1; 1], [1; 2]];

            cov = [10, 1; 1, 2];

            obj = GaussianContainer(type, size_, cols, dim, mu, cov);

            % Test 1
            testCase.verifyEqual(size(obj.E), [dim, obj.Size]);
            testCase.verifyTrue(Utility.isSymmetricMatrix(obj.E_XtX));
            testCase.verifyTrue(Utility.isSymmetricMatrix(obj.E_XXt));

            testCase.verifyEqual(obj.E, mu);
            testCase.verifyEqual(obj.E_Xt, mu');

            testCase.verifyEqual(obj.E_SNC, [12; 14; 17]);
            testCase.verifyEqual(obj.E_XtX, [[12, 0, 0]; [0, 14, 3]; [0, 3, 17]]);
            testCase.verifyEqual(obj.E_XXt, [32, 6; 6, 11]);
            testCase.verifyEqual(obj.E_TrXtX, trace(obj.E_XtX));
        end

        % Test 2: 
        % type: "DS"
        % cols = false
        function testFormatDependentProperties_2(testCase)
            % Setup
            %   cov: [10, 1; 1, 2];
            type = "DS";
            size_ = 3;
            dim = 2;
            cols = false;
            mu = [[0; 0], [1; 1], [1; 2]];

            cov = [10, 1; 1, 2];

            obj = GaussianContainer(type, size_, cols, dim, mu, cov);

            % -------------------------------------------------------
            % Test 1
            testCase.verifyEqual(size(obj.E), [obj.Size, dim]);
            testCase.verifyTrue(Utility.isSymmetricMatrix(obj.E_XtX));
            testCase.verifyTrue(Utility.isSymmetricMatrix(obj.E_XXt));

            testCase.verifyEqual(obj.E, mu');
            testCase.verifyEqual(obj.E_Xt, mu);

            testCase.verifyEqual(obj.E_SNC, [32; 11]);
            testCase.verifyEqual(obj.E_XtX, [[32, 6]; [6, 11]]);
            testCase.verifyEqual(obj.E_XXt, [12, 0, 0; 0, 14, 3; 0, 3, 17]);
            testCase.verifyEqual(obj.E_TrXtX, trace(obj.E_XtX));
        end

        % Test 3: 
        % type: "DD"
        % cols = true
        function testFormatDependentProperties_3(testCase)
            % Setup
            %   dist1: standard normal
            %   dist2: mu = 1, cov: diag([2, 2])
            %   dist3: mu = [1; 2], cov: newCov
            type = "DD";
            size_ = 3;
            dim = 2;
            cols = true;
            mu = [[0; 0], [1; 1], [1; 2]];

            obj = GaussianContainer(type, size_, cols, dim, mu);
           
            cov = repmat(eye(dim), 1, 1, size_);
            cov(:, :, 2) = diag([2, 2]);
            cov(:, :, 3) = [10, 1; 1, 2];

            obj.updateDistributionsCovariance(cov);
            
            % Test 1
            testCase.verifyEqual(size(obj.E), [dim, obj.Size]);
            testCase.verifyTrue(Utility.isSymmetricMatrix(obj.E_XtX));
            testCase.verifyTrue(Utility.isSymmetricMatrix(obj.E_XXt));

            testCase.verifyEqual(obj.E, mu);
            testCase.verifyEqual(obj.E_Xt, mu');

            testCase.verifyEqual(obj.E_SNC, [2; 6; 17]);
            testCase.verifyEqual(obj.E_XtX, [[2, 0, 0]; [0, 6, 3]; [0, 3, 17]]);
            testCase.verifyEqual(obj.E_XXt, [15, 4; 4, 10])
            testCase.verifyEqual(obj.E_TrXtX, trace(obj.E_XtX));

            % Test 2
            % Update mu and cov
            newMu = [[5; 7], [1; 1], [1; 2]];
            obj.updateDistributionsMu(newMu);
            cov(:, :, 1) = [3, 4; 4, 7];
            obj.updateDistributionsCovariance(cov);

            testCase.verifyEqual(size(obj.E), [dim, obj.Size]);
            testCase.verifyTrue(Utility.isSymmetricMatrix(obj.E_XtX));
            testCase.verifyTrue(Utility.isSymmetricMatrix(obj.E_XXt));

            testCase.verifyEqual(obj.E, newMu);
            testCase.verifyEqual(obj.E_Xt, newMu');

            testCase.verifyEqual(obj.E_SNC, [84; 6; 17]);
            testCase.verifyEqual(obj.E_XtX, [[84, 12, 19]; [12, 6, 3]; [19, 3, 17]]);
            testCase.verifyEqual(obj.E_XXt, [42, 43; 43, 65]);
            testCase.verifyEqual(obj.E_TrXtX, trace(obj.E_XtX));
        end

        % Test 4: 
        % type: "DD"
        % cols = false
        function testFormatDependentProperties_4(testCase)
            % Setup
            %   dist1: standard normal
            %   dist2: mu = 1, cov: diag([2, 2])
            %   dist3: mu = [1; 2], cov: newCov
            type = "DD";
            size_ = 3;
            dim = 2;
            cols = false;
            mu = [[0; 0], [1; 1], [1; 2]];

            obj = GaussianContainer(type, size_, cols, dim, mu);
           
            cov = repmat(eye(dim), 1, 1, size_);
            cov(:, :, 2) = diag([2, 2]);
            cov(:, :, 3) = [10, 1; 1, 2];

            obj.updateDistributionsCovariance(cov);

            % Test 1
            testCase.verifyEqual(size(obj.E), [obj.Size, dim]);
            testCase.verifyTrue(Utility.isSymmetricMatrix(obj.E_XtX));
            testCase.verifyTrue(Utility.isSymmetricMatrix(obj.E_XXt));

            testCase.verifyEqual(obj.E, mu');
            testCase.verifyEqual(obj.E_Xt, mu);
 
            testCase.verifyEqual(obj.E_SNC, [15; 10]);
            testCase.verifyEqual(obj.E_XtX, [[15, 4]; [4, 10]]);
            testCase.verifyEqual(obj.E_XXt, [[2, 0, 0]; [0, 6, 3]; [0, 3, 17]]);
            testCase.verifyEqual(obj.E_TrXtX, trace(obj.E_XtX));


            % Test 2
            newMu = [[5; 7], [1; 1], [1; 2]];
            obj.updateDistributionsMu(newMu);
            cov(:, :, 1) = [3, 4; 4, 7];
            obj.updateDistributionsCovariance(cov);
           
            testCase.verifyEqual(size(obj.E), [obj.Size, dim]);
            testCase.verifyTrue(Utility.isSymmetricMatrix(obj.E_XtX));
            testCase.verifyTrue(Utility.isSymmetricMatrix(obj.E_XXt));

            testCase.verifyEqual(obj.E, newMu');
            testCase.verifyEqual(obj.E_Xt, newMu);

            testCase.verifyEqual(obj.E_SNC, [42; 65]);
            testCase.verifyEqual(obj.E_XtX, [[42, 43]; [43, 65]]);
            testCase.verifyEqual(obj.E_XXt, [84, 12, 19; 12, 6, 3; 19, 3, 17]);
            testCase.verifyEqual(obj.E_TrXtX, trace(obj.E_XtX));
        end
    end
end
