classdef BoundsTest < matlab.unittest.TestCase
    methods (Static)
        function verifyObject(testCase, bS, bM)
            % bS - bound when argument is a scalar
            % bM - bound when argument is a matrix
            c = bM.c();
            g = bM.g();
            h = bM.h();

            testCase.verifyTrue(all(c(:) == bS.c()));
            testCase.verifyTrue(all(g(:) == bS.g()));
            testCase.verifyTrue(all(h(:) == bS.h()));
        end
    end

    methods(Test)
        function testLaplaceApproximation(testCase)
            N = 10;
            a = randi(15);
            A = a * ones(N, N);
            BoundsTest.verifyObject(testCase, LaplaceApproximation(a), ...
                 LaplaceApproximation(A));
        end

        function testJaakkolaBound(testCase)
            N = 10;
            a = randi(15);
            A = a * ones(N, N);
            BoundsTest.verifyObject(testCase, JaakkolaBound(a), ...
                 JaakkolaBound(A));
        end

        function testBohningBound(testCase)
            N = 10;
            a = randi(15);
            A = a * ones(N, N);
            BoundsTest.verifyObject(testCase, BohningBound(a), ...
                 BohningBound(A));
        end
    end
end
  