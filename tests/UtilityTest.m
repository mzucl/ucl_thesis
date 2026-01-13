classdef UtilityTest < matlab.unittest.TestCase
    methods (Test)
        function testIsSingleNumber(testCase)
            testCase.verifyTrue(NumericValidation.isFiniteNumericScalar(5));
            testCase.verifyTrue(~NumericValidation.isFiniteNumericScalar(NaN));
            testCase.verifyTrue(~NumericValidation.isFiniteNumericScalar([1, 2]));
            testCase.verifyTrue(NumericValidation.isFiniteNumericScalar(5));
            testCase.verifyTrue(~NumericValidation.isFiniteNumericScalar([1, 2; 2, 3]));
            testCase.verifyTrue(~NumericValidation.isFiniteNumericScalar( Gamma() ));
        end

        function testIsMonotonicIncreasing(testCase)
            testCase.verifyTrue(MatrixValidation.isMonotonicIncreasing([1, 2, 3]));
            testCase.verifyTrue(MatrixValidation.isMonotonicIncreasing([-1.23, -0.56, 3]));
            testCase.verifyTrue(~MatrixValidation.isMonotonicIncreasing([-1.23, -0.56, 0, 1, 5, 2]));
            testCase.verifyTrue(~MatrixValidation.isMonotonicIncreasing([1.23, 0.56, 0.11]));
        end

        function testAreEqual(testCase)
            testCase.verifyTrue(TypeValidation.areEquivalent(NaN, NaN));
            testCase.verifyTrue(~TypeValidation.areEquivalent(Gamma(), NaN));
            testCase.verifyTrue(~TypeValidation.areEquivalent(NaN, Gamma()));
            testCase.verifyTrue(TypeValidation.areEquivalent(Gamma(1, 2), Gamma(1, 2)));
            testCase.verifyTrue(~TypeValidation.areEquivalent(Gamma(1, 2), Gamma(3, 2)));
            testCase.verifyTrue(~TypeValidation.areEquivalent(Gamma(1, 2), Gaussian())); % Different objects
        end
    end
end
