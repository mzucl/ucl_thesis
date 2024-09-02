classdef UtilityTest < matlab.unittest.TestCase
    methods (Test)
        function testIsSingleNumber(testCase)
            testCase.verifyTrue(Utility.isSingleNumber(5));
            testCase.verifyTrue(~Utility.isSingleNumber(NaN));
            testCase.verifyTrue(~Utility.isSingleNumber([1, 2]));
            testCase.verifyTrue(Utility.isSingleNumber(5));
            testCase.verifyTrue(~Utility.isSingleNumber([1, 2; 2, 3]));
            testCase.verifyTrue(~Utility.isSingleNumber( Gamma() ));
        end

        function testIsMonotonicIncreasing(testCase)
            testCase.verifyTrue(Utility.isMonotonicIncreasing([1, 2, 3]));
            testCase.verifyTrue(Utility.isMonotonicIncreasing([-1.23, -0.56, 3]));
            testCase.verifyTrue(~Utility.isMonotonicIncreasing([-1.23, -0.56, 0, 1, 5, 2]));
            testCase.verifyTrue(~Utility.isMonotonicIncreasing([1.23, 0.56, 0.11]));
        end

        function testAreEqual(testCase)
            testCase.verifyTrue(Utility.areEqual(NaN, NaN));
            testCase.verifyTrue(~Utility.areEqual(Gamma(), NaN));
            testCase.verifyTrue(~Utility.areEqual(NaN, Gamma()));
            testCase.verifyTrue(Utility.areEqual(Gamma(1, 2), Gamma(1, 2)));
            testCase.verifyTrue(~Utility.areEqual(Gamma(1, 2), Gamma(3, 2)));
            testCase.verifyTrue(~Utility.areEqual(Gamma(1, 2), Gaussian())); % Different objects
        end
    end
end
