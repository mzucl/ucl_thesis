classdef ViewHandlerTest < matlab.unittest.TestCase
    methods (Test)
        function testConstructor(testCase)
            data = [
                [1, 2, 3, 4, 5]; 
                [2, 1, 2, 1, 2];
                [1, 3, 5, 7, 9]
                ];
            
            % Test 1: featuresInCols = true
            vh = ViewHandler(data);

            testCase.verifyEqual(vh.X, data');

            % Test 2: featuresInCols = false
            vh = ViewHandler(data, false);

            testCase.verifyEqual(vh.X, data);
        end



        %% Getters/Retrieval methods
        function testRetrievalMethods(testCase)
            data = [
                [1, 2, 3, 4, 5]; 
                [2, 1, 2, 1, 2];
                [1, 3, 5, 7, 9]
                ];
            
            % Test 1: featuresInCols = true
            vh = ViewHandler(data);
            idx = 2;
            testCase.verifyEqual(vh.getObservation(idx), [2; 1; 2; 1; 2]);
            testCase.verifyEqual(vh.getObservation(idx, true), [2 1 2 1 2]);
            testCase.verifyEqual(vh.getRow(idx), [2 1 3]);
            testCase.verifyEqual(vh.getRow(idx, true), [2; 1; 3]);
            testCase.verifyEqual(vh.getObservationNormSq(idx), 14);
            testCase.verifyEqual(vh.getObservationEntry(idx + 1, 3), 5);
            
            % Test 2: featuresInCols = false
            vh = ViewHandler(data, false);
            idx = 4;
            testCase.verifyEqual(vh.getObservation(idx), [4; 1; 7]);
            testCase.verifyEqual(vh.getObservation(idx, true), [4, 1, 7]);
            testCase.verifyEqual(vh.getRow(idx - 1), [1 3 5 7 9]);
            testCase.verifyEqual(vh.getRow(idx - 1, true), [1; 3; 5; 7; 9]);
            testCase.verifyEqual(vh.getObservationNormSq(idx), 66);
            testCase.verifyEqual(vh.getObservationEntry(idx + 1, 3), 9);
        end

        function testGetters(testCase)
            data = [
                [1, 2, 3, 4, 5]; 
                [2, 1, 2, 1, 2];
                [1, 3, 5, 7, 9]
                ];
            
            % Test 1: featuresInCols = true
            vh = ViewHandler(data);
            testCase.verifyEqual(vh.N, 3);
            testCase.verifyEqual(vh.D, 5);
            testCase.verifyEqual(vh.Tr_XtX, trace(data * data'));

            % Test 2: featuresInCols = false
            vh = ViewHandler(data, false);
            testCase.verifyEqual(vh.N, 5);
            testCase.verifyEqual(vh.D, 3);
            testCase.verifyEqual(vh.Tr_XtX, trace(data' * data));
        end
    end
end
