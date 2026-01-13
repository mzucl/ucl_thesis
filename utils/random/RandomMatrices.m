classdef RandomMatrices
    % RANDOMMATRICES Utility class to generate random matrices
    %
    % Provides static methods to generate:
    %   - Binary matrices
    %   - Random integer matrices
    %   - Symmetric positive definite (SPD) matrices
    %   - Random rotation matrices
    %
    % Example usage:
    %   A = RandomMatrices.binaryMatrix(3,4);
    %   B = RandomMatrices.intMatrix(5,5,1,20);
    %   C = RandomMatrices.spdMatrix(4);
    %   R = RandomMatrices.rotationMatrix(3);

    methods(Static)
        %% Binary matrix
        function matrix = binaryMatrix(m, n)
            % BINARYMATRIX Generate m-by-n matrix of random 0s and 1s
            arguments
                m (1,1) {mustBeInteger, mustBePositive}
                n (1,1) {mustBeInteger, mustBePositive} = m
            end
            matrix = randi([0, 1], m, n);
        end

        %% Random integer matrix
        function matrix = intMatrix(m, n, minValue, maxValue)
            % INTMATRIX Generate m-by-n matrix of random integers
            %
            % Usage:
            %   intMatrix(m)            -> m-by-m matrix, default 1..10
            %   intMatrix(m,n)          -> m-by-n matrix, default 1..10
            %   intMatrix(m,n,min,max)  -> specify range

            arguments
                m (1,1) {mustBeInteger, mustBePositive}
                n (1,1) {mustBeInteger, mustBePositive} = m
                minValue (1,1) {mustBeInteger} = 1
                maxValue (1,1) {mustBeInteger} = 10
            end

            matrix = randi([minValue, maxValue], m, n);
        end

        %% Symmetric positive definite matrix
        function matrix = spdMatrix(n)
            % SPDMATRIX Generate n-by-n symmetric positive definite matrix
            arguments
                n (1,1) {mustBeInteger, mustBePositive}
            end

            R = randn(n);
            matrix = R' * R;
        end

        %% Random rotation matrix
        function R = rotationMatrix(n)
            % ROTATIONMATRIX Generate a random n-by-n rotation matrix (det = 1)
            % Efficiently ensures det(R) = 1 without computing full determinant
            arguments
                n (1,1) {mustBeInteger, mustBePositive}
            end
        
            % [Q, Rtri] = qr(Drandn(n));
            [Q, ~] = qr(Drandn(n));
            
            % Determine sign of determinant from upper-triangular R
            % TODO: Check if the condition below is a valid (the goal is to
            % avoid computing det(Q))
            % if prod(sign(diag(Rtri))) < 0
            if det(Q) < 0
                Q(:,1) = -Q(:,1);  % Flip first column to correct determinant
            end
        
            R = Q;
        end
    end
end