classdef DoubleWrapper < handle
    % DoubleWrapper
    % A handle class that wraps a double value, enabling pass-by-reference semantics.
    %
    % Usage:
    %   w = DoubleWrapper(5);      % Create wrapper with initial value 5

    properties
        Val
    end
    
    methods
        function obj = DoubleWrapper(val)
            if nargin > 0
                obj.Val = val;
            end
        end
        
        function setValue(obj, newValue)
            obj.Val = newValue;
        end
    end
end