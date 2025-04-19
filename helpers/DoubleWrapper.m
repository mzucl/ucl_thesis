% A wrapper class to enable passing double value by reference
classdef DoubleWrapper < handle
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