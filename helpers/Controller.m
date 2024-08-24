classdef Controller
    properties (Access = private)
        dirtyFlags % Store flags for each dependent property
    end

    methods
        function obj = Controller()
            obj.dirtyFlags = containers.Map();
        end

        function setDirtyFlag(obj, propertyName)
            obj.dirtyFlags(propertyName) = true;
        end

        function clearDirtyFlag(obj, propertyName)
            obj.dirtyFlags(propertyName) = false;
        end

        function isDirty = isDirty(obj, propertyName)
            if isKey(obj.dirtyFlags, propertyName)
                isDirty = obj.dirtyFlags(propertyName);
            else
                isDirty = true; % Assume dirty if not initialized
            end
        end
    end
end
