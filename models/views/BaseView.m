% Base class for View
classdef BaseView < handle
    properties
        D
    end

    methods
        function obj = BaseView(varargin)
        end

        function obj = qAlphaUpdate(obj)
            bNew = obj.alpha.prior.b + 1/2 * obj.W.E_SNC;
            obj.alpha.updateAllDistributionsB(bNew);
        end
    end
end