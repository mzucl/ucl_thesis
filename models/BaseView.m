% Base class for View
classdef BaseView < handle
    properties
        D
    end

    methods
        function obj = BaseView(varargin)
            disp('BaseView');
            % Constructor for DataView
            % Optionally, you can handle input arguments here
        end
    end
end