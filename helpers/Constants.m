% Constants are defined on the level of the framework
classdef Constants
    properties (Constant)
        %% Visualization.m
        BLUE      = '#DAE8FC';
        PURPLE    = '#E1D5E7';
        RED       = '#F8CECC';
        DARK_BLUE = '#6C8EBF';
        GREEN     = '#D5E8D4';

        % Export
        FIGURES_FOLDER = 'figures';
        EXPORT_TO_PDF  = false;



        %% General
        % Tolerance for comparison of floating point values
        % if(a == b) -> if(abs(a - b) < Constants.TOL)
        TOL = 1e-9;  
    end
end