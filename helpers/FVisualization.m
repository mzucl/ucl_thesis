classdef FVisualization
    methods (Static, Access=public)
        function cmap = divergingMap(start_hue, end_hue, n)
            % Create a diverging colormap between two hues using HSV interpolation
            % start_hue, end_hue: angles in degrees (e.g., 40 for yellow, 240 for blue)
            % n: number of colors
        
            % Normalize hue values between 0 and 1 for HSV
            h1 = mod(start_hue / 360, 1);
            h2 = mod(end_hue / 360, 1);
        
            % Linearly interpolate in HSV space
            h = linspace(h1, h2, n)';
            s = ones(n, 1);  % full saturation
            v = ones(n, 1);  % full value
            hsv = [h, s, v];
        
            % Convert to RGB
            cmap = hsv2rgb(hsv);
        end

        
        function colors = getColors(values, cmapName)
            % Generate a color mapping for a list of values.
            % values: vector of values
            % cmapName: 'yellow' (default) or a MATLAB colormap name (e.g., 'parula', 'jet', etc.)
        
            if nargin < 3
                cmapName = 'yellow';
            end
        
            values = unique(values); % Ensure unique for mapping
            n = numel(values);
        
            % Determine colormap and normalization
            if strcmp(cmapName, 'yellow')
                % Yellow-like diverging colormap (custom logic)
                cmap = flipud(FVisualization.divergingMap(240, 40, n));  % create diverging palette
            else
                cmap = colormap(cmapName);
                cmap = interp1(linspace(0, 1, size(cmap, 1)), cmap, linspace(0, 1, n));
            end
        
            if min(values) < 0 && max(values) > 0
                normVals = (values - min(values)) / (max(values) - min(values));
            elseif min(values) >= 0
                normVals = (values - min(values)) / (max(values) - min(values));
            elseif max(values) <= 0
                normVals = (values - min(values)) / (max(values) - min(values));
            end
        
            cmapScaled = interp1(linspace(0, 1, size(cmap, 1)), cmap, normVals);
            colors = containers.Map('KeyType', 'double', 'ValueType', 'any');
            for i = 1:n
                colors(values(i)) = cmapScaled(i, :);
            end
        end


        function [cmap, norm, colorMapping] = generateColorMap(Ws)
            % Based on factor weights in Ws, generate a gradient colormap.
            % Ws: cell array of matrices
        
            allValues = sort(cell2mat(cellfun(@(w) w(:), Ws, 'UniformOutput', false)));
            colorMapping = FVisualization.getColors(allValues);
        
            % Extract color values from map
            values = cell2mat(colorMapping.keys);
            colorList = cell2mat(colorMapping.values');
            cmap = colorList;  % Matrix of RGB rows
        
            % Use MATLAB normalization for image display
            norm = struct('vmin', min(values), 'vmax', max(values));
        end


        function plotParams = getFactorsMatrixPlotParams(factors_Ws, d, n_features)
            % Compiles and returns all parameters needed to plot the factors matrix.
            
            % Flatten all the factors and calculate the global min and max values
            allVals = cell2mat(cellfun(@(W) W(:), factors_Ws, 'UniformOutput', false));
            plotParams.vmin = min(min(allVals));
            plotParams.vmax = max(max(allVals));
            
            % Parameters related to plot formatting
            plotParams.axes_sep = cumsum(d(1:end-1)); % Separate modalities based on dimensionalities
            plotParams.label_size = 10;
            plotParams.top_margin = 0.04;
            plotParams.bottom_margin = 0.04;
        
            % Calculate figure height in inches based on the number of features
            dpi = 72.27; % Dots per inch for plotting
            matrix_height_pt = plotParams.label_size * n_features; % Height in points
            matrix_height_in = matrix_height_pt / dpi; % Height in inches
            
            % Calculate the total figure height adjusted by margins
            plotParams.figure_height = matrix_height_in / (1 - plotParams.top_margin - plotParams.bottom_margin);
        end

        % function plotParams = getFactorsMatrixPlotParams(factorsWs, d, nFeatures)
        %     % Compile parameters for plotting stable factor matrix
        % 
        %     allData = cell2mat(cellfun(@(w) w(:), factorsWs, 'UniformOutput', false));
        %     plotParams.vmin = min(allData);
        %     plotParams.vmax = max(allData);
        %     plotParams.axes_sep = cumsum(d(1:end-1));
        %     plotParams.label_size = 10;
        %     plotParams.top_margin = 0.04;
        %     plotParams.bottom_margin = 0.04;
        % 
        %     dpi = 72.27;
        %     matrixHeightPt = plotParams.label_size * nFeatures;
        %     matrixHeightIn = matrixHeightPt / dpi;
        % 
        %     plotParams.figure_height = matrixHeightIn / ...
        %         (1 - plotParams.top_margin - plotParams.bottom_margin);
        % end

        function plotStableFactorsMatrix(factorsWs, d, featureLabels, saveFig)
            % Function that plots the stable factors matrix (heatmap) for each stable factor.
            % 
            % Inputs:
            %   - factorsWs: cell array where each cell is a matrix of size [n_features x n_models]
            %   - d: array of modality dimensionalities
            %   - featureLabels: cell array of strings for feature labels
            %   - saveFig: (optional) boolean, whether to save the figure (default: false)
        
            if nargin < 5
                saveFig = false;
            end
        
            % Generate colormap and normalization values
            [cmap, norm, ~] = FVisualization.generateColorMap(factorsWs);
        
            % Get plotting parameters
            plotParams = FVisualization.getFactorsMatrixPlotParams(factorsWs, d, numel(featureLabels));
        
            % Setup figure
            nFactors = numel(factorsWs);
            figure('Units', 'normalized', 'Position', [0.1 0.1 0.1 + 0.2 * nFactors, 0.6]);
            t = tiledlayout(1, nFactors, ...
                'TileSpacing', 'none', ...
                'Padding', 'none');
        
            for i = 1:nFactors
                nexttile;
        
                factor = factorsWs{i};
                imagesc(factor');  % transpose to match orientation
                colormap(cmap);
                caxis([plotParams.vmin, plotParams.vmax]);
        
                if i == 1
                    set(gca, 'YTick', 1:numel(featureLabels), 'YTickLabel', featureLabels, ...
                        'FontSize', plotParams.label_size);
                else
                    set(gca, 'YTick', [], 'FontSize', plotParams.label_size);
                end
                set(gca, 'XTick', []);
        
                if i == nFactors
                    colorbar;
                end
        
                % Draw horizontal lines to separate modalities
                hold on;
                for sep = plotParams.axes_sep
                    yline(sep + 0.5, 'k', 'LineWidth', 1.2);  % yline at boundary between modalities
                end
                hold off;
            end
        
            % Set global labels
            xlabel(t, 'Stable factor', 'FontSize', 3 * plotParams.label_size);
            ylabel(t, 'Feature', 'FontSize', 3 * plotParams.label_size);
            title(t, 'Stable factors', 'FontSize', 3 * plotParams.label_size);
        
            % Save or show
            % exportgraphics(gcf, fullfile(obj.plots_path, 'stable_factors_matrix.svg'), 'ContentType', 'vector');
            if saveFig
                exportgraphics(gcf, fullfile('plots', 'stable_factors_matrix.svg'), 'ContentType', 'vector');
                exportgraphics(gcf, fullfile('plots', 'stable_factors_matrix.jpg'), 'ContentType', 'image');
                close;
            end
        end


    end
end
