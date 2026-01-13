classdef Visualization
    methods (Static, Access=private)
        function formatFigure(hfig)
            %formatFigure Formats a MATLAB figure for publication-quality appearance.
            %
            % Description:
            % This function applies a consistent set of formatting options to the 
            % specified MATLAB figure handle. It adjusts font size, removes boxes 
            % around axes, and sets the text interpreters for LaTeX-style rendering.
            %
            % Syntax:
            %   formatFigure(hfig)
            %
            % Inputs:
            %   hfig - Handle to a MATLAB figure. The figure to be formatted.
            %
            % Outputs:
            %   None. The formatting is applied directly to the provided figure handle.
            CustomError.validateNumberOfParameters(nargin, 1, 1);

            set(findall(hfig, '-property','FontSize'),'FontSize', 17);
            set(findall(hfig, '-property', 'Box'), 'Box', 'off');
            set(findall(hfig, '-property', 'Interpreter'), 'Interpreter', 'latex');
            set(findall(hfig, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex');
        end



        function rgb = hexToRGB(hex)
            %hexToRGB - Converts a hexadecimal color code to an RGB vector.
            %
            % Description:
            %   This function takes a hexadecimal color code (e.g., '#FF5733') as input 
            %   and converts it into an RGB vector with values normalized between 0 and 1.
            %
            % Input:
            %   hex - A string representing the hexadecimal color code. The input must 
            %         start with a '#' followed by six hexadecimal digits (e.g., '#RRGGBB').
            %
            % Output:
            %   rgb - A 1x3 vector containing the red, green, and blue color components 
            %         of the input hex code. Each component is a value between 0 and 1.
            CustomError.validateNumberOfParameters(nargin, 1, 1);

            rgb = sscanf(hex(2:end),'%2x%2x%2x',[1 3]) / 255;
        end



        function exportFigure(hfig, figName, subfolderName)
            %exportFigure - Saves or exports a MATLAB figure to a file.
            %
            % Description:
            %   This function saves the specified MATLAB figure to a file with the 
            %   given name. Optionally, the figure can be saved to a specific subfolder.
            %   The figure is exported as a high-quality PNG image. If `config.EXPORT_TO_PDF`
            %   is `true`, it is also exported as a PDF.
            %
            % Input:
            %   hfig         - Handle to the figure to be saved.
            %   figName      - Name of the output file. If empty, a timestamp is used.
            %   subfolderName (optional) - A string specifying the subfolder where 
            %                  the figure should be saved. If not provided or empty, the figure 
            %                  is saved in `config.FIGURES_FOLDER` directory.
            %
            % Output:
            %   None. The function saves the figure to a file.
            CustomError.validateNumberOfParameters(nargin, 1, 3);

            if nargin < 2 || isempty(figName)
                figName = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
            end

            isSubfolderSpecified = false;
            if nargin > 2 && ~isempty(subfolderName)
                isSubfolderSpecified = true;
            end
            

            width = 20;
            ratio = 0.65; % height/weight ratio

            set(hfig, 'Units', 'centimeters', ...
                      'Position', [3 3 width ratio * width]);
            
            pos = get(hfig, 'Position');
            
            set(hfig, 'PaperPositionMode', 'Auto', ...
                      'PaperUnits', 'centimeters', ...
                      'PaperSize', [pos(3), pos(4)]);

            folderName = ConfigUtils.getValue('Export', 'FIGURES_FOLDER');

            % If subfolder name is specified
            if isSubfolderSpecified
                folderName = [folderName, '/', subfolderName];
            end

            % Save figure
            if ~exist(folderName, 'dir')
                mkdir(folderName);
            end

            % Export to PNG
            figNamePNG = [figName, '.png'];
            filePath = fullfile(folderName, figNamePNG);
            set(gcf, 'PaperPositionMode', 'auto');
            exportgraphics(hfig, filePath, 'Resolution', 300);

            % Export to PDF
            if ConfigUtils.getValue('Export', 'EXPORT_TO_PDF')
                figNamePDF = [figName, '.pdf'];
                filePath = fullfile(folderName, figNamePDF);
                set(gcf, 'PaperPositionMode', 'auto');
                exportgraphics(hfig, filePath, 'ContentType', 'vector');
            end
        end



        function renderHintonDiagram(matrix, ax, figTitle, backgroundColor)
            % renderHintonDiagram - Generates a Hinton diagram to visualize matrix values.
            %
            % Description:
            %   This function creates a Hinton diagram to visualize the values of 
            %   the input matrix. The size of each square in the diagram is 
            %   proportional to the magnitude of the corresponding matrix value, 
            %   with positive values represented by white squares and negative 
            %   values by black squares.
            %
            % Input:
            %   matrix         - A numeric matrix whose values will be visualized 
            %                    in the Hinton diagram.
            %   ax (optional)  - Handle to the axis on which the diagram will be drawn. 
            %                    If not provided, the current axis (gca) is used.
            %   figTitle (optional) - A string specifying the title of the figure. 
            %                         Defaults to an empty string if not provided.
            %   backgroundColor (optional) - A boolean value indicating whether 
            %                                to set the background color of the axis. 
            %                                Defaults to true, which sets the background 
            %                                color to a predefined (in config.txt file) blue.
            %
            % Output:
            %   None. The function visualizes the matrix values as a Hinton diagram
            %   on the specified or current axis.
            CustomError.validateNumberOfParameters(nargin, 1, 4);

            % Set default values for optional parameters
            if nargin < 4
                backgroundColor = true;
                if nargin < 3
                    figTitle = '';
                    if nargin < 2
                        ax = gca;
                    end
                end
            end
            
            % Set the current axis to ax
            axes(ax);
            axis(ax, 'equal'); % Set the same length in every direction
            
            if backgroundColor
                set(ax, 'Color', Visualization.hexToRGB(ConfigUtils.getValue('Colors', 'BLUE')));
            end

            maxWeight = max(abs(matrix(:)));

            for i = 1:size(matrix, 1) % y
                for j = 1:size(matrix, 2) % x
                    % Determine the size of the square
                    weight = matrix(i, j);
                    height = sqrt(abs(weight) / maxWeight);
                    width = height;
        
                    % White for positive, black for negative
                    color = [1 1 1] * (weight >= 0);
        
                    rectangle('Position', [j - width / 2, i - height / 2, width, height], 'FaceColor', color, 'EdgeColor', 'none');
                end
            end

            ylim(ax, [0, size(matrix, 1) + 1]);
            xlim(ax, [0, size(matrix, 2) + 1]);

            set(ax, 'YDir', 'reverse', 'XAxisLocation', 'top');            
            % axis(ax, 'off');
            set(ax, 'XColor', 'none', 'YColor', 'none');

            if ~isempty(figTitle)
                title(ax, figTitle, 'Interpreter', 'latex');
            end
        end



        function renderFactorLoadings(W, dimList, labelPos, ax, figTitle)
            CustomError.validateNumberOfParameters(nargin, 2, 5);
            
            % Set default values for optional parameters
            if nargin < 5
                figTitle = '';
                if nargin < 4
                    ax = gca;
                    if nargin < 3
                        labelPos = 'top';
                    end
                end
            end
            
            % Set the current axis to ax
            axes(ax);

            [D, K] = size(W);
        
            % Max height per each row
            maxHeight = 0;
            for k = 1:K
                h = max(W(:, k)) - min(W(:, k));
                if h > maxHeight
                    maxHeight = h;
                end
            end
            offset = ceil(1.25 * maxHeight);
            x = linspace(1, D, D);

            for k = 1:K
                y = W(:, k) + (K - k + 1) * offset;
                plot(x, y, 'Color', 'black', 'LineWidth', 1.5); %, 'Label', ['ind: ', k]);
                hold on;
            end
            
            ax.XAxis.Visible = 'off';
            ax.YAxis.Visible = 'off';
            
            linePos = 0;
            for k = 1:length(dimList) - 1 % Don't plot vertical line for the last view!
                linePos = linePos + dimList(k);
                labelText = ['$D_{', num2str(k), '}$'];
                xline(linePos, 'Color', ConfigUtils.getValue('Colors', 'DARK_BLUE'), 'LineWidth', 2, ...
                    'Label', labelText, 'LabelVerticalAlignment', labelPos, ...
                    'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');
            end

            if ~isempty(figTitle)
                title(ax, figTitle, 'Interpreter', 'latex');
            end
        end

        
    end



    methods (Static, Access=public)
        function plotHintonDiagrams(arrW, subplotTitles, figTitle, figName, subfolderName)
            % plotHintonDiagrams - Plots Hinton diagrams for multiple matrices side by side.
            %
            % Description:
            %   This function generates Hinton diagrams for each matrix in the input array 
            %   `arrW`. The diagrams are displayed side by side in a single figure. Optionally, 
            %   if a `figName` is provided, the figure will be saved to a file with the specified name. 
            %   The function also allows specifying a subfolder for saving the figure, using the 
            %   `subfolderName` argument. 
            %
            %   `subplotTitles` provides the titles for each subplot. If not enough titles are provided, 
            %   only the first subplots will have titles, and the remaining subplots will not have titles. 
            %   If more titles are provided than there are subplots, only the first `number of subplots` 
            %   will be used. The main figure title (`figTitle`) will be applied only if more than one subplot is created.
            %
            % Input:
            %   - arrW (required): A cell array of matrices to be visualized as Hinton diagrams.
            %   - subplotTitles (optional): A cell array of titles for each subplot. 
            %   - figTitle (optional): The title for the entire figure. It will be applied only 
            %                          if multiple subplots are created.
            %   - figName (optional): The name of the file to which the figure will be exported. 
            %                         If not provided, the figure will not be saved.
            %                         If '' provided, current date and time
            %                         will be used.
            %   - subfolderName (optional): The subfolder where the figure will be saved. 
            %                                If not provided, the figure will be saved in the current folder.
            %
            % Output:
            %   - None. The function visualizes the Hinton diagrams on the screen and optionally exports the figure.
            CustomError.validateNumberOfParameters(nargin, 1, 5);
            
            isSubplotTitlesSpecified = false;
            addTitle = false;
            saveFig = false;
            if nargin > 1
                if ~isempty(subplotTitles)
                    isSubplotTitlesSpecified = true;
                end
                if nargin > 2
                    if ~isempty(figTitle)
                        addTitle = true;
                    end
                    if nargin > 3
                        saveFig = true;
                        if nargin == 4 % `subfolderName` is not provided!
                            subfolderName = '';
                        end
                    end
                end
            end
            

            % Create figure
            hfig = figure;

            numPlots = length(arrW);

            hfigAx = axes(hfig); % Moved here for optimization purposes
            for i = 1:numPlots
                ax = LogicUtils.ternary(numPlots == 1, hfigAx, subplot(1, numPlots, i));
                Visualization.renderHintonDiagram(arrW{i}, ax, ...
                    LogicUtils.ternaryOpt(isSubplotTitlesSpecified && i <= length(subplotTitles), @() subplotTitles{i}, @() ''));
            end

            if addTitle && numPlots > 1
                sgtitle(figTitle);
            end
           
            % Format figure
            Visualization.formatFigure(hfig);
            
            % Save figure
            if saveFig
                Visualization.exportFigure(hfig, figName, subfolderName);
            end
        end



        % [NOTE] The number of inferred latent factors may differ from the true
        % number of latent factors. As a result, they are not plotted as subplots
        % in a single figure.
        function plotLatentFactors(Z, figTitle, figName, subfolderName)
            % plotLatentFactors - Plots each latent factor from a matrix Z as a separate subplot.
            %
            % Description:
            %   This function visualizes the rows of the matrix `Z` as individual time series, 
            %   plotting each row (latent factor) in a separate subplot stacked vertically. 
            %   An optional figure title can be added using `figTitle`, and the figure can be 
            %   exported to a file using `figName`. If a `subfolderName` is specified, the figure 
            %   will be saved in that subfolder; otherwise, it will be saved in the `config.FIGURES_FOLDER`.
            %   If `figName` is provided but it is empty, the figure will be saved 
            %   using the current timestamp as its name.
            %
            % Input:
            %   - Z (required): A matrix of size [K x N], where each row corresponds to a latent factor.
            %   - figTitle (optional): A string specifying the overall figure title. If empty or not provided, no title is added.
            %   - figName (optional): The name of the file to save the figure as. If empty, a timestamp is used.
            %   - subfolderName (optional): The subfolder where the figure will be saved.
            %
            % Output:
            %   - None. The function displays the latent factors in a figure and optionally saves it to a file.
            CustomError.validateNumberOfParameters(nargin, 1, 4);

            addTitle = false;
            saveFig = false;
            if nargin > 1
                if ~isempty(figTitle)
                    addTitle = true;
                end
                if nargin > 2
                    saveFig = true;
                    if nargin == 3 % `subfolderName` is not specified
                        subfolderName = '';
                    end
                end
            end

            hfig = figure;
            numFactors = size(Z, 1); % Z = [K x N], where K is the number of latent factors!
            for i = 1:numFactors
                subplot(numFactors, 1, i);
                factor = Z(i, :);
                plot(factor, '.', 'MarkerSize', 4);
                hold on;
            end

            if addTitle
                sgtitle(figTitle);
            end

            Visualization.formatFigure(hfig);

            % Save figure
            if saveFig
                Visualization.exportFigure(hfig, figName, subfolderName);
            end
        end



        function plotFactorLoadings(W, dimList, figTitle, figName, subfolderName)
            % plotLoadings - Visualizes the factor loadings for all views as a heatmap.
            %
            % Description:
            %   This function plots the transpose of the factor loading matrix `W`, such that
            %   each row corresponds to a single latent factor and each column corresponds to 
            %   a feature. The matrix `W` is assumed to concatenate 
            %   loadings from multiple views, and the boundaries between views are indicated 
            %   with vertical lines based on the `dimList` input.
            %
            %   An optional title can be added using `figTitle`, and the figure can be saved 
            %   by specifying `figName`. If `subfolderName` is also provided, the figure 
            %   will be saved to the specified subfolder; otherwise, it will be saved in 
            %   `config.FIGURES_FOLDER` directory.
            %
            % Input:
            %   - W (required): A matrix of size [D_total x K], where D_total is the sum of 
            %                  feature dimensions across all views and K is the number of latent factors.
            %   - dimList (required): A vector specifying the number of features in each view. 
            %                         Used to draw vertical lines that separate views in the plot.
            %   - figTitle (optional): Title for the entire figure. If empty or not provided, 
            %                          no title is added.
            %   - figName (optional): Name of the file to which the figure will be saved. 
            %                         If not provided, the figure will not be saved. 
            %                         If provided as an empty value, a timestamp will be used as the filename.
            %   - subfolderName (optional): Subfolder in which the figure will be saved.
            %   - ax (optional):      Handle to the axis on which the diagram will be drawn. 
            %                         If not provided, the current axis (gca) is used.
            % Output:
            %   - None.
            CustomError.validateNumberOfParameters(nargin, 2, 5);
            
            addTitle = false;
            saveFig = false;
            if nargin > 2
                if ~isempty(figTitle)
                    addTitle = true;
                end
                if nargin > 3
                    saveFig = true;
                    if nargin == 4 % `subfolderName` is not provided!
                        subfolderName = '';
                    end
                end
            end

            % Create figure
            hfig = figure;

            % [NOTE] We don't have subplots here, so we don't pass
            % `figTitle` to `renderFactorLoadings`, we simple add it below.
            Visualization.renderFactorLoadings(W, dimList);

            if addTitle
                title(figTitle, 'Interpreter', 'latex');
            end

            % Format figure
            Visualization.formatFigure(hfig);
            
            % Save figure
            if saveFig
                Visualization.exportFigure(hfig, figName, subfolderName);
            end
        end
    

   
        % [NOTE] This function can be extended with `subplotTitles` param.
        function plotFactorLoadingsAndAlpha(W, dimList, alpha, labelPos, height, widthToHeightRatio, ...
                figTitle, figName, subfolderName)
            % Parameters
            % ----------
            % W : matrix, [D_total x K]
            % dimList: number of features in each view 
            CustomError.validateNumberOfParameters(nargin, 3, 9);

            % Assign default values to optional parameters if they are missing or empty ('')
            defaultHeight = 400;
            defaultWidthToHeightRatio = 2.5;
            defaultLabelPos = 'top';
            if nargin < 6
                widthToHeightRatio = defaultWidthToHeightRatio; % Optimal value for 2 views (at least for synth data)
                if nargin < 5
                    height = defaultHeight;
                    if nargin < 4
                        labelPos = defaultLabelPos;
                    end
                end
            end
            
            addTitle = false;
            saveFig = false;
            if nargin > 6
                if isempty(labelPos)
                    labelPos = defaultLabelPos;
                end
                if isempty(height)
                    height = defaultHeight;
                end
                if isempty(widthToHeightRatio)
                    widthToHeightRatio = defaultWidthToHeightRatio;
                end

                if ~isempty(figTitle)
                    addTitle = true;
                end
                if nargin > 7
                    saveFig = true;
                    if nargin == 8 % `subfolderName` is not provided!
                        subfolderName = '';
                    end
                end
            end

            % Create figure
            hfig = figure('Position', [100, 100, widthToHeightRatio * height, height]);
            subplot('Position', [0.1, 0.1, 0.6, 0.8]);
            Visualization.renderFactorLoadings(W, dimList, labelPos);
            
            subplot('Position', [0.75, 0.1, 0.2, 0.8]);
            Visualization.renderHintonDiagram(alpha);

            if addTitle
                sgtitle(figTitle, 'Interpreter', 'latex');
            end

            % Format figure
            Visualization.formatFigure(hfig);
            
            % Save figure
            if saveFig
                Visualization.exportFigure(hfig, figName, subfolderName);
            end
        end
    end
end
