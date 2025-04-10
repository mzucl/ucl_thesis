classdef Visualization
    methods (Static)
        function formatFigure(hfig)
        % formatFigure Formats a MATLAB figure for publication-quality appearance.
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
            set(findall(hfig, '-property','FontSize'),'FontSize', 17);
            set(findall(hfig, '-property', 'Box'), 'Box', 'off');
            set(findall(hfig, '-property', 'Interpreter'), 'Interpreter', 'latex');
            set(findall(hfig, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex');
        end



        function rgb = hexToRGB(hex)
        % hexToRGB - Converts a hexadecimal color code to an RGB vector.
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
            rgb = sscanf(hex(2:end),'%2x%2x%2x',[1 3]) / 255;
        end



        function exportFigure(hfig, figName, subfolderName)
            % exportFigure - Saves or exports a MATLAB figure to a file.
            %
            % Description:
            %   This function saves the specified MATLAB figure to a file with the 
            %   given name. Optionally, the figure can be saved to a specific subfolder.
            %   The figure is exported as a high-quality PNG image.
            %
            % Input:
            %   hfig         - Handle to the figure to be saved.
            %   figName      - A string specifying the name of the output file 
            %                  (e.g., 'myfigure').
            %   subfolderName (optional) - A string specifying the subfolder where 
            %                  the figure should be saved. If not provided, the figure 
            %                  is saved in the current working directory.
            %
            % Output:
            %   None. The function saves the figure to a file.
            if nargin < 2
                CustomError.raiseError('InputCheck', CustomError.ERR_NOT_ENOUGH_INPUT_ARG);
            elseif nargin > 3
                CustomError.raiseError('InputCheck', CustomError.ERR_TOO_MANY_INPUT_ARG);
            end
            
            width = 20;
            ratio = 0.65; % height/weight rati

            set(hfig, 'Units', 'centimeters', ...
                      'Position', [3 3 width ratio * width]);
            
            pos = get(hfig, 'Position');
            
            set(hfig, 'PaperPositionMode', 'Auto', ...
                      'PaperUnits', 'centimeters', ...
                      'PaperSize', [pos(3), pos(4)]);

            folderName = Constants.FIGURES_FOLDER;

            % If a subfolder name is specified
            if nargin > 2 && ~isempty(subfolderName)
                folderName = [Constants.FIGURES_FOLDER, '/', subfolderName];
            end

            % Save figure
            if ~exist(folderName, 'dir')
                mkdir(folderName);
            end

            % Export to .pdf
            if Constants.EXPORT_TO_PDF
                figNamePDF = [figName, '.pdf'];
                filePath = fullfile(folderName, figNamePDF);
                set(gcf, 'PaperPositionMode', 'auto');
                exportgraphics(hfig, filePath, 'ContentType', 'vector');
            end

            % Export to .png
            figNamePNG = [figName, '.png'];
            filePath = fullfile(folderName, figNamePNG);
            set(gcf, 'PaperPositionMode', 'auto');
            exportgraphics(hfig, filePath, 'Resolution', 300);
        end

    

        function hintonDiagram(matrix, ax, figTitle, backgroundColor)
            % hintonDiagram - Generates a Hinton diagram to visualize matrix values.
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
            %                                color to a predefined (in Constants.m file) blue.
            %
            % Output:
            %   None. The function visualizes the matrix values as a Hinton diagram
            %   on the specified or current axis.
            if (nargin < 1)
                CustomError.raiseError('InputCheck', CustomError.ERR_NOT_ENOUGH_INPUT_ARG);
            end

            % Set default values for optional parameters
            if nargin < 2
                ax = gca;
            end
            if nargin < 3
                figTitle = '';
            end
            if nargin < 4
                backgroundColor = true;
            end
            
            % Set the current axis to ax
            axes(ax);
            axis(ax, 'equal'); % Set the same length in every direction
            
            if backgroundColor
                set(ax, 'Color', Visualization.hexToRGB(Constants.BLUE));
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
                title(ax, figTitle);
            end
        end

        
        
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
            %   - subfolderName (optional): The subfolder where the figure will be saved. 
            %                                If not provided, the figure will be saved in the current folder.
            %
            % Output:
            %   - None. The function visualizes the Hinton diagrams on the screen and optionally exports the figure.

            if (nargin < 1)
                CustomError.raiseError('InputCheck', CustomError.ERR_NOT_ENOUGH_INPUT_ARG);
            end

            saveFig = false;
            isSubfolderSpecified = false;
            if nargin > 3
                saveFig = true;
                if nargin > 4
                    isSubfolderSpecified = true;
                end
            end

            % Create figure
            hfig = figure;

            numPlots = length(arrW);

            hfigAx = axes(hfig); % Moved here for optimization purposes
            for i = 1:numPlots
                ax = Utility.ternary(numPlots == 1, hfigAx, subplot(1, numPlots, i));
                Visualization.hintonDiagram(arrW{i}, ax, ...
                    Utility.ternaryOpt(i <= length(subplotTitles), @() subplotTitles{i}, @() ''));
            end

            if numPlots > 1 && ~isempty(figTitle)
                sgtitle(figTitle);
            end
           
            % Format figure
            Visualization.formatFigure(hfig);

            % Save figure
            if saveFig
                Visualization.exportFigure(hfig, figName, ...
                    Utility.ternaryOpt(isSubfolderSpecified, @() subfolderName, @() ''));
            end
        end

        
        
        
        
        
        
        function plotLatentFactors(Z, figTitle, figName, subfolderName)
            % Optional parameters: figTitle, figName, subfolderName
            if nargin < 1
                error(['##### ERROR IN THE CLASS ' mfilename('class') ': Not enough input arguments provided.']);
            elseif nargin < 2
                figTitle = '';
            end

            hfig = figure;
            numFactors = size(Z, 1); % Z = [K x N]
            for i = 1:numFactors
                subplot(numFactors, 1, i);
                factor = Z(i, :); % Factors are in the rows of Z
                plot(factor, '.', 'MarkerSize', 4);
                hold on;
            end
            
            if ~isempty(figTitle)
                sgtitle(figTitle);
            end

            Visualization.formatFigure(hfig);

            % Save figure: if 'figTitle' is provided
            if nargin > 2
                Visualization.exportFigure(hfig, figName, ...
                    Utility.ternaryOpt(nargin == 2, @() '', @() subfolderName));
            end
        end
        
        function plotLoadings(W, dimList, figTitle, figName, subfolderName)
            % Parameters
            % ----------
            % W : matrix, [D_total x K]
            % dimList: number of features in each view 
            if nargin < 2
                error(['##### ERROR IN THE CLASS ' mfilename('class') ': Not enough input arguments provided.']);
            elseif nargin < 3
                figTitle = '';
            end

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
        
            hfig = figure;
            hold on;
        
            for k = 1:K
                y = W(:, k) + (K - k + 1) * offset;
        
                plot(x, y, 'Color', 'black', 'LineWidth', 1.5); %, 'Label', ['ind: ', k]);
            end
            
            ax = gca;
            ax.XAxis.Visible = 'off';
            ax.YAxis.Visible = 'off';
            
            linePos = 0;
            for k = 1:length(dimList) - 1 % Don't plot vertical line for the last view
                linePos = linePos + dimList(k);
                labelText = ['$D_{', num2str(k), '}$'];
                xline(linePos, 'Color', Constants.DARK_BLUE, 'LineWidth', 2, ...
                    'Label', labelText, 'LabelVerticalAlignment', 'top', ...
                    'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');
            end

            if ~isempty(figTitle)
                title(figTitle);
            end

            Visualization.formatFigure(hfig);

            % Save figure
            if nargin > 2
                Visualization.exportFigure(hfig, figName, ...
                    Utility.ternaryOpt(nargin == 3, @() '', @() subfolderName));
            end
        end
    
    

        % TODO: DRY this: check the function above
        function plotLoadingsAndAlpha(W, dimList, alpha, figsTitle, labelPos, figName, subfolderName)
            % Parameters
            % ----------
            % W : matrix, [D_total x K]
            % dimList: number of features in each view 
            if nargin < 2
                error(['##### ERROR IN THE CLASS ' mfilename('class') ': Not enough input arguments provided.']);
            % elseif nargin < 3
            %     figsTitle = '';
            end

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
        
            height = 400;
            hfig = figure('Position', [100, 100, 3.5 * height, height]); % TODO: the number * height must be a param, for 2G 2.5 was a good value
            subplot('Position', [0.1, 0.1, 0.6, 0.8]);
            hold on;
        
            for k = 1:K
                y = W(:, k) + (K - k + 1) * offset;
        
                plot(x, y, 'Color', 'black', 'LineWidth', 1.5); %, 'Label', ['ind: ', k]);
            end
            
            ax = gca;
            ax.XAxis.Visible = 'off';
            ax.YAxis.Visible = 'off';
            
            linePos = 0;
            for k = 1:length(dimList) - 1 % Don't plot vertical line for the last view
                linePos = linePos + dimList(k);
                labelText = ['$D_{', num2str(k), '}$'];
                xline(linePos, 'Color', Constants.DARK_BLUE, 'LineWidth', 2, ...
                    'Label', labelText, 'LabelVerticalAlignment', labelPos, ...
                    'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');
            end

            % if ~isempty(figTitle)
            %     title(figTitle);
            % end
            
            % title(figsTitle{1}, 'Interpreter', 'latex');

            subplot('Position', [0.75, 0.1, 0.2, 0.8]);
            Visualization.hintonDiagram(alpha, gca, '');

            % title(figsTitle{2}, 'Interpreter', 'latex');

            sgtitle(figsTitle, 'Interpreter', 'latex');

            Visualization.formatFigure(hfig);

            % % Save figure
            % if nargin > 3
            %     Visualization.exportFigure(hfig, figName, ...
            %         Utility.ternaryOpt(nargin == 3, @() '', @() subfolderName));
            % end
        end
        
    end
end
