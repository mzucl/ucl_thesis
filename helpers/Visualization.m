classdef Visualization
    properties(Constant)
        FIGURES_FOLDER = 'figures';
    end

    methods (Static)
        function formatFigure(hfig)
            set(findall(hfig, '-property','FontSize'),'FontSize', 17);
            set(findall(hfig, '-property', 'Box'), 'Box', 'off');
            set(findall(hfig, '-property', 'Interpreter'), 'Interpreter', 'latex');
            set(findall(hfig, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex');
        end

        function rgb = hexToRGB(hex)
            rgb = sscanf(hex(2:end),'%2x%2x%2x',[1 3]) / 255;
        end

        function exportFigure(hfig, figName, subfolderName)
            % EXPORTFIGURE Save or export a MATLAB figure to a file.
            % hfig      : handle to the figure
            % filename  : name of the file (e.g., 'myfigure')

            % Optional parameters: subfolderName
            if nargin < 2
                error(['##### ERROR IN THE CLASS ' mfilename('class') ': Not enough input arguments provided.']);
            elseif nargin > 3
                error(['##### ERROR IN THE CLASS ' mfilename('class') ': Too many input arguments provided.']);
            end
            
            width = 20;
            ratio = 0.65; % height/weight ratio

            set(hfig, 'Units', 'centimeters', ...
                      'Position', [3 3 width ratio * width]);
            
            pos = get(hfig, 'Position');
            
            set(hfig, 'PaperPositionMode', 'Auto', ...
                      'PaperUnits', 'centimeters', ...
                      'PaperSize', [pos(3), pos(4)]);

            folderName = Visualization.FIGURES_FOLDER;
            if nargin > 2 && ~isempty(subfolderName)
                folderName = [Visualization.FIGURES_FOLDER, '/', subfolderName];
            end

            % Save figure
            if ~exist(folderName, 'dir')
                mkdir(folderName);
            end

            % % Export to .pdf
            % figNamePDF = [figName, '.pdf'];
            % filePath = fullfile(folderName, figNamePDF);
            % set(gcf, 'PaperPositionMode', 'auto');
            % exportgraphics(hfig, filePath, 'ContentType', 'vector');

            % Export to .png
            figNamePNG = [figName, '.png'];
            filePath = fullfile(folderName, figNamePNG);
            set(gcf, 'PaperPositionMode', 'auto');
            exportgraphics(hfig, filePath, 'Resolution', 300);
        end     
    
        function hintonDiagram(matrix, ax, figTitle, backgroundColor)
            % Optional parameters: ax, figTitle, backgroundColor
            if nargin < 2
                % Use the current axis if none is provided
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
            axis(ax, 'equal'); % same length in every direction
            
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

        function plotHintonDiagrams(arrW, figName, subfolderName)
            % Optional parameters: figName, subfolderName
            if nargin < 1
                error(['##### ERROR IN THE CLASS ' mfilename('class') ': Not enough input arguments provided.']);
            end

            hfig = figure;

            numPlots = length(arrW);

            for i = 1:numPlots
                ax = subplot(1, numPlots, i);
                Visualization.hintonDiagram(arrW{i}, ax);
            end

            Visualization.formatFigure(hfig);

            % Save figure
            if nargin > 1
                Visualization.exportFigure(hfig, figName, ...
                    Utility.ternaryOpt(nargin == 2, @() '', @() subfolderName));
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
