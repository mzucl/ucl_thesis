classdef Visualization
    methods (Static)
        function hintonDiagram(matrix, ax, figureTitle)
            if nargin < 2
                % Use the current axis if none is provided
                ax = gca;
            end
            if nargin < 3
                figureTitle = '';
            end
            
            % Set the current axis to ax
            axes(ax);
            axis(ax, 'equal'); % same length in every direction

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
            axis(ax, 'off');
            if ~isempty(figureTitle)
                title(ax, figureTitle);
            end
        end

        function hintonDiagramPlot(arrW, titles, folderName, figName)
            if nargin < 4
                error(['##### ERROR IN THE CLASS ' mfilename('class') ': Not enough input arguments provided.']);
            end
            hfig = figure;

            ax1 = subplot(1, 2, 1);
            Visualization.hintonDiagram(arrW{1}, ax1, titles{1});
            
            ax2 = subplot(1, 2, 2);
            Visualization.hintonDiagram(arrW{2}, ax2, titles{2});
           
            picturewidth = 20;
            hw_ratio = 0.65;
            set(findall(hfig, '-property','FontSize'),'FontSize', 17);
            % set(findall(hfig, '-property', 'Box'), 'Box', 'off');
            set(findall(hfig, '-property', 'Interpreter'), 'Interpreter', 'latex');
            set(findall(hfig, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex');
            
            set(hfig, 'Units', 'centimeters', ...
                      'Position', [3 3 picturewidth hw_ratio * picturewidth]);
            
            pos = get(hfig, 'Position');
            
            set(hfig, 'PaperPositionMode', 'Auto', ...
                      'PaperUnits', 'centimeters', ...
                      'PaperSize', [pos(3), pos(4)]);

            % Save figure
            if ~exist(folderName, 'dir')
                mkdir(folderName);
            end

            filePath = fullfile(folderName, figName);

            print(hfig, filePath, '-dpdf', '-vector', '-fillpage');
            print(hfig, filePath, '-dpng', '-vector');

        end

        function plotStructVariables(resArr, offset)
            if nargin < 2
                offset = 1;
            end
            if offset < 1
                error(['##### ERROR IN THE CLASS Visualization' ': Offset is the starting iteration, it can not be ' ...
                    'less than 1']);
            end
            numIterations = length(resArr);
            
            % Check the first struct to get the field names
            firstRes = resArr{1};
            fieldNames = fieldnames(firstRes);
            
            numFields = length(fieldNames);
            
            % Determine the number of rows needed for subplots with 2 columns
            numRows = ceil(numFields / 2);
            
            figure;
            
            for i = 1:numFields
                subplot(numRows, 2, i);
                
                data = zeros(1, numIterations - (offset - 1));
                
                % Collect data across all iterations
                for j = offset:numIterations
                    data(j - (offset - 1)) = resArr{j}.(fieldNames{i});
                end
                
                % Plot the data
                plot(offset:numIterations, data, 'LineWidth', 1.5);
                
                % Set the title and labels
                title(['Variable: ', fieldNames{i}]);
                xlabel('Iteration');
                ylabel(fieldNames{i});
                grid on;
            end
            
            sgtitle('Evolution of ELBO variables over iteration');
        end

        function plotLoadings(W, dimList, figureTitle)
            % Parameters
            % ----------
            % W : matrix, [D_total x K]
            % dimList: number of features in each view 
            if nargin < 3
                figureTitle = '';
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
        
            figure;
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
                xline(linePos, 'Color', 'blue', 'LineWidth', 1.5, ...
                    'Label', labelText, 'LabelVerticalAlignment', 'top', ...
                    'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');
            end

            if ~isempty(figureTitle)
                title(figureTitle);
            end
        end
    end
end
