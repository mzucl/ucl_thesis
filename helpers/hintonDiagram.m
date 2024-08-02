function hintonDiagram(W)
    % TODO: Add title as an optional parameter
    % TODO: Limit y and change the background color -> Bishop!
    maxWeight = max(abs(W(:)));
    figure;
    hold on;
    for i = 1:size(W, 1)
        for j = 1:size(W, 2)
            % Determine the size of the square
            weight = W(i, j);
            height = sqrt(abs(weight) / maxWeight);
            width = height;

            % White for positive, black for negative
            color = [1 1 1] * (weight >= 0);

            rectangle('Position', [j - width / 2, i - height / 2, width, height], 'FaceColor', color, 'EdgeColor', 'k');
        end
    end
    hold off;
    set(gca, 'YDir', 'reverse', 'XAxisLocation', 'top');
    title('Hinton Diagram of Principal Components');
end