xi = 2.5;
x = linspace(-6, 12, 1000);

bounds = {BohningBound(xi), JaakkolaBound(xi), LaplaceApproximation(xi)};
titles = {'$B\ddot{o}hning$', '$Jaakkola$', '$Laplace$'};

% Define sigma function
sigma = 1 ./ (1 + exp(-x));

hfig = figure;
plot(x, sigma, 'LineWidth', 2, 'DisplayName', '$\sigma(a)$'); 
hold on;

xlim([-6 12]);
% xlabel('$a$', 'Interpreter', 'latex');

% title('Quadratic bounds and approximation of $\sigma(a)$ at $\xi = 2.5$', Interpreter='latex');

% Define bound functions
for i = 1:length(bounds)
    bound = bounds{i};
    boundFun = bound.c() + bound.g() * (x - xi) + 1/2 * bound.h() * (x - xi).^2;
    plot(x, exp(x - boundFun), 'LineWidth', 1, 'DisplayName', titles{i});
    legend('Interpreter', 'latex');
    hold on;
end

verticalLinesPos = {-xi, xi};
vecticalLinesLabels = {'$-\xi$', '$\xi$'};

for i = 1:length(verticalLinesPos)
    linePos = verticalLinesPos{i};
    xline(linePos, 'Color', Constants.DARK_BLUE, 'LineWidth', 1, ...
        'Label', vecticalLinesLabels{i}, 'LabelVerticalAlignment', 'top', ...
        'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex', 'HandleVisibility', 'off');
end

legendHandle = legend('show');
set(legendHandle, 'Location', 'eastoutside');
grid on;

Visualization.formatFigure(hfig);
Visualization.exportFigure(hfig, 'bounds', '');
