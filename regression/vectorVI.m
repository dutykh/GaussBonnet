% Script file named "vectorVI.m"

close all

% Initialize folders
if ~exist('shots', 'dir')
    mkdir('shots');
end

if ~exist('data', 'dir')
    mkdir('data');
end

% Initialize the table
dataTable = table('Size', [0, 3], 'VariableTypes', {'double', 'double', 'cell'}, ...
    'VariableNames', {'alpha', 'L', 'deltaOmega'});

% α = 0.2
dataTable = addSpectralGapDataVec(dataTable, 0.2, 1, [0.7271, 0.7270, 0.7270, 0.7271, 0.7270]);
dataTable = addSpectralGapDataVec(dataTable, 0.2, 2, [0.7275, 0.7274, 0.7274, 0.7274, 0.7274]);
dataTable = addSpectralGapDataVec(dataTable, 0.2, 3, [0.7279, 0.7279, 0.7278, 0.7279, 0.7278]);

% α = 0.5
dataTable = addSpectralGapDataVec(dataTable, 0.5, 1, [0.7260, 0.7261, 0.7260, 0.7261, 0.7262]);
dataTable = addSpectralGapDataVec(dataTable, 0.5, 2, [0.7270, 0.7269, 0.7269, 0.7269, 0.7269]);
dataTable = addSpectralGapDataVec(dataTable, 0.5, 3, [0.7281, 0.7280, 0.7280, 0.7278, 0.7278]);

% α = 1.0
dataTable = addSpectralGapDataVec(dataTable, 1.0, 1, [0.7163, 0.7164, 0.7164, 0.7165, 0.7165]);
dataTable = addSpectralGapDataVec(dataTable, 1.0, 2, [0.7199, 0.7194, 0.7191, 0.7189, 0.7186]);
dataTable = addSpectralGapDataVec(dataTable, 1.0, 3, [0.7248, 0.7237, 0.7229, 0.7221, 0.7216]);

% α = 0.2, L = 4
dataTable = addSpectralGapDataVec(dataTable, 0.2, 4, [0.7286, 0.7284, 0.7285, 0.7284, 0.7283]);
% α = 0.5, L = 4
dataTable = addSpectralGapDataVec(dataTable, 0.5, 4, [0.7299, 0.7297, 0.7295, 0.7293, 0.7291]);
% α = 1.0, L = 4 (using the lowest group as requested)
dataTable = addSpectralGapDataVec(dataTable, 1.0, 4, [0.7268, 0.7287, 0.7327, 0.7370]);

% Calculate means for each combination of alpha and L
numRows = height(dataTable);
meanDeltaOmega = zeros(numRows, 1);
stdDeltaOmega = zeros(numRows, 1);

for i = 1:numRows
    measurements = dataTable.deltaOmega{i};
    meanDeltaOmega(i) = mean(measurements);
    stdDeltaOmega(i) = std(measurements);
end

% Add these statistics to the table
dataTable.meanDeltaOmega = meanDeltaOmega;
dataTable.stdDeltaOmega = stdDeltaOmega;

% Display the table with means
disp(dataTable);

% Linear regression analysis for each alpha value using means
uniqueAlpha = unique(dataTable.alpha);
% Colors that work well in both color and B&W - using darker green instead of light green
colors = {[0.8500, 0.3250, 0.0980],  % Strong orange
          [0.0, 0.4470, 0.7410],     % Deep blue
          [0.1500, 0.4000, 0.0800]}; % Dark green instead of light green
% Line styles
linestyles = {'-', '--', ':'};
% Markers for individual points
dataMarkers = {'o', 's', 'd'};
% Special markers for average points
avgMarkers = {'*', 'p', 'h'};

% Create figure with white background
figure('Position', [100, 100, 800, 600], 'Color', 'w');
hold on;

% Store line objects for legend
lineObjects = gobjects(length(uniqueAlpha), 1);

% Process each alpha value
for i = 1:length(uniqueAlpha)
    alpha_value = uniqueAlpha(i);
    disp(['Analysis for \hat{\alpha} = ', num2str(alpha_value)]);

    % Extract data points for this alpha
    indices = dataTable.alpha == alpha_value;
    alpha_data = dataTable(indices, :);

    % Plot individual data points (with solid color but smaller)
    for j = 1:height(alpha_data)
        L_j = alpha_data.L(j);
        deltaOmega_j = alpha_data.deltaOmega{j};
        scatter(repmat(L_j, length(deltaOmega_j), 1), deltaOmega_j, 25, ...
            colors{i}, dataMarkers{i}, 'filled', 'MarkerFaceAlpha', 0.5);
    end

    % Extract L values and mean deltaOmega values for regression
    L_values = alpha_data.L;
    mean_deltaOmega = alpha_data.meanDeltaOmega;

    % Plot average points (solid/larger)
    scatter(L_values, mean_deltaOmega, 100, colors{i}, avgMarkers{i}, 'filled', ...
        'LineWidth', 1.5);

    % Linear regression using only the averages
    X = [ones(length(L_values), 1), L_values];  % Add column of ones for intercept
    b = X \ mean_deltaOmega;  % Solve linear system

    % Calculate fitted values
    y_fit = X * b;

    % Calculate R²
    SSresid = sum((mean_deltaOmega - y_fit).^2);
    SStotal = sum((mean_deltaOmega - mean(mean_deltaOmega)).^2);
    rsq = 1 - SSresid/SStotal;

    % Display regression coefficients and R²
    fprintf('\\hat{\\alpha} = %.1f: \\Delta = %.6f·L + %.6f, R² = %.6f\n', alpha_value, b(2), b(1), rsq);

    % Plot regression line
    L_range = linspace(min(L_values)-0.2, max(L_values)+0.2, 100);
    delta_fit = b(1) + b(2) * L_range;
    lineObjects(i) = plot(L_range, delta_fit, 'Color', colors{i}, 'LineStyle', linestyles{i}, 'LineWidth', 2);
end

% Configure plot
xlabel('$\ell_3$', 'FontSize', 18, 'Interpreter', 'latex');
ylabel('$\Delta\Omega$', 'FontSize', 16, 'Interpreter', 'latex'); % Removed "(imaginary part)"

% Create simplified legend
legend(lineObjects, {'$\ \hat{\alpha} = 0.2$', '$\ \hat{\alpha} = 0.5$', '$\ \hat{\alpha} = 1.0$'}, ...
    'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 14, 'Box', 'off');

% Remove grid
grid off;
box on;
set(gca, 'FontSize', 14);

% Adjust axis limits for better visualization
xlim([0.8, 4.2]);
ylim_current = ylim;
% Add a small margin around the data
ylim([ylim_current(1)-0.0005, ylim_current(2)+0.0005]);

% Set tight layout to save only the figure area
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 6];
fig.PaperSize = [8 6];

% Save the figure in both PNG and PDF formats with tight layout
saveas(gcf, 'shots/vector-linreg-01.png');
print('shots/vector-linreg-01', '-dpdf', '-r300', '-bestfit');

% Save the data to a .mat file in the 'data' subfolder
save('data/vector-linreg-01.mat', 'dataTable');

fprintf('Analysis complete. Figures saved in the "shots" folder and data saved in the "data" folder.\n');