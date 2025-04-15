% Script file named "scalarL.m"

close all

% Initialize folders:
if ~exist('shots', 'dir')
    mkdir('shots');
end

if ~exist('data', 'dir')
    mkdir('data');
end

% Initialize the table:
dataTable = table('Size', [0, 3], 'VariableTypes', {'double', 'double', 'cell'}, ...
    'VariableNames', {'alpha', 'L', 'deltaOmega'});

% Add data only for α = 0.2 and α = 0.5
% α = 0.2, L = 0
dataTable = addSpectralGapData(dataTable, 0.2, 0, [0.7271, 0.7274, 0.7272, 0.7273, 0.7273]);
% α = 0.2, L = 1
dataTable = addSpectralGapData(dataTable, 0.2, 1, [0.7275, 0.7277, 0.7277, 0.7276, 0.7276]);
% α = 0.2, L = 2
dataTable = addSpectralGapData(dataTable, 0.2, 2, [0.7282, 0.7282, 0.7281, 0.7281, 0.7281]);
% α = 0.5, L = 0
dataTable = addSpectralGapData(dataTable, 0.5, 0, [0.7264, 0.7263, 0.7264, 0.7265, 0.7264]);
% α = 0.5, L = 1
dataTable = addSpectralGapData(dataTable, 0.5, 1, [0.7272, 0.7271, 0.7271, 0.7271, 0.7271]);
% α = 0.5, L = 2
dataTable = addSpectralGapData(dataTable, 0.5, 2, [0.7286, 0.7284, 0.7282, 0.7282, 0.7281, 0.7284]);

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

% Linear regression analysis for the two alpha values
uniqueAlpha = unique(dataTable.alpha);
% Use colors that work well in both color and B&W
colors = {[0.8500, 0.3250, 0.0980],  % Strong orange for α = 0.2
          [0.0, 0.4470, 0.7410]};    % Deep blue for α = 0.5
% Line styles
linestyles = {'-', '--'};
% Markers for individual points
dataMarkers = {'o', 's'};
% Special markers for average points
avgMarkers = {'*', 'p'};

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
        scatter(repmat(L_j, length(deltaOmega_j), 1), deltaOmega_j, 25,...
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
    fprintf('\\hat{\\alpha} = %.2f: \\Delta = %.6f·L + %.6f, R² = %.6f\n', alpha_value, b(2), b(1), rsq);

    % Plot regression line
    L_range = linspace(0, 2.2, 100);  % Extend slightly beyond the data range
    delta_fit = b(1) + b(2) * L_range;
    lineObjects(i) = plot(L_range, delta_fit, 'Color', colors{i}, 'LineStyle', linestyles{i}, 'LineWidth', 2);
end

% Configure plot
xlabel('$\ell_3$', 'FontSize', 18, 'Interpreter', 'latex');
ylabel('$\Delta\Omega$', 'FontSize', 16, 'Interpreter', 'latex');

% Create simplified legend
legend(lineObjects, {'$\ \hat{\alpha} = 0.2$', '$\ \hat{\alpha} = 0.5$'}, ...
    'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 14, 'Box', 'off');

% Remove grid
grid off;
box on;
set(gca, 'FontSize', 14);

% Adjust axis limits for better visualization
xlim([-0.1, 2.2]);
% Keep y-axis focused on the data range with a small margin
ymin = min(dataTable.meanDeltaOmega) - 0.0005;
ymax = max(dataTable.meanDeltaOmega) + 0.0005;
ylim([ymin, ymax]);

% Set tight layout to save only the figure area
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 6];
fig.PaperSize = [8 6];

% Save the figure in both PNG and PDF formats with tight layout
saveas(gcf, 'shots/scalar-linreg-01.png');
print('shots/scalar-linreg-01', '-dpdf', '-r300', '-bestfit');

% Save the data to a .mat file
save('data/scalar-linreg-01.mat', 'dataTable');

fprintf('Analysis complete. Figures saved in the "shots" folder.\n');
---------- Forwarded message ---------
From: Lu Lu via SIAM <Mail@connectedcommunity.org>
Date: Mon, Apr 14, 2025 at 7:43 PM
Subject: SIAG on Computational Science and Engineering Community : Join the Yale/UNC-CH Geophysical Waveform Inversion Kaggle Competition!
To: <ahmaddeeb25@gmail.com>


￼
SIAG on Computational Science and Engineering Community
Post New Message
Join the Yale/UNC-CH Geophysical Waveform Inversion Kaggle Competition!
Reply to Group	Reply to Sender
￼	
Apr 14, 2025 11:43 AM  |  ￼  view attached
Lu Lu
Dear all, 

We are excited to announce the Yale/UNC-CH Geophysical Waveform Inversion Competition on Kaggle! This competition invites participants to develop machine learning algorithms for estimating subsurface properties, such as velocity maps, from seismic waveform data. 

By entering this challenge, you'll have the opportunity to bridge the gaps in seismic analysis by combining physics and machine learning. Your participation could lead to significant advancements in full waveform inversion, which has applications in subsurface energy exploration, medical diagnostics, non-destructive material testing, and more. 

Entry Deadline: June 23, 2025 
Prizes: 
1st Place - $12,000 
2nd Place - $10,000 
3rd Place - $10,000 
4th Place - $10,000 
5th Place - $8,000 
This is a unique opportunity to showcase your skills, contribute to a critical area of research, and win substantial prizes. We look forward to seeing your innovative solutions and wish you the best of luck in the competition. 

Best regards, 

Lu Lu, Assistant Professor, Yale University
Youzuo Lin, Associate Professor, University of North Carolina at Chapel Hill 



------------------------------
Lu Lu
Assistant Professor of Statistics and Data Science
Faculty of Institute for Foundations of Data Science
Faculty of Wu Tsai Institute
Faculty of Center for Algorithms, Data, and Market Design at Yale
Faculty of Center for Materials Innovation
Yale University
https://lugroup.yale.edu

------------------------------
  Reply to Group Online   View Thread   Recommend   Forward   Flag as Inappropriate  



 
You are subscribed to "SIAG on Computational Science and Engineering Community" as ahmaddeeb25@gmail.com. To change your subscriptions, go to My Subscriptions. To unsubscribe from this community discussion, go to Unsubscribe.
￼
