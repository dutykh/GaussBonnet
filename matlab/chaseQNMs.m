close all;
clear all;
format longE;

% --- User Configuration ---
% !!! IMPORTANT: Set the path to your Advanpix Multiprecision Computing Toolbox !!!
advanpixPath = '/path/to/advanpix/'; 
addpath(advanpixPath);

% Define the resolutions to use for the calculation
list = [100, 120, 140]; % Example resolutions, adjust as needed
L = length(list); 

% Define the directory containing the data files (M0, M1, M2 matrices)
dataDir = 'data/'; 
% Define the output file name for computed QNMs
outputFile = 'qnms.mat';
% --- End User Configuration ---


% --- QNM Computation ---
qnms = cell(1, L); % Cell array to store QNMs for each resolution

fprintf('Starting QNM computation for resolutions: [%s]\n', num2str(list));

for idx = 1:L
  n = list(idx);
  nstr = num2str(n);

  fprintf('  Processing resolution n = %d ... ', n);

  % Set precision for Advanpix
  mp.Digits(n);

  % Construct file paths for the matrices
  m0name = fullfile(dataDir, ['M0_', nstr, '.mat']);
  m1name = fullfile(dataDir, ['M1_', nstr, '.mat']);
  m2name = fullfile(dataDir, ['M2_', nstr, '.mat']);

  % Check if files exist
  if ~exist(m0name, 'file') || ~exist(m1name, 'file') || ~exist(m2name, 'file')
      fprintf('error: Matrix files not found for n = %d. Skipping.\n', n);
      qnms{idx} = []; % Store empty to indicate failure
      continue; % Skip to the next resolution
  end

  % Load matrices using Advanpix
  try
      M0 = mp.read(m0name);
      M1 = mp.read(m1name);
      M2 = mp.read(m2name);
  catch ME
      fprintf('error loading matrix files for n = %d: %s. Skipping.\n', n, ME.message);
      qnms{idx} = []; % Store empty to indicate failure
      continue; % Skip to the next resolution
  end

  % Solve the polynomial eigenvalue problem: (M2*lambda^2 + 1i*M1*lambda + M0) * x = 0
  try
      e = polyeig(M0, mp('1i')*M1, M2);
      qnms{idx} = e; % Store the computed eigenvalues (QNMs)
      fprintf('done.\n');
  catch ME
      fprintf('error during polyeig for n = %d: %s. Skipping.\n', n, ME.message);
      qnms{idx} = []; % Store empty to indicate failure
      continue; % Skip to the next resolution
  end

end % for idx

fprintf('QNM computation finished.\n');


% --- Save Results ---
fprintf('Saving computed QNMs to %s ... ', outputFile);
try
    save(outputFile, 'qnms', 'list'); % Save QNMs and the resolutions used
    fprintf('done.\n');
catch ME
    fprintf('error saving results: %s\n', ME.message);
end


% --- Basic Plotting ---
fprintf('Generating plot of computed QNMs...\n');
figure();
set(gcf, 'Color', 'w'); % White background
hold on; 
grid on;

colors = lines(L); % Get distinct colors for each resolution

plotted_indices = []; % Keep track of which resolutions were successfully plotted
legends_list = {}; % Store legend entries

for idx = 1:L
    if ~isempty(qnms{idx}) % Check if QNMs were computed and stored for this resolution
        plot(real(qnms{idx}), imag(qnms{idx}), 'o', ...
             'MarkerFaceColor', colors(idx,:), ...
             'MarkerEdgeColor', 'k', ... 
             'MarkerSize', 6, ...
             'LineStyle', 'none');
        plotted_indices(end+1) = idx; % Record that this index was plotted
        legends_list{end+1} = ['n = ' num2str(list(idx))]; % Create legend entry
    else
        fprintf('  Skipping plot for resolution n = %d (no data or computation error).\n', list(idx));
    end
end % for idx

if ~isempty(plotted_indices)
    legend(legends_list, 'Location', 'best'); % Add legend only if something was plotted
else
    fprintf('  No data successfully computed or plotted.\n');
end

xlabel('Re(\omega)');
ylabel('Im(\omega)');
title('Computed Quasinormal Modes (QNMs)');
hold off;

fprintf('Plotting finished. Figure displayed.\n');

% --- End of Script ---