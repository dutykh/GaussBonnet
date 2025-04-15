% Unified helper function for adding spectral gap data
% Usage: dt = addSpectralGapData(dt, alpha_value, L_value, measurements)
% The order of alpha and L does not matter, but the table must use the same order for all rows.

function dt = addSpectralGapData(dt, alpha_value, L_value, measurements)
    % Detect table variable order and assign appropriately
    varNames = dt.Properties.VariableNames;
    newRow = table();
    if ismember('alpha', varNames) && ismember('L', varNames)
        newRow.alpha = alpha_value;
        newRow.L = L_value;
    elseif ismember('L', varNames) && ismember('alpha', varNames)
        newRow.L = L_value;
        newRow.alpha = alpha_value;
    else
        error('Table must have variables named "alpha" and "L"');
    end
    newRow.deltaOmega = {measurements}; % Store as cell array
    dt = [dt; newRow];
end
