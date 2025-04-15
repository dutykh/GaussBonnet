%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: addSpectralGapData.m
%
% Authors: Dr Davide Batic & Dr Denys Dutykh
%          Khalifa University of Science and Technology, Abu Dhabi, UAE
%
% License: GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
%          (see https://www.gnu.org/licenses/gpl-3.0.html)
%
% Description:
%   Helper function to add a row of spectral gap data to a MATLAB table.
%   Used by the regression scripts for scalar, tensor, and vector perturbations.
%   The order of alpha and L does not matter, but the table must use the same
%   order for all rows.
%
% Usage:
%   dt = addSpectralGapData(dt, alpha_value, L_value, measurements)
%   where dt is a MATLAB table with 'alpha', 'L', and 'deltaOmega' columns.
%
% Contact:
%   For questions, contact the authors via their institutional pages.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
