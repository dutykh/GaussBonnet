% Define helper function for adding data:

function dt = addSpectralGapDataTens(dt, L_value, alpha_value, measurements)

    newRow = table();
    newRow.L = L_value;
    newRow.alpha = alpha_value;
    newRow.deltaOmega = {measurements}; % Store as cell array
    dt = [dt; newRow];

end