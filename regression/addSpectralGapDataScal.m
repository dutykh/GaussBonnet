% Define helper function for adding data:

function dt = addSpectralGapDataScal(dt, alpha_value, L_value, measurements)

    newRow = table();
    newRow.alpha = alpha_value;
    newRow.L = L_value;
    newRow.deltaOmega = {measurements}; % Store as cell array
    dt = [dt; newRow];

end