# Gauss-Bonnet and Schwarzschild-Tanghelini's QNMs spectral analysis project

## Authors
Dr Davide Batic & Dr Denys Dutykh  
Khalifa University of Science and Technology, Abu Dhabi, UAE

## Analytical Computations (Maple)

This repository also contains three Maple worksheet files with the analytical computations supporting our study:

- `maple/01-NONEXTREME CASE-SUPPLEMENTAL MATERIAL-SCALAR CASE.mw`  
  Analytical calculations for the scalar perturbation case.
- `maple/02-NONEXTREME CASE-SUPPLEMENTAL MATERIAL-TENSORIAL CASE.mw`  
  Analytical calculations for the tensorial perturbation case.
- `maple/03-NONEXTREME CASE-SUPPLEMENTAL MATERIAL-VECTORIAL CASE.mw`  
  Analytical calculations for the vectorial perturbation case.

These files contain symbolic derivations and supplemental material relevant to the regression analysis performed in the MATLAB scripts.

## Linear Regression Analysis

The `regression` folder contains MATLAB scripts that perform linear regression analysis on numerically computed spectral gap data for quasinormal modes (QNMs) of black holes in Gauss-Bonnet and Schwarzschild-Tangherlini spacetimes. Each script corresponds to a different type of perturbation (scalar, tensor, or vector) and analyzes the dependence of the spectral gap on the coupling parameter $\alpha$ and the mode number $L$.

### Main MATLAB Scripts
- `scalarL.m` — Regression for scalar perturbations
- `tensorL.m` — Regression for tensor perturbations
- `vectorVI.m`, `vectorX.m`, `vectorXV.m` — Regression for various vector perturbations

### What the Scripts Do
- Read and organize numerical data for different $(\alpha, L)$ pairs
- Compute the mean and standard deviation of the spectral gap for each pair
- Perform linear regression to fit the spectral gap as a function of $L$
- Visualize the results, saving figures (PDF/PNG) in the `shots` directory
- Save processed data in the `data` directory as `.mat` files

### How to Run
1. Open MATLAB or GNU Octave in the `regression` directory.
2. Run any of the scripts (e.g., `scalarL.m`) by typing its name:
   ```matlab
   scalarL
   ```
   The script will automatically create `shots` and `data` subfolders if they do not exist.

### Connection to Analytical Results
The regression analysis is informed by the analytical computations in the Maple worksheet files. The symbolic derivations in the Maple files provide theoretical predictions and context for the numerical trends observed and quantified in the regression scripts.

For more details, see the comments in each script and the analytical Maple files.

## License
This project is licensed under the GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007.  
See [LICENSE](LICENSE) for details.

---
