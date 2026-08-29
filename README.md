# Quasinormal Modes of Gauss–Bonnet Black Holes via the Spectral Method

Supplementary computational material (Maple, MATLAB/Octave) for the study of
scalar, vector and tensor quasinormal modes (QNMs) of higher-dimensional
Schwarzschild–Tangherlini black holes corrected by a Gauss–Bonnet (GB) term.

## Authors

**Dr Davide Batic** ([ORCID 0000-0002-7619-4900](https://orcid.org/0000-0002-7619-4900))
and **Dr Denys Dutykh** ([ORCID 0000-0001-5247-2788](https://orcid.org/0000-0001-5247-2788))

Mathematics Department, Khalifa University of Science and Technology,
Abu Dhabi, United Arab Emirates.

## Reference

> D. Batic and D. Dutykh,
> *Quasinormal modes of Gauss–Bonnet black holes via the spectral method:
> Scalar, vector, and tensor perturbations*,
> **Phys. Rev. D 114**, 044015 (2026).

- Journal: <https://doi.org/10.1103/91q6-r3jd>
- Preprint: [arXiv:2608.06083 [gr-qc]](https://arxiv.org/abs/2608.06083)
- Code repository: <https://github.com/dutykh/GaussBonnet>
- Catalogue entry: [Quasi-Normal Modes Hall of Fame](https://qnms.denys-dutykh.com/)

BibTeX:

```bibtex
@article{Batic2026,
  author        = {Batic, Davide and Dutykh, Denys},
  title         = {Quasinormal modes of {G}auss--{B}onnet black holes via the
                   spectral method: {S}calar, vector, and tensor perturbations},
  journal       = {Phys. Rev. D},
  volume        = {114},
  number        = {4},
  pages         = {044015},
  year          = {2026},
  doi           = {10.1103/91q6-r3jd},
  eprint        = {2608.06083},
  archivePrefix = {arXiv},
  primaryClass  = {gr-qc}
}
```

## Scope of the study

The paper presents a unified treatment of the three perturbation sectors (spin
$s = 0, 1, 2$) of Einstein–Gauss–Bonnet black holes in $D \in \{5, 6, 7, 8, 10,
11, 12, 26\}$, using a high-precision Chebyshev spectral method that remains
reliable well beyond the range of sixth-order WKB and characteristic-integration
techniques. The main findings are the appearance of overdamped purely imaginary
modes, non-monotonic real parts of the higher overtones, a strong amplification
of the dimensionless QNM frequencies in string-motivated dimensions, an exact
scalar-monopole/vector-dipole isospectrality at vanishing GB coupling (proved
analytically through a Darboux factorisation), and the first numerical
confirmation of the six-dimensional tensor instability, with the mass threshold
$G M \leqslant 158.1\,\alpha^{3/2}$.

## Computational pipeline

The material in this repository implements the numerical part of the study in
four stages:

1. **Symbolic derivation** (`maple/`). The perturbation equation of each sector
   is reduced, in a compactified radial variable, to a second-order ODE whose
   coefficients are computed in closed form, together with the effective
   potentials $U_S$, $U_V$, $U_T$ and the auxiliary functions $K$, $q$.
2. **Matrix assembly** (`matlab/*.mpl`). The ODE is discretised on a Chebyshev
   collocation grid, producing the three matrices $M_0$, $M_1$, $M_2$ of a
   quadratic eigenvalue problem, exported in MATLAB ASCII format.
3. **Eigenvalue solution** (`matlab/chaseQNMs.m`). The quadratic pencil is
   solved in multiprecision arithmetic, and the QNM frequencies are extracted
   and compared across spectral resolutions.
4. **Post-processing** (`regression/`). The spectral gap of the purely imaginary
   branch is regressed linearly against the angular momentum number, for several
   values of the GB coupling.

### Conventions

Throughout the Maple sources, `n` denotes the dimension of the horizon sphere
$S^n$, so that the spacetime dimension is $D = n + 2$: the mass term of the
metric function carries the Schwarzschild–Tangherlini exponent
$n - 1 = D - 3$, and the angular eigenvalues are $\lambda = \ell(\ell + n - 1)$
in the scalar and tensor sectors and $\lambda = (\ell - 1)(\ell + n)$ in the
vector sector. The dimensionless coupling is written `ah` ($\hat{\alpha}$) and
the horizon radius `rhoh` ($\rho_h$). The radial coordinate is compactified onto
$y \in [-1, 1]$, with $y = -1$ at the event horizon, where the metric function
vanishes, and $y = +1$ at spatial infinity, so that the exterior region is
exactly the interval on which the Chebyshev basis is defined.

## Repository layout

```
GaussBonnet/
├── maple/                             Symbolic derivations (Maple 2024)
│   ├── 01-...-SCALAR CASE.mw
│   ├── 02-...-TENSORIAL CASE.mw
│   └── 03-...-VECTORIAL CASE.mw
├── matlab/                            Matrix assembly and QNM computation
│   ├── matrixassembler-scalar.mpl     Collocation matrices, spin s = 0
│   ├── matrixassembler-vector.mpl     Collocation matrices, spin s = 1
│   ├── matrixassembler-tensor.mpl     Collocation matrices, spin s = 2
│   └── chaseQNMs.m                    Quadratic eigenvalue solver
├── regression/                        Spectral-gap regression analysis
│   ├── addSpectralGapData.m           Table-building helper
│   ├── scalarL.m, tensorL.m
│   ├── vectorVI.m, vectorX.m, vectorXV.m
│   ├── data/                          Saved tables (.mat)
│   └── shots/                         Generated figures (PDF/PNG)
├── LICENSE                            GNU GPL, version 3
└── README.md
```

## Analytical computations (Maple)

Three Maple 2024 worksheets carry the symbolic derivations that support the
numerical work, one per perturbation sector, all for the non-extreme case:

| Worksheet | Sector | Angular eigenvalue |
|---|---|---|
| `maple/01-NONEXTREME CASE-SUPPLEMENTAL MATERIAL-SCALAR CASE.mw` | scalar, $s = 0$ | $\lambda = \ell(\ell + n - 1)$ |
| `maple/02-NONEXTREME CASE-SUPPLEMENTAL MATERIAL-TENSORIAL CASE.mw` | tensor, $s = 2$ | $\lambda = \ell(\ell + n - 1)$ |
| `maple/03-NONEXTREME CASE-SUPPLEMENTAL MATERIAL-VECTORIAL CASE.mw` | vector, $s = 1$ | $\lambda = (\ell - 1)(\ell + n)$ |

Each worksheet locates the event horizon as the relevant root of
$\rho_h^{\,n-1} + \hat{\alpha}\,\rho_h^{\,n-3} - 2 = 0$, forms the coupling
$\kappa = 8\hat{\alpha}/\rho_h^{\,n+1}$, builds the metric function $f$ together
with the derivatives needed by the sector at hand, and assembles the effective
potential and the coefficients of the second-order ODE. The worksheets open in
Maple 2024 or later; the parameters at the top (`n`, `alphah`, `L`) are set to a
representative case and are meant to be edited.

## Matrix assembly (Maple)

The three `matlab/matrixassembler-*.mpl` files each define a single procedure:

```maple
MatrixAssembler(d, N, L, n, ah, p)
```

| Argument | Meaning |
|---|---|
| `d`  | number of digits used in the symbolic/numeric evaluation |
| `N`  | number of Chebyshev modes (spectral resolution) |
| `L`  | angular momentum number ($\ell \geqslant 0, 1, 2$ for $s = 0, 1, 2$) |
| `n`  | dimension of the horizon sphere, $D = n + 2$ |
| `ah` | dimensionless Gauss–Bonnet coupling $\hat{\alpha}$ |
| `p`  | path whose `data/` subfolder receives the exported matrices |

The procedure expands the unknown in the Chebyshev basis,
$P(y) = \sum_{j=0}^{N-1} a_j T_j (y)$, collocates the ODE at the Chebyshev roots
$y_i = \cos\!\big((2i-1)\pi/(2N)\big)$, and extracts the coefficient of each
$a_j$ to build

$$\big(M_2\,\omega^2 + \mathrm{i}\,M_1\,\omega + M_0\big)\,\mathbf{a} = 0 .$$

The scalar and tensor versions extract the coefficients term by term; the vector
version precomputes the Chebyshev value and derivative matrices and caches the
coefficient functions with `option remember`, which makes it appreciably faster
at large `N`. All three write `M0_<N>.mat`, `M1_<N>.mat` and `M2_<N>.mat` into
`<p>/data/` through `ExportMatrix(..., target = MATLAB, mode = ascii)`.

Typical session:

```maple
read "matrixassembler-scalar.mpl":
MatrixAssembler(150, 100, 0, 3, 0.1, "/path/to/GaussBonnet/matlab"):
```

## QNM computation (MATLAB)

`matlab/chaseQNMs.m` reads the assembled matrices and solves the quadratic
eigenvalue problem in multiprecision arithmetic.

What the script does:

- loops over a list of spectral resolutions (`list = [100, 120, 140]` by
  default), setting the working precision with `mp.Digits(n)` for each;
- loads `data/M0_<n>.mat`, `data/M1_<n>.mat`, `data/M2_<n>.mat` with `mp.read`,
  skipping any resolution whose files are missing;
- solves $(M_2\lambda^2 + \mathrm{i}M_1\lambda + M_0)\mathbf{x} = 0$ with
  `polyeig`, the eigenvalues $\lambda$ being the QNM frequencies;
- stores the spectra of all resolutions in a cell array saved to `qnms.mat`,
  which allows the physical modes (those stable under refinement) to be told
  apart from the spurious ones;
- plots the computed frequencies in the complex plane,
  $\mathrm{Re}(\omega)$ against $\mathrm{Im}(\omega)$, colour-coded by
  resolution.

Requirements:

- MATLAB (R2016a or later is recommended for `polyeig`);
- the [Advanpix Multiprecision Computing Toolbox](https://www.advanpix.com/),
  which is required: the pencil is severely ill-conditioned at high resolution
  and double precision is not sufficient.

How to run:

1. Edit the `advanpixPath` variable at the top of `chaseQNMs.m` so that it
   points at your Advanpix installation.
2. Generate the matrices for the resolutions of interest with the Maple
   assembler of the relevant sector, so that `matlab/data/` contains the
   matching `M0_*.mat`, `M1_*.mat` and `M2_*.mat` files.
3. Adjust `list` to the resolutions you generated.
4. Start MATLAB in the `matlab` directory and run

   ```matlab
   chaseQNMs
   ```

The script reports its progress in the command window, writes `matlab/qnms.mat`
and displays the QNM spectrum.

## Linear regression analysis (MATLAB/Octave)

The `regression` folder studies how the spectral gap $\Delta\Omega$ of the
purely imaginary branch varies with the angular momentum number $\ell_n$, for
several values of the coupling $\hat{\alpha}$. Each script embeds the frequencies
obtained at successive spectral resolutions, averages them, performs a linear
least-squares fit of $\Delta\Omega$ against $\ell_n$, reports the slope, the
intercept and the coefficient of determination $R^2$, and produces the
corresponding figure.

| Script | Sector | Horizon sphere | Couplings $\hat{\alpha}$ | Output |
|---|---|---|---|---|
| `scalarL.m`  | scalar | $n = 3$ ($D = 5$) | 0.2, 0.5      | `scalar-linreg-01` |
| `tensorL.m`  | tensor | $n = 4$ ($D = 6$) | 0.1, 0.2, 0.5 | `tensor-linreg-01` |
| `vectorVI.m` | vector | $n = 3$ ($D = 5$) | 0.2, 0.5, 1.0 | `vector-linreg-01` |
| `vectorX.m`  | vector | $n = 4$ ($D = 6$) | 0.1, 0.2, 0.5 | `vector-linreg-02` |
| `vectorXV.m` | vector | $n = 5$ ($D = 7$) | 0.1, 0.2      | `vector-linreg-03` |

The shared helper `addSpectralGapData.m` appends one row
$(\hat{\alpha}, \ell, \{\Delta\Omega\})$ to the data table, keeping the repeated
measurements in a cell so that both the mean and the standard deviation remain
available.

Each script writes its figure to `regression/shots/` in PNG and PDF form and the
underlying table to `regression/data/` as a `.mat` file; both subfolders are
created if absent. The versions committed here reproduce the figures of the
paper.

How to run:

```matlab
cd regression
scalarL     % or tensorL, vectorVI, vectorX, vectorXV
```

These scripts use only core language features and run under MATLAB or GNU
Octave; no toolbox is needed.

## Reproducing a spectrum end to end

1. Choose a sector, a dimension $D = n + 2$ and a coupling $\hat{\alpha}$.
2. Run the matching `matrixassembler-*.mpl` in Maple for a few increasing values
   of `N`, with the working precision `d` at least equal to `N`, which is the
   precision at which `chaseQNMs.m` then solves the eigenvalue problem.
3. Run `chaseQNMs.m` on that list of resolutions.
4. Retain the eigenvalues that are stable across resolutions; these are the
   physical QNMs. The remaining ones are discretisation artefacts.
5. For the purely imaginary branch, feed the resulting spectral gaps into the
   relevant regression script.

## Quasi-Normal Modes Hall of Fame

<https://qnms.denys-dutykh.com/>

The *Quasi-Normal Modes Hall of Fame* ("Celestial Music") is the public
catalogue of the QNM research group at Khalifa University. It gathers the
group's peer-reviewed manuscripts on quasinormal modes of black holes,
wormholes, quantum-corrected and non-commutative geometry inspired compact
objects, and lists for each of them the full citation, the DOI, the journal and
arXiv links, and the companion code repository. The present article, together
with the codes collected here, is part of that catalogue.

## License

Released under the GNU General Public License, version 3, 29 June 2007.
See [LICENSE](LICENSE) for the full text.
