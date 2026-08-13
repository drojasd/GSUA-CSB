# GSUA-CSB

[![MATLAB](https://img.shields.io/badge/MATLAB-R2019b%2B-blue.svg)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/205731654.svg)](https://zenodo.org/badge/latestdoi/205731654)

## 🐍 New: GSUA-CSB is now available for Python

Every non-Simulink capability of this toolbox — global sensitivity analysis (Sobol, Jansen,
Saltelli, and the toolbox's own Xiao method), uncertainty analysis and Monte Carlo filtering,
multistart parameter estimation, practical identifiability analysis (including spectral-clustering
detection of multiple global minima), range refinement, the namesake Confidence Sub-contour Box
algorithm, profile likelihood, and composable Matplotlib plotting — has a Python port in
[`python/`](python/), built on NumPy/SciPy/scikit-learn/SymPy. No MATLAB license required.

```bash
cd python && pip install -e ".[all]"
```

```python
import gsua_csb as gc
# gc.UserFunctionModel / gc.SymbolicODEModel, gc.design_matrix, gc.sensitivity_analysis,
# gc.parameter_estimation, gc.identifiability_analysis, gc.confidence_subcontour_box, ...
```

- **[`python/README.md`](python/README.md)** — what's implemented, install options, design notes
- **[`python/USERGUIDE.md`](python/USERGUIDE.md)** — worked examples for every capability
- **[`python/examples/system_identification_cycle.py`](python/examples/system_identification_cycle.py)**
  — the full semi-automated identification workflow, runnable end to end

Simulink-backed models are intentionally not ported (no Python equivalent — call MATLAB directly
via its Python Engine API instead); everything else targets full parity with the MATLAB toolbox
below, including symbolic-ODE models via SymPy.

**GSUA-CSB** (Global Sensitivity and Uncertainty Analysis — Confidence Sub-contour Box) is a MATLAB toolbox for
validating mathematical models implemented with Symbolic Math Toolbox or Simulink. It works uniformly across
three kinds of model definitions:

- **Simulink** models
- **Symbolic Math Toolbox** systems of ODEs
- **User-defined MATLAB functions** (black-box models)

The toolbox supports:

- variance-based global sensitivity analysis
- uncertainty analysis
- parameter estimation
- practical identifiability analysis (including detection of multiple global minima)
- confidence sub-contour box estimation for fitted parameters
- visualization workflows for model-validation studies

GSUA-CSB was developed at Universidad EAFIT and builds on previous work by Carlos Mario Vélez on GSUA for
dynamical systems using variance-based methods.

A single "GSUA table" (a MATLAB `table` carrying custom properties) records parameter ranges, nominal values,
and everything the toolbox needs to evaluate the underlying model, so the same functions
(`gsua_sa`, `gsua_ua`, `gsua_pe`, `gsua_ia`, ...) work regardless of which kind of model produced the table.

## Why It Matters

Mathematical models used in epidemiology, public health, engineering, and biological systems often depend on
uncertain parameters. GSUA-CSB helps researchers understand which parameters matter, how uncertainty propagates
through the model, and how identifiable fitted parameters are under available data.

## Links

- User guide and examples: [GSUA-CSB documentation](https://drojasd.github.io/GSUA-CSB/gsua_userguide)
- MATLAB File Exchange: [View GSUA-CSB on File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/72637-gsua-csb)
- Latest GitHub release: [releases](https://github.com/drojasd/GSUA-CSB/releases)
- DOI: [Zenodo DOI](https://zenodo.org/badge/latestdoi/205731654)

## Installation

### From File Exchange / Add-On Explorer
Search for "GSUA-CSB" in MATLAB's Add-On Explorer, or install directly from
[File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/72637-gsua-csb).
Every tagged GitHub release automatically syncs to File Exchange, so the Add-On Explorer always tracks the
latest release.

### From the packaged toolbox
Install `GSUA-CSB.mltbx` by double-clicking it in MATLAB.

### From Source
1. Clone this repository
2. Add the source folders to your MATLAB path:
   ```matlab
   addpath(genpath('Functions'))
   addpath(genpath('Add Funcs'))
   addpath('progressbar')
   ```
3. Or open `GSUA-UCI.prj` in MATLAB to load the project with its path already configured

## Getting Started

Open the getting-started guide for an interactive introduction:

```matlab
open doc/GettingStarted.mlx
```

or the full user guide:

```matlab
open gsua_userguide
```

## Typical workflow

```matlab
% 1. Describe the model (Simulink / symbolic ODE / user-defined function)
T = gsua_dataprep(modelName, Ranges, ParNames);

% 2. Sample the parameter space
M = gsua_dmatrix(T, 1000);

% 3. Analyze
[T,Y,J] = gsua_sa(M,T);           % global sensitivity analysis
Y       = gsua_ua(M,T);           % uncertainty analysis
[T,res] = gsua_pe(T,xdata,ydata); % parameter estimation

% 4. Check practical identifiability across repeated estimations,
%    including whether the estimations converge to more than one
%    distinct global minimum:
[T,clusterInfo] = gsua_ia(T,T.Estlsqc,false,false,true,true);
```

## Functions

### Model setup
| Function | Description |
|----------|-------------|
| `gsua_dataprep` | Consolidation of data and environment (dispatches to the right setup routine for Simulink / symbolic / user-defined models) |
| `gsua_dpmat` | Consolidation of data and symbolic MATLAB environment |
| `sens_dataprep` | Consolidation of data and Simulink environment |
| `gsua_userdefined` | Build a GSUA table around a user-defined MATLAB function |

### Sampling
| Function | Description |
|----------|-------------|
| `gsua_dmatrix` | Design of experiments (factor space sampling: Latin Hypercube, Uniform, Sobol) |

### Analysis
| Function | Description |
|----------|-------------|
| `gsua_sa` | Global sensitivity analysis (Xiao, Sobol, Jansen, Saltelli, brute-force, OAT methods) |
| `gsua_ua` | Uncertainty analysis via Monte-Carlo simulation, with automatic Monte-Carlo filtering |
| `gsua_pe` | Parameter estimation (lsqcurvefit, lsqnonlin, ga, particleswarm, patternsearch, surrogateopt, simulannealbnd, fmincon). The correlation-penalized multi-objective cost path (`'margin'` ≠ 0/1, i.e. `gsua_costf`) is NaN-tolerant per output row — a signal with its own gaps (`NaN` in `ydata`) no longer blinds other, fully-populated rows in the same joint fit. When that path is used, `T`'s `Margin`/`Alpha` `CustomProperties` are recorded automatically for later recovery by `gsua_noisefloor` |
| `gsua_likelihood` | Profile-likelihood confidence intervals for each parameter |

### Practical identifiability & confidence ranges
| Function | Description |
|----------|-------------|
| `gsua_ia` | Practical identifiability analysis with diagnostic plots — including optional spectral-clustering detection of **multiple global minima** among repeated estimation runs and optional fit-quality (`'cost'`) filtering of failed multistart runs before any statistic is computed |
| `gsua_dia` | Headless (no-plot) counterpart of `gsua_ia`, same clustering and fit-quality-filtering support |
| `gsua_costcutoff` | Shared fit-quality filter behind `gsua_ia`/`gsua_dia`'s `'cost'` option (fixed-tolerance or automatic-gap cutoff, with a minimum-runs floor) |
| `gsua_noisefloor` | Noise-calibrated fit-acceptance threshold via parametric bootstrap — replaces the scale-dependent `res < 1.5*res(1)` idiom with a cutoff on the `gsua_costf` cost scale derived from the data's own observation noise (Poisson / quasi-Poisson (default) / global-NB dispersion models, chosen for count-like or cumulative fitted outputs). Recovers `'margin'`/`'alpha'` from `gsua_pe`'s `T` automatically, or accepts explicit overrides |
| `gsua_medianCI` | Distribution-free confidence interval for the median |
| `gsua_oatr` / `gsua_oatr2` | Once-at-a-time range expansion/reduction |
| `gsua_csb` | Uncertainty-based confidence sub-contour box estimation |

### Evaluation & visualization
| Function | Description |
|----------|-------------|
| `gsua_eval` | Evaluate the model for a handful of parameter sets, with optional plotting |
| `gsua_deval` | Evaluate the model behind a GSUA table for one parameter set (internal dispatcher) |
| `gsua_pardeval` | Evaluate the model behind a GSUA table for a batch of parameter sets |
| `gsua_plot` | Visualization function for sensitivity/uncertainty/identifiability results |
| `gsua_MCF` | Monte-Carlo filtering plots (behavioral vs. non-behavioral parameter sets) |
| `gsua_mtest` | Animate model output as one parameter is swept, saved as a GIF |
| `gsua_covmetric` | Percentile-band (P5/P50/P95) coverage-accuracy and band-tightness metrics for `gsua_ua` results, on the same margin-normalized `gsua_costf` scale (below 1 = within tolerance) — the stopping signal for semi-automated identification |

### Persistence
| Function | Description |
|----------|-------------|
| `gsua_save` | Save a GSUA table to disk in a portable form |
| `gsua_load` | Rebuild a GSUA table's Solver handle after loading from disk |

### Internal helpers
`gsua_costf`, `gsua_rcostf`, `gsua_costfMulti`, `gsua_likecost` (objective/cost functions for optimizers),
`gsua_corr2omitnan` (NaN-tolerant `corr2`, used by `gsua_costf`),
`gsua_intrp` (ODE solution interpolation), `gsua_timer` (progress reporting), `gsua_depth` (band depth ranking),
`gsua_odefun` (symbolic-to-function-handle ODE conversion), `gsua_ref` (range-adequacy refinement),
`sens_montecarlo` (Simulink Monte-Carlo evaluation engine).

## Examples

See the `Examples/` folder:
- `default_example.m` — minimal end-to-end run
- `gsua_main.m` — full workflow walkthrough
- `user_free.m` / `user_dependent.m` — user-defined-function model examples

## Requirements

MATLAB with Statistics and Machine Learning Toolbox, Optimization Toolbox, and Global Optimization Toolbox
recommended for full functionality. Symbolic Math Toolbox and/or Simulink are required depending on which
kind of model you work with.

## Methods Implemented

Sensitivity index estimators implemented in this toolbox are based on:

1. Saltelli, A., Annoni, P., Azzini, I., Campolongo, F., Ratto, M., and Tarantola, S. (2010). Variance based
   sensitivity analysis of model output: design and estimator for the total sensitivity index. *Computer
   Physics Communications*, 181(2), 259-270.
2. Xiao, S., Lu, Z., and Wang, P. (2018). Multivariate global sensitivity analysis based on distance
   components decomposition. *Risk Analysis*, 38(12), 2703-2721.

## Citation

If you use this toolbox, cite:

Rojas-Diaz, Daniel and Velez-Sanchez, Carlos Mario (2019). GSUA-CSB (https://www.github.com/drojasd/GSUA-CSB),
GitHub. doi:10.5755281/zenodo.3383316.

## License

This project is distributed under the MIT License — see [LICENSE](LICENSE). Third-party code bundled under
`Add Funcs/` retains its own license (see the `LICENSE.md`/`README.md` inside that folder).

---

## For Contributors

Source lives under `Functions/` (public API), `Add Funcs/` (bundled third-party utilities), `Examples/`, and
`doc/`. Static analysis is run with MATLAB's Code Analyzer (`checkcode`); there is currently no automated test
suite — contributions adding `matlab.unittest` coverage under a `tests/` folder are welcome.

**Note on releases:** this repository is linked to MATLAB File Exchange — publishing a GitHub Release
automatically updates the File Exchange listing. Only tag a release when the toolbox is genuinely ready to
publish.
