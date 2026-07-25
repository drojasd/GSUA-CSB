# gsua-csb (Python)

Python port of the [GSUA-CSB MATLAB toolbox](https://github.com/drojasd/GSUA-CSB): sensitivity
analysis, uncertainty analysis, parameter estimation, and practical identifiability analysis
(including multiple-global-minima detection) for symbolic-ODE and user-defined models.

**Status: early/partial.** This is an in-progress port, not yet feature-complete with the MATLAB
toolbox. See "What's implemented" below for the current state.

## Why a separate Python package, same repo

The MATLAB source lives at the repository root (`Functions/`, `Examples/`, ...); this Python
package lives in `python/` alongside it. Both track the same GSUA-CSB project and release
cadence. Simulink-backed models are intentionally not ported — there is no Python equivalent, and a
Python user with a Simulink model should call MATLAB directly via the MATLAB Engine API rather than
have Simulink stepping reimplemented here.

## Installation (development)

```bash
cd python
pip install -e ".[dev]"
```

Optional extras:
- `.[symbolic]` — SymPy, for symbolic-ODE models (not yet implemented, see below)
- `.[sensitivity]` — SALib, for the standard Sobol/Saltelli/Jansen sensitivity estimators
- `.[all]` — everything above

## Design

The MATLAB toolbox's central data structure is a `table` with MATLAB-specific custom properties
(`Kind`, `Solver`, `Domain`, `Fixed`, ...) bolted on via `addprop`, dispatched on an integer `Kind`
inside nearly every function. That was a workaround for MATLAB not having a clean way to attach a
callable + metadata to tabular data. This package replaces it with a small `Model` class hierarchy
(`gsua_csb.Model`, `gsua_csb.UserFunctionModel`, ...) — real polymorphism instead of a switch
statement repeated in every function.

Each capability is exposed two ways:
- **Idiomatic name** (the real implementation): `gsua_csb.costf(...)`, `gsua_csb.median_ci(...)`
- **MATLAB-parity alias** (same function, for cross-reference with the MATLAB toolbox or its
  citations): `gsua_csb.gsua_costf(...)`, `gsua_csb.gsua_medianci(...)`

## What's implemented

| Capability | Status | MATLAB equivalent |
|---|---|---|
| `Model` / `UserFunctionModel` | Done | table + CustomProperties, `gsua_userdefined` |
| Cost functions (`costf`, `rcostf`, `costf_multi`, `likecost`) | Done | `gsua_costf`, `gsua_rcostf`, `gsua_costfMulti`, `gsua_likecost` |
| Coverage/tightness metric (`coverage_metric`) | Done | `gsua_covmetric` |
| Distribution-free median CI, band depth | Done | `gsua_medianCI`, `gsua_depth` |
| Symbolic-ODE models (SymPy + `scipy.integrate.solve_ivp`) | Done | `gsua_dpmat`, `gsua_odefun` |
| Sampling / design matrix (`scipy.stats.qmc`) | Done | `gsua_dmatrix` |
| Sensitivity analysis (SALib + custom Xiao method) | Planned | `gsua_sa` |
| Uncertainty analysis + Monte Carlo filtering | Planned | `gsua_ua`, `gsua_MCF` |
| Parameter estimation (`scipy.optimize`) | Planned | `gsua_pe` |
| Identifiability + multi-minima clustering (`sklearn`) | Planned | `gsua_ia`, `gsua_dia` |
| Range refinement, confidence sub-contour box | Planned | `gsua_oatr`, `gsua_oatr2`, `gsua_csb` |
| Profile likelihood | Planned | `gsua_likelihood` |
| Plotting (Matplotlib) | Planned | `gsua_plot` |
| Simulink models | **Not planned** — use the MATLAB Engine API instead | `sens_dataprep`, `sens_montecarlo` |

## Testing

```bash
pytest
```

## License

MIT — see [../LICENSE](../LICENSE) (same license as the MATLAB toolbox).
