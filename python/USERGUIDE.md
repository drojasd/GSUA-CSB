# gsua-csb Python user guide

This is a practical, example-driven walkthrough of every capability in the Python port. For a
short summary of what's implemented and how it maps to the MATLAB toolbox, see
[README.md](README.md); for the full narrated semi-automation workflow, see
[`examples/system_identification_cycle.py`](examples/system_identification_cycle.py).

Every code block below is runnable as shown (they're drawn from the package's own test suite, so
they're known to work against the current release).

## Installation

```bash
cd python
pip install -e ".[dev]"        # core + symbolic-ODE support requires the extra below
pip install -e ".[all]"        # core + SymPy (symbolic-ODE models) + dev tools
```

```python
import gsua_csb as gc
```

## 1. Defining a model

Everything in the toolbox operates on a `Model`: parameter names, bounds, a nominal point, and a
way to evaluate the model at a given parameter vector. There are two concrete `Model` classes.

### `UserFunctionModel` — wrap any Python callable

```python
import numpy as np
from gsua_csb import UserFunctionModel

def exp_decay(params, xdata):
    amplitude, rate = params
    return amplitude * np.exp(-rate * xdata)

model = UserFunctionModel(
    func=exp_decay,
    names=["amplitude", "rate"],
    range=np.array([[0.5, 4.0], [0.05, 1.5]]),   # [lower, upper] per parameter
    nominal=np.array([2.0, 0.5]),                 # defaults to the midpoint of range if omitted
    domain=np.linspace(0, 5, 25),                 # the x-axis func is evaluated over
)

model.evaluate(model.nominal, model.domain)   # -> array of shape (25,)
```

A parameter is **fixed** when its lower and upper bound are equal (`range[i, 0] == range[i, 1]`).
Every function below automatically excludes fixed parameters from sampling/estimation/analysis.
`model.fix("rate")` collapses a parameter to its current nominal value in place — the standard move
when practical identifiability analysis says a parameter can't be estimated (see §6).

For a *scalar* (domain-less) model, just omit `domain` and write `func(params)`.

### `SymbolicODEModel` — a system of ODEs defined with SymPy

Requires the `symbolic` extra (`pip install -e ".[symbolic]"`).

```python
import sympy as sp
from gsua_csb import SymbolicODEModel

t = sp.Symbol("t")
x, y = sp.symbols("x y")             # state variables
a, b, c, d = sp.symbols("a b c d")   # true parameters

model = SymbolicODEModel(
    odes=[a * x - b * x * y, -c * y + d * x * y],   # Lotka-Volterra
    state_vars=[x, y],
    t=t,
    params=[a, b, c, d],
    domain=np.linspace(0, 10, 50),
    names=["x0", "y0", "a", "b", "c", "d"],           # initial conditions + true params
    range=np.tile([[0.5, 1.5]], (6, 1)),
    solver_kwargs={"rtol": 1e-8, "atol": 1e-10},      # passed to scipy.integrate.solve_ivp
)
```

The parameter vector is `[initial_conditions..., true_params...]`, matching the `names`/`range`
order above. `model.evaluate(params, xdata)` returns an `(n_states, len(xdata))` array.

## 2. Sampling a design matrix

```python
from gsua_csb import design_matrix

M = design_matrix(model, n=500, method="latin_hypercube", seed=0)  # (500, n_params)
```

`method` is `"latin_hypercube"` (default, stratified), `"sobol"` (low-discrepancy), or `"uniform"`.
Fixed parameters are held constant in every row rather than sampled.

## 3. Global sensitivity analysis

```python
from gsua_csb import sensitivity_analysis

result = sensitivity_analysis(model, M, method="xiao")  # "xiao" (default), "sobol", "jansen",
                                                          # "saltelli", or "oat"
result.Si       # (Np,) first-order index on the scalar SSE cost vs. the nominal run
result.STi      # (Np,) total-order index
result.Si_vec   # (Np, Nd) time-dependent first-order index (None for method="oat")
result.STi_vec  # (Np, Nd) time-dependent total-order index
```

`"xiao"` is the toolbox's own robust estimator [Xiao, Lu & Wang, 2018] and the only one of the
five with no equivalent in `SALib`. Pass `y_exp=` to score against real data instead of the
model's own nominal run, and `output_index=` to pick which state to analyze for a multi-state
model (e.g. a `SymbolicODEModel`).

## 4. Uncertainty analysis and Monte Carlo filtering

```python
from gsua_csb import uncertainty_analysis, monte_carlo_filter

ua = uncertainty_analysis(model, M, y_exp=ydata)   # ua.Y: (N, Nd) ensemble; ua.y_nom: reference
mcf = monte_carlo_filter(model, M, ua.Y, ua.y_nom)  # behavioral vs. non-behavioral split
mcf.names   # free parameter names
mcf.low     # (N_low, Np_free) parameter sets whose cost fell below the reference
mcf.high    # (N_high, Np_free) parameter sets whose cost fell above the reference
```

## 5. Parameter estimation (multistart)

```python
from gsua_csb import parameter_estimation

pe = parameter_estimation(
    model, xdata, ydata,
    n=8,                       # 8 independent starting points
    solver="least_squares",    # "least_squares", "minimize", "differential_evolution", "dual_annealing"
    seed=0,
)
pe.x[0]      # best-fit parameter vector (rows sorted by fit quality, x[0] is the best)
pe.cost[0]   # its cost
```

`n > 1` is not optional overhead — it's what makes the next step (identifiability analysis)
meaningful. A single estimate can't tell you whether the optimizer found *the* answer or *an*
answer.

## 6. Practical identifiability analysis

```python
from gsua_csb import identifiability_analysis

ia = identifiability_analysis(model, pe.x, cost=pe.cost, correction=True)
ia.range              # (Np, 2) new confidence interval per parameter
ia.index              # (Np,) identifiability index in [0, 1] -- 0 well identified, 1 poorly identified
ia.correlation        # (Np, Np) correlation matrix among the repeated estimates
ia.n_bad_fit_removed  # runs dropped for converging to a much worse cost than the best one
```

Pass `cost=pe.cost` whenever the repeated estimates came from `parameter_estimation` (or any other
scored multistart run). A run that converged to a bad local optimum contributes essentially
arbitrary parameter values — without this filter, those failed runs can make a genuinely
well-identified parameter *look* poorly identified. `identifiability_analysis` drops any run whose
cost exceeds `best_cost * (1 + cost_rtol)` (default 10% tolerance) before computing anything else.

A high `index` or a strong pairwise `correlation` (e.g. `|r| > 0.9`) usually means two parameters
are structurally entangled — the data constrains their *combination*, not either one individually.
The standard fix is `model.fix(name)` on whichever of the pair has the wider interval, then
re-running estimation on what's left. See
[`examples/system_identification_cycle.py`](examples/system_identification_cycle.py) for this
played out end to end on a model where only a product of two parameters is identifiable.

### Detecting multiple global minima

```python
ia = identifiability_analysis(model, pe.x, cluster=True, max_k=5, sil_threshold=0.5)
ia.cluster.num_clusters   # 1 if the repeated estimates form a single basin
ia.cluster.centers        # (num_clusters, Np) candidate global minima
```

Repeated estimations from a non-convex problem can land in distinct basins of attraction instead
of scattering around one point. `cluster=True` runs spectral clustering across candidate cluster
counts and keeps the split only if its mean silhouette score clears `sil_threshold` — a low score
means the "best" split found isn't actually well separated, so the estimates are treated as one
basin.

## 7. Range refinement and the Confidence Sub-contour Box

```python
from gsua_csb import range_refinement, confidence_subcontour_box

rr = range_refinement(model, lim=0.3)   # one-at-a-time expansion/reduction per free parameter
rr.range

csb = confidence_subcontour_box(model, n=300, y_exp=ydata, reps=40, lim=0.3, stop=0.5, seed=0)
csb.range               # final narrowed box
csb.converged           # whether `stop`'s behavioral fraction was reached within `reps`
csb.behavioral_fraction # (n_iterations,) fraction of runs within tolerance at each iteration
```

`confidence_subcontour_box` is the toolbox's namesake algorithm: it iteratively narrows the *whole*
parameter box toward the region consistent with `ydata`, rather than one parameter at a time.

## 8. Profile-likelihood confidence intervals

```python
from gsua_csb import profile_likelihood

pl = profile_likelihood(model, xdata, ydata=ydata, alpha=0.95, margin=0.1)
pl.range   # (Np, 2) confidence interval from the likelihood-ratio test
```

An alternative to the sampling-based intervals above, grounded directly in the likelihood-ratio
test: each parameter is stepped away from its nominal value (refitting every other free parameter
at each step) until the fit degrades past a chi-squared threshold.

## 9. Plotting

Every `plot_*` function takes a result object directly and an optional `ax`, so you can compose
them into your own figure layout (standard Matplotlib/pandas convention):

```python
import matplotlib.pyplot as plt
from gsua_csb import plot_uncertainty, plot_sensitivity_bar, plot_identifiability_correlation

fig, axes = plt.subplots(1, 3, figsize=(15, 4))
plot_uncertainty(ua, ax=axes[0])
plot_sensitivity_bar(result, index="Si", ax=axes[1])
plot_identifiability_correlation(ia, ax=axes[2])
fig.tight_layout()
```

Also available: `plot_sensitivity_area` (time-dependent indices, stacked), `plot_identifiability_index`
(sorted index bar chart), and `plot_mcf` (prior/low/high ECDF panels, one per free parameter).

## 10. The semi-automation identification cycle

The pattern above — sample, estimate, check identifiability, act on what you find, repeat — is the
toolbox's core workflow. [`examples/system_identification_cycle.py`](examples/system_identification_cycle.py)
runs it end to end:

1. **Reachability check** — before spending any optimizer budget, confirm the real data actually
   falls inside the model's plausible output range.
2. **Multistart estimation** — fit repeatedly, not once.
3. **Identifiability analysis** — turn repeated estimates into confidence intervals, a
   per-parameter index, and a correlation matrix.
4. **Confidence-interval thinness check** (`coverage_metric`) — two independent questions: does the
   identified region fit the data (`cost_data`), and is it tight (`cost_band`)? A model can fit
   well while remaining unidentifiable, and vice versa.
5. **Corrective action** — fix whichever of an entangled pair is worse-constrained, and repeat.

Run it directly to see the numbers: `python examples/system_identification_cycle.py`.

## MATLAB → Python name reference

Every function below is also available under its idiomatic Python name (e.g. `costf` for
`gsua_costf`) — see [README.md](README.md#design) for the naming convention.

| MATLAB | Python |
|---|---|
| `gsua_costf`, `gsua_rcostf`, `gsua_costfMulti`, `gsua_likecost` | `costf`, `rcostf`, `costf_multi`, `likecost` |
| `gsua_covmetric` | `coverage_metric` |
| `gsua_medianCI`, `gsua_depth` | `median_ci`, `band_depth` |
| `gsua_dmatrix` | `design_matrix` |
| `gsua_sa` | `sensitivity_analysis` |
| `gsua_ua`, `gsua_MCF` | `uncertainty_analysis`, `monte_carlo_filter` |
| `gsua_pe` | `parameter_estimation` |
| `gsua_ia`, `gsua_dia` | `identifiability_analysis` |
| `gsua_oatr`, `gsua_oatr2` | `range_refinement` |
| `gsua_csb` | `confidence_subcontour_box` |
| `gsua_likelihood` | `profile_likelihood` |
| `gsua_plot` | `plot_uncertainty`, `plot_sensitivity_bar`, `plot_sensitivity_area`, `plot_identifiability_correlation`, `plot_identifiability_index`, `plot_mcf` |
| `gsua_dpmat`, `gsua_odefun` | `SymbolicODEModel` |
| `gsua_userdefined` | `UserFunctionModel` |

## Getting help

Found a bug, or a place where this port's behavior diverges from the MATLAB original in a way
that isn't documented in the relevant module's docstring? Please open an issue at
[github.com/drojasd/GSUA-CSB](https://github.com/drojasd/GSUA-CSB/issues).
