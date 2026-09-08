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

### Choosing which output to work with

A model that returns several signals — the three compartments of an SIR system, say — usually has
only one of them measured. `model.output` selects which, and is the counterpart of
`T.Properties.CustomProperties.output` on a MATLAB summary table:

```python
model.set_output("infected")     # by name
model.set_output(1)              # or by index
model.set_output([0, 2])         # or several
model.set_output(None)           # back to all outputs
model.active_output_names        # -> what evaluate() currently returns
```

It can also be given at construction: `SymbolicODEModel(..., output=1)`.

Set it once and everything follows, because every routine in the package reaches the model
through `model.evaluate`. That includes `parameter_estimation`, `profile_likelihood`,
`confidence_subcontour_box` and `noise_floor`, none of which take an `output_index` argument — so
before this existed, fitting one state of a multi-state model meant wrapping it in a second model
just to slice the output. The `output_index` arguments on `sensitivity_analysis` and
`uncertainty_analysis` still work and index the *active* outputs, i.e. whatever is left after this
selection.

Prefer `model.evaluate(params, xdata)` over calling your own function directly: it is the
toolbox's evaluation path, it applies this selection, and it is what every routine uses
internally. It is the counterpart of MATLAB's `gsua_eval`.

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

### Importing a model from PEtab

Requires the `petab` extra (`pip install -e ".[petab]"`, or `.[all]`). [PEtab](https://petab.readthedocs.io/)
is a community format (SBML + a few TSV tables) for specifying "fit this model to this data"
problems, used across the systems-biology tooling ecosystem (AMICI/pyPESTO, COPASI, Data2Dynamics),
with a curated collection of real published models and data:
[Benchmark-Models-PEtab](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab).

```python
from gsua_csb import load_petab

problem = load_petab("Boehm_JProteomeRes2014/Boehm_JProteomeRes2014.yaml")
problem.model             # a Model, ready for sensitivity_analysis/parameter_estimation/etc.
problem.xdata              # measurement time points
problem.ydata              # (n_observables, len(xdata)) measurement values
problem.observable_names   # e.g. ['pSTAT5A_rel', 'pSTAT5B_rel', 'rSTAT5A_rel']

sim = problem.model.evaluate(problem.model.nominal, problem.xdata)  # (n_observables, len(xdata))
```

The returned `model` handles everything transparently: species initial conditions and SBML
parameters become `names`/`range`/`nominal` (free if PEtab marks them `estimate=1`, fixed
otherwise), and `evaluate()` internally simulates the ODE system and applies each observable's
formula (which can be a derived expression of multiple species, not just a bare species name) —
so it plugs directly into `parameter_estimation`, `sensitivity_analysis`, `identifiability_analysis`,
and everything else in this guide with no special-casing.

This is a scoped importer, not a general PEtab implementation — one problem, one SBML model. A
problem with more than one simulation condition (e.g. the same model fit jointly across several
regions or intervention scenarios, common in epidemiology PEtab files) returns a
`{conditionId: PEtabProblem}` dict instead of a bare `PEtabProblem`; a condition with no
measurements at all (a forward-projection "what if" scenario alongside the one condition actually
fit to data) is skipped with a warning rather than raising. Condition-table overrides are
supported: a literal numeric cell fixes that target for the condition, and a cell naming another
parameter table row makes *that* parameter — free or fixed, with its own bounds — the thing
estimated for that target in that condition (e.g. a region-specific transmission rate). See the
`parse_sbml`/`load_petab` docstrings for the exact SBML feature subset supported (compartments,
assignment rules including chained ones, function definitions, `piecewise`/`&&`/`||`, and initial
assignments — including one deriving a constant parameter from others — are all handled; rate
rules, events, and preequilibration are not).

Verified against five real benchmarks: `Boehm_JProteomeRes2014` (every supported systems-biology
SBML feature at once) and `Giordano_Nature2020` (an epidemiology model exercising function
definitions, time-dependent piecewise assignment rules, and an assignment-rule-defined "reporter"
species used directly as an observable) both match their own published reference simulation to
~1e-4 relative error or better; `Bertozzi_PNAS2020` (two regions, each with condition-table
parameter-reference overrides) matches to ~1e-7. `Okuonghae_ChaosSolitonsFractals2020` (literal
condition-table overrides) has no published reference simulation to check against but produces
finite, plausibly-scaled output. `Perelson_Science1996`'s own reference simulation does not match a
simulation at its `nominalValue` parameters (~45% error, wrong qualitative shape) even though the
import mechanism is the same one validated to high precision on the other four — most likely that
reference file corresponds to different (e.g. fitted) parameter values, not an importer bug, but
this is unconfirmed; see the test suite for the full investigation.

## 2. Sampling a design matrix

```python
from gsua_csb import design_matrix

M = design_matrix(model, n=500, method="latin_hypercube", seed=0)  # (500, n_params)
```

`method` is `"latin_hypercube"` (default, stratified), `"sobol"` (low-discrepancy), or `"uniform"`.
Fixed parameters are held constant in every row rather than sampled.

A free parameter spanning many orders of magnitude (common for rate constants in epidemiology/
systems-biology models — e.g. one PEtab benchmark's bound is `[1e-13, 1000]` with a true value of
`2e-12`) needs `Model.log_scale` set for that position: linear-uniform sampling over such a bound
puts effectively all its mass in the top decade, so the true value is essentially never sampled
anywhere close to. `load_petab` sets this automatically from PEtab's `parameterScale` column; for a
model built by hand, pass `log_scale=[...]` to `UserFunctionModel`/`SymbolicODEModel`. Both
`design_matrix` and `parameter_estimation` respect it (search happens in log10 space internally;
every value you pass in or get back — `initial_points`, `pe.x` — stays in natural/linear units).

### Joint (correlation-preserving) sampling

The three methods above are **marginal**: each parameter is drawn independently inside its own
interval. For a model with parameter confounding that discards the joint structure, so draws leave
the identified manifold and any band built from them reports uncertainty the inference does not
actually contain. `method="joint"` instead draws whole parameter vectors from an accepted-estimate
ensemble:

```python
ia = identifiability_analysis(model, pe.x, cost=pe.cost)
M  = design_matrix(model, 2000, "joint", seed=0, pool=ia.estimates_used)
```

`IdentifiabilityResult.estimates_used` is the pool that produced `ia.range`, after fit-quality
filtering, dominant-cluster restriction and outlier removal. Passing it also sidesteps a second
problem: `ia.range` is a confidence interval of the **median** of that pool, which narrows as the
pool grows and is not the spread of values consistent with the data — joint draws carry the pool's
own spread instead. For that reason `clip` defaults to no clipping; pass the original model's
physical bounds if you need clipping, not `ia.range`.

`joint_type` selects `"bootstrap"` (default — exact correlations, cannot leave the manifold, but
limited to the vectors already in the pool), `"smooth_bootstrap"` (adds a variance-corrected kernel
for new points near the manifold) or `"gaussian"` (mean + covariance; assumes the manifold is
linear, so on a curved ridge it places mass beside it). A pool smaller than `min_pool_n` (default 5,
matching `min_corr_n`) is refused with a warning and falls back to marginal sampling — correlation
from two points is exactly ±1 regardless of any real relationship, so such a band would be
confidently wrong rather than approximately right.

Note that a parameter-uncertainty band is not a prediction interval: observation noise is a separate
term, modelled by [`noise_floor`](#9-noise-calibrated-fit-acceptance-threshold). See
`paper_experiments/joint_vs_marginal_bands.py` for a worked comparison on real data.

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
well-identified parameter *look* poorly identified. By default (`cost_method="rtol"`),
`identifiability_analysis` drops any run whose cost exceeds `best_cost * (1 + cost_rtol)` (default
10%) before computing anything else. A fixed tolerance is a somewhat arbitrary answer to a
genuinely data-dependent question, so an alternative is available:

```python
ia = identifiability_analysis(model, pe.x, cost=pe.cost, cost_method="gap")
```

`cost_method="gap"` instead finds the largest jump in the *sorted* costs (in log space) and cuts
there — the automated version of reading a "waterfall plot" by eye: a cluster of converged runs,
a jump, then a scattered tail of failures. It only cuts if that jump is at least `cost_gap_ratio`
(default 3x) larger than the typical gap between sorted costs, so a smooth continuum of similar
costs is left untouched rather than sliced arbitrarily. Either method respects `min_keep` (default
3): whatever the cutoff decides, at least `min_keep` of the best-cost runs always survive, so a
too-aggressive filter can't silently collapse the analysis down to a single degenerate point.

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

When a genuine multi-cluster split is found, `ia.range`/`ia.correlation`/`ia.index`/`ia.nominal`
are computed from the **dominant (largest) cluster's points only** — pooling separated basins into
one confidence interval or correlation isn't a meaningful summary. The other basins are still
fully available via `ia.cluster.centers`. If `outlier=True` is combined with `cluster=True`,
outlier removal runs *after* clustering and is scoped to that dominant cluster only: running
Mahalanobis-distance outlier detection *before* clustering would treat a genuine second basin as
"outliers" relative to the pooled mean/covariance and risk deleting it before clustering ever gets
a chance to find it.

A genuine multi-cluster split can leave very few points in the dominant cluster — two well-
separated basins with two runs each is a perfectly legitimate outcome, not a hypothetical. Any 2
points correlate at exactly ±1 regardless of any real relationship, so `ia.correlation`/`ia.index`
would otherwise report every parameter as maximally, spuriously "strongly correlated" with every
other. Below `min_corr_n` points (default 5), `ia.correlation` is left `NaN`,
`ia.correlation_reliable` is `False`, and `ia.index` falls back to interval width alone rather than
averaging in a number that looks meaningful but isn't:

```python
ia = identifiability_analysis(model, pe.x, cluster=True, min_corr_n=5)
ia.correlation_reliable   # False if the dominant cluster is too small to trust a correlation from
```

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

## 9. Noise-calibrated fit-acceptance threshold

```python
from gsua_csb import noise_floor

pe = parameter_estimation(model, xdata, ydata, n=2000, solver="minimize", margin=0.1)
out = noise_floor(model, xdata, ydata, pe.x, pe.cost, margin=pe.margin, alpha=pe.alpha)
out.threshold    # calibrated cutoff on the same cost scale as pe.cost
out.accepted     # (N,) bool -- pe.cost < out.threshold

accepted_x = pe.x[out.accepted]
ia = identifiability_analysis(model, accepted_x, cost=pe.cost[out.accepted])
```

Replaces the common but scale-dependent idiom `lims = sum(pe.cost < 1.5*pe.cost.min())`: that rule
asks "how close is this to the best fit found", not "is this statistically distinguishable from a
perfect model given how noisy the data actually is" — it gets *narrower* as the fit improves
(backwards), and has no mechanism to guarantee the accepted band actually covers the data it was
fit to. `noise_floor` instead simulates the best pool member, estimates an observation-noise model
from its own residuals, bootstraps synthetic datasets under that noise, and scores the *true*
model against each synthetic dataset with the same cost function `pe.cost` is on — the resulting
distribution's upper quantile (default 95%) is the largest cost still consistent with a correct
model, calibrated against the data's own noise rather than against how well the optimizer happened
to do.

`margin`/`alpha` must be passed through from whatever produced the pool (typically `pe.margin`/
`pe.alpha` from a `PEResult`) — they are not auto-recovered, since a threshold computed with the
wrong ones would be silently meaningless. Three selectable noise models (`noise=`): `"poisson"`
(almost always too tight for real data, never the default), `"quasipoisson"` (recommended default,
`var = phi*mu`), and `"nbglobal"` (`var = mu + mu**2/k`, report as a bound only — a single global
`k` is dominated by the largest-`mu` points). `out.by_model` reports all three regardless, as a
sensitivity check. For a fitted output that is a *cumulative* series, pass `cumulative=[...]` (one
flag per output row) — noise is estimated/generated on the incident (first-differenced) series and
re-accumulated, not applied to the cumulative series directly.

This is a goodness-of-fit band-sizing calibration, not a formal confidence region — it doesn't
replace `profile_likelihood`/`confidence_subcontour_box` above, it's complementary and needs no
refitting. It's a *parametric* bootstrap conditional on the best fit being approximately correct;
under model misspecification the floor is optimistic. One nuance worth knowing: the calibrated
*threshold* is intentionally not invariant to a pure rescale of the data's units (it's tied to a
physically meaningful count-noise model), unlike `costf`'s own regulator normalization, which is
scale-invariant by construction — that's what makes `pe.cost` values comparable across a model's
absolute cost units in the first place.

## 10. Plotting

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

`plot_identifiability_graph(ia)` draws the same correlation structure as a network instead of a
heatmap: nodes are parameters (colored by `ia.index`), and an edge connects two parameters whose
`|correlation| > 0.5` — a cluster of mutually-connected parameters reads as "the data constrains
some combination of these, not each individually." When `ia.correlation_reliable` is `False`, it
renders as isolated nodes rather than a misleading fully-connected graph.

Also available: `plot_sensitivity_area` (time-dependent indices, stacked), `plot_identifiability_index`
(sorted index bar chart), and `plot_mcf` (prior/low/high ECDF panels, one per free parameter).

## 11. The semi-automation identification cycle

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

## Moving between MATLAB and Python

The two ports track the same project, but they are not drop-in translations of each other. This
section collects every difference that has actually caught someone out. If you only read one
thing here, read the first table.

### Traps — same call, different meaning

| | MATLAB | Python | What happens if you carry it across |
|---|---|---|---|
| **Argument order** | `gsua_sa(M, T)`, `gsua_ua(M, T)` | `sensitivity_analysis(model, M)`, `uncertainty_analysis(model, M)` | Model and design matrix are **swapped**. A swap now raises a clear `TypeError` naming the correct order, rather than silently returning nonsense. |
| **`margin` in profile likelihood** | `gsua_likelihood(..., margin=1.08, ...)` — the assumed relative std is `margin - 1` | `profile_likelihood(..., margin=0.08)` — the assumed relative std *is* `margin` | **Offset by one.** Passing MATLAB's `1.1` to Python asserts 110% noise; passing Python's `0.1` to MATLAB asserts 90% noise *and* flips `gsua_pe`'s internal `+1` offset positive, silently switching its inner refit from the likelihood to plain least squares. Use `1+m` in MATLAB where Python takes `m`. |
| **`margin` sign in estimation** | `gsua_pe` with `margin < 0` (after its `+1` offset) selects a Gaussian negative-log-likelihood objective | `parameter_estimation` takes `abs(margin)` and always uses the regulator cost | Same input, **different objective**, no error. |
| **Repeated-estimate orientation** | `gsua_ia(T, T_est)` with `T_est` as Np × N — one **column** per run | `identifiability_analysis(model, estimates)` with `(N, Np)` — one **row** per run | Transposed. Consistent with `PEResult.x`, but a square case (Np == N) fails silently. |
| **Joint-sampling pool orientation** | `Tia.Est` is Np × nPool | `pool=` is `(n_pool, Np)` | Transposed, for the same reason. |

### Conventions that differ by design

- **Fixed factors.** MATLAB drops them from `T`, so the table shrinks and indices shift. Python
  keeps them, with `model.fixed` marking them, so indices stay stable across a `fix()` call.
- **Output selection.** MATLAB sets `T.Properties.CustomProperties.output`; Python sets
  `model.output` / `model.set_output(...)`. Same idea, and both are honoured everywhere.
- **Results.** MATLAB writes back into the table (`T.Range`, `T.Est`, `T.Si`, ...), so calls chain
  naturally. Python returns a dataclass per routine, so to feed one result into the next you
  assign it back yourself (`model.range = ia.range`).
- **Plotting.** MATLAB's analysis functions open figures as a side effect — `gsua_ua` always opens
  at least two, `gsua_ia` four, and `gsua_eval` plots by default unless you pass its sixth
  argument as `false`. Python never plots implicitly; the `plot_*` functions take a result object.
- **Log-scale search** (`model.log_scale`) is Python-only. **Simulink models** are MATLAB-only.

### Defaults that differ

| Concept | MATLAB | Python |
|---|---|---|
| Save results to disk | `gsua_pe` writes `Estimations.mat` to the working directory (`'save'` defaults true) | opt-in: `parameter_estimation(..., save_path=...)` writes an `.npz`; nothing is written otherwise |
| Progress printing | `gsua_pe`/`gsua_ia` print unconditionally; `gsua_eval` prints a progress line | silent, except `warnings.warn` |
| Parallelism | `gsua_sa`/`gsua_ua` default `'parallel', true` | `sensitivity_analysis`/`uncertainty_analysis`/`Model.evaluate_batch` take `n_jobs` (default `1`, serial; `-1` all cores) |
| Default evaluation grid | `gsua_eval` expands `domain` to a unit-step grid | `model.domain` is used verbatim — give it the points you want, not just the endpoints |
| Plot on sampling | `gsua_dmatrix(..., 'Show','on')` draws a scatter | `design_matrix` returns `M` only |

### Options with no Python equivalent

`gsua_pe`: `A`, `B`, `Aeq`, `Beq`, `nonlcon` (linear and nonlinear constraints have no path at
all), `Multistart`, `Show`, `timer`; solvers `particle`, `psearch`, `surrogate`. (`save` is
covered by `save_path`, and negative `margin` selects the same Gaussian-NLL objective it does in
MATLAB.)
`gsua_likelihood`: `reps`, `parallel`, `show`, `saver` (so a long Python profile is unresumable).
`gsua_sa`: `bandwidth`, and the `brute-force` method. `gsua_oatr`/`gsua_csb`: `titerlimit`,
`parallel`, `show`, `breaking`, `stretch`. `gsua_costcutoff` has no standalone Python function —
the logic is reachable only through `identifiability_analysis`'s `cost=` argument.

### Solver names

| MATLAB `'solver'` | Python `solver=` |
|---|---|
| `lsqc` (default), `lsqn` | `least_squares` (default) |
| `fmincon` | `minimize` |
| `ga` | `differential_evolution` |
| `annealing` | `dual_annealing` |
| `particle`, `psearch`, `surrogate` | *not ported* |

### Numerical differences to expect

- **Profile likelihood.** Both use the same acceptance threshold (chi-squared at one degree of
  freedom, halved). Once `margin` is matched per the table above, the two agree closely. Python
  refits with a single local `minimize` where MATLAB runs a `reps`-restart `gsua_pe`.
- **Range refinement.** MATLAB accepts a bound inside a 10% dead band; Python root-finds the exact
  crossing, so its bounds are systematically tighter by up to that slack.
- **Percentile convention.** `gsua_covmetric` uses MATLAB's Hazen percentiles; `coverage_metric`
  uses NumPy's linear default. The 5/95 bands differ slightly for small ensembles.
- **NaN handling.** `rcostf` is NaN-tolerant in Python but not in MATLAB, so gappy data gives
  different objectives across the two.

### Name reference

Every function below is also available under its idiomatic Python name (e.g. `costf` for
`gsua_costf`) — see [README.md](README.md#design) for the naming convention. The `gsua_*` aliases
match the MATLAB spelling exactly, capitals included: `gsua_MCF`, `gsua_medianCI`,
`gsua_costfMulti`.

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
| `gsua_noisefloor` | `noise_floor` |
| `gsua_plot` | `plot_uncertainty`, `plot_sensitivity_bar`, `plot_sensitivity_area`, `plot_identifiability_correlation`, `plot_identifiability_graph`, `plot_identifiability_index`, `plot_mcf` |
| `gsua_dpmat`, `gsua_odefun` | `SymbolicODEModel` |
| `gsua_userdefined` | `UserFunctionModel` |
| `gsua_eval`, `gsua_deval`, `gsua_pardeval` | `Model.evaluate`, `Model.evaluate_batch` |
| `gsua_dataprep` | the `Model` constructors above |
| `gsua_costcutoff` | *(no standalone equivalent; use `identifiability_analysis(cost=...)`)* |

## Getting help

Found a bug, or a place where this port's behavior diverges from the MATLAB original in a way
that isn't documented in the relevant module's docstring? Please open an issue at
[github.com/drojasd/GSUA-CSB](https://github.com/drojasd/GSUA-CSB/issues).
