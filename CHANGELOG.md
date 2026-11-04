# Changelog

## 1.17.1 (MATLAB) / 0.3.1 (Python) — 2026-11-04

Patch release: three fixes, with MATLAB/Python parity and tests in both languages.

### Fixed

- **Transposed uncertainty plots when N == numel(xdata).** `gsua_plot('UncertaintyAnalysis', ...)`
  called `plot(t, Y(:,:,j))` with a simulations×time matrix. When that matrix is square, MATLAB
  plots its columns (time instants), so the ensemble came out transposed — each drawn curve was one
  time instant sweeping the simulations. Each `Y(:,:,j)` is now transposed explicitly so every line
  is one simulation regardless of shape, and the curves-vs-boxplot choice uses `numel(t)` (not
  `size(t,2)`), so a column `xdata` still draws curves. Affects `gsua_eval`/`gsua_ua`/`gsua_csb`,
  which share this branch. (Python's `plot_uncertainty` iterates rows explicitly and was already
  immune; a square-ensemble regression test was added on both sides.)
- **`.mltbx` pinned a maximum MATLAB release.** The package declared `releaseCompatibility
  end="R2025b"` (from `param.release.end="latest"`, resolved to the packaging release), so installing
  on a newer release or on MATLAB Online reported "not supported for this MATLAB release" even though
  the toolbox runs. `buildToolbox` now leaves `MaximumMatlabRelease` empty — the documented "no
  maximum restriction" — which serializes to `end="latest"` (MATLAB's unbounded sentinel, not a
  concrete release). The R2017a minimum and MATLAB Online compatibility are kept; a test inspects the
  built artifact's `systemRequirements.xml`.
- **`gsua_bootstrap` was slow by default.** With no `'opt'`, the B internal refits ran the solver
  with Display on, so `fmincon`/`ga`/… printed a result block per start; on an interactive desktop
  that console rendering dominated the runtime (a 40×8 fmincon bootstrap measured ~234 s versus
  ~54 s with Display off). `gsua_bootstrap` now builds Display-off options for the chosen solver by
  default; Display is cosmetic, so the ensemble for a given seed is unchanged. (Python's refits use
  SciPy, which is silent by default — already at parity.)

## 1.17.0 (MATLAB) / 0.3.0 (Python) — 2026-11-04

Overhaul of the `gsua_pe` → `gsua_noisefloor` → `gsua_ia` workflow (course-driven), with MATLAB and
Python parity, tests in both languages, and updated docs and course module 5.

### New

- **`gsua_bootstrap` / `parametric_bootstrap`** — parametric-bootstrap ensemble of parameter
  estimates for identifiability analysis: from the best fit it draws `B` synthetic datasets under the
  data's own count-noise model (`poisson`/`quasipoisson`/`nbglobal`, with cumulative-series support),
  re-fits each with `gsua_pe`, and returns the ensemble (seed, progress, optional parallel). Unlike a
  multistart pool (which collapses to one optimum), its spread is honest parameter uncertainty — pass
  it to `gsua_ia`.
- **`gsua_clustersummary` / `cluster_summary`** — per-cluster, per-parameter table: n, median and
  distribution-free median CI, mean and Student-t CI, 2.5–97.5 percentile interval, that interval's
  width relative to the original range, and an identifiability verdict (<25% identifiable, >50% not,
  else dubious).
- **Per-cluster results in `clusterInfo`** — `Range{c}`/`Index{c}`/`Est{c}`/`Median{c}` for every
  cluster (computed like the dominant), not only the dominant one.
- **Cost-filter mask in `clusterInfo`** — `Kept` (survivors over the original `T_est` columns) and
  `Idx`/`labels_full` over the original columns (`NaN` for cost-dropped runs).
- **`SilhouetteByK`/`SilhouetteK` in `clusterInfo`** and a silhouette-vs-K plot with the
  `silThreshold` line.
- **Per-cluster visualizations** — pairwise scatter matrix colored by cluster, and per-cluster
  curves vs the data (`gsua_ia` `'xdata'`/`'ydata'`; Python `plot_cluster_pairs`/`plot_cluster_curves`).
- **`'truth'` overlay** in `gsua_ia` (and Python `plot_parameter_cdfs`/`plot_cluster_pairs`) marking
  reference values on CDFs, boxplots and scatter — for synthetic/teaching data.
- **Threshold-selection plots** — `gsua_noisefloor(..., 'Show', true)` / `plot_noise_floor`: sorted
  pool costs (log scale) with the threshold and accepted/rejected runs, plus the bootstrap cost
  distribution with the chosen quantile.
- Shared helpers `gsua_noisemodel`, `gsua_noisedraw`, `gsua_iastats` (one implementation behind
  `gsua_noisefloor`/`gsua_bootstrap` and `gsua_ia`/`gsua_dia`).

### Changed

- **Identifiability index for flat parameters** — the interval-width term now defaults to the
  estimate 2.5–97.5 percentile dispersion over the **original** range (`IndexMethod='spread'`), so a
  parameter scattered across its whole prior is no longer rated as identified just because its median
  is well pinned. The previous median-CI behaviour is kept as `IndexMethod='ci'`.
- Course module 5 (MATLAB + Python) now uses `gsua_bootstrap` and `gsua_clustersummary`.
- Documented the `margin` convention (a relative residual band; used internally as `m+1`; rescales
  the cost, not the optimum) in `gsua_pe`/`gsua_costf`/`gsua_noisefloor` and the Python equivalents.

### Fixed

- **`gsua_pe` solver dispatch** — `'psearch'` passed `Np` instead of a start point, and `'annealing'`
  omitted `x0` (shifting every argument); both now pass `IP(i,:)` (and `IP` is generated for them).
  Added a smoke test for all eight MATLAB solvers / four Python solvers.
- **Single-output orientation** — `gsua_noisefloor`/`gsua_bootstrap`/`gsua_noisemodel` now orient
  `ydata` to outputs×samples, so a single output passed as a column (`len×1`, as produced by
  `squeeze(gsua_eval(...)).'`) no longer mistakes `len` for the output count. (Python was already
  immune via `atleast_2d`.)
