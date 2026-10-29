%[text] # GSUA-CSB Toolbox — user guide
%[text] **Global Sensitivity and Uncertainty Analysis — Confidence Subcontour Box (GSUA-CSB)** is a command-line toolbox, developed at Universidad EAFIT, for validating mathematical models built in Simulink, with the Symbolic Math Toolbox, or as plain MATLAB functions. It performs variance-based sensitivity analysis (SA) and uncertainty analysis (UA), estimates model parameters (PE), analyzes their practical identifiability, and estimates confidence subcontour boxes (CSB) for the estimates. It builds on the earlier work of Carlos Mario Vélez, *GSUA of dynamical systems using variance-based methods*.
%[text] **Developers:** Daniel Rojas Díaz (drojasd@eafit.edu.co) · Carlos Mario Vélez Sánchez (cmvelez@eafit.edu.co).
%[text] **How this guide is organized.** It opens with the one function every workflow starts from — `gsua_dataprep`, which turns a model into the summary table `T` — and then walks the pipeline: simulate, sample, UA, SA, estimate, analyze identifiability, and size confidence intervals. Each function has its own section with its purpose, its call syntax, and a complete list of its **name–value arguments with defaults**. The three model-preparation workflows (Simulink, symbolic, user-defined) are detailed in the appendix.
%[text] Every function also has standard help: type `help gsua_dmatrix` (or any function name) at the command line for the same reference, and `open gsua_userguide` to reopen this guide.
%[text:tableOfContents]{"heading":"**Contents**"}
%%
%[text] ## 0. Preparing a model — `gsua_dataprep`
%[text] `gsua_dataprep` is the entry point: it consolidates a model and its factor ranges into the summary table `T` that every other function consumes. It dispatches on the first argument, so the call differs by model kind:
%[text] - **Simulink model** — `[T,tout] = gsua_dataprep(model_name, Ranges, Parameter_names)` \
%[text] - **Symbolic ODE system** — `[T,solver] = gsua_dataprep(odes, vars, domain, model_name)` \
%[text] - **User-defined function** — `[T,sol] = gsua_dataprep(func, Ranges, 'domain', domain)` \
%[text] A *factor* is any input, parameter, or initial condition under analysis. `Ranges` is an `Np x 2` array of `[lower upper]` bounds, one row per factor; give a factor a degenerate range (`lower == upper`) to **fix** it — fixed factors drop out of `T`. The appendix shows each workflow end to end, including the name–value arguments specific to each kind (`'rMethod'`, `'nominal'`, `'output'`, `'names'`, `'out_names'`, `'solver'`, `'opt'`, `'vectorized'`, `'modelkind'`, `'Step'`).
%[text] For the worked examples below we use the bundled symbolic model, loaded with:
default_example;
T
%%
%[text] ## 1. Simulate the model — `gsua_eval`
%[text] `gsua_eval` runs the model for one or a few factor vectors and returns the output, interpolated onto `xdata` when given. It is the quickest way to check that a freshly prepared `T` actually simulates, and to overlay candidate factor families against data. For many simulations (Monte-Carlo sweeps) use `gsua_pardeval` instead.
%[text] Call syntax: `[y,xdata] = gsua_eval(values, Table, xdata, ydata, parallel, show, verbose)`. `values` is `Np x N` (factors by simulations); the trailing arguments are **positional**, not name–value:
%[text] - **`xdata`** (default `[]`) — time points to interpolate the output onto. Empty uses the model's own domain grid. \
%[text] - **`ydata`** (default `[]`) — reference/experimental data drawn alongside the model output for comparison. \
%[text] - **`parallel`** (default `false`) — evaluate the columns of `values` in parallel. \
%[text] - **`show`** (default `true`) — draw the result with `gsua_plot`; pass `false` to suppress the figure. \
%[text] - **`verbose`** (default `true`) — print the `Progress: N %` line; pass `false` for a quiet run. \
% A single simulation at the nominal factor values stored in T.Nominal:
y = gsua_eval(T.Nominal, T, [], [], false, false, false);
% Build a noisy synthetic dataset to reuse in later sections:
xdata = linspace(0, 264, 265);
yclean = gsua_eval([1000 200 0.1 0.06]', T, xdata, [], false, false, false);
yexp = yclean .* (1 + 0.2*randn(size(yclean)));
%[text] To change how many model outputs are exposed (symbolic and user-defined models), set the `output` custom property before evaluating, e.g. `T.Properties.CustomProperties.output = 1:3;`.
%%
%[text] ## 2. Sample the factor space — `gsua_dmatrix`
%[text] `gsua_dmatrix` draws the design matrix `M` (`N x Np`) that feeds UA and SA. Its second output is a copy of `T` with any newly fixed factors updated.
%[text] Call syntax: `[M,T2] = gsua_dmatrix(T, N, ...)`.
%[text] **Name–value arguments**
%[text] - **`'Method'`** (default `'LatinHypercube'`) — sampling scheme: `'LatinHypercube'`, `'Uniform'`, `'Sobol'`, or `'Joint'`. The first three are *marginal*: each factor is drawn independently inside its own range. `'Joint'` draws whole factor vectors from an ensemble of accepted estimates (see below). \
%[text] - **`'Show'`** (default `'off'`) — set to `'on'` to draw the scatter matrix of sampled factor pairs. \
%[text] - **`'Seed'`** (default `[]`) — RNG seed for reproducible draws. \
%[text] The next four arguments apply only to `'Method','Joint'`:
%[text] - **`'JointType'`** (default `'Bootstrap'`) — how the ensemble becomes samples. `'Bootstrap'` resamples whole pool vectors with replacement (exact correlations, no distributional assumption, cannot leave the identified manifold, but limited to the distinct vectors in the pool). `'SmoothBootstrap'` adds a variance-corrected Gaussian kernel so draws are new points near the manifold. `'Gaussian'` uses the ensemble mean and covariance (handles a rank-deficient covariance, but assumes a locally linear manifold — on a curved ridge it places mass off the manifold). \
%[text] - **`'Pool'`** (default `[]`) — supply the ensemble explicitly as an `Np x nPool` matrix instead of reading `Tia.Est`. \
%[text] - **`'MinPoolN'`** (default `5`) — smallest pool accepted; a smaller pool is refused (correlation from a handful of points is degenerate — two points always correlate at ±1) and the call falls back to marginal Latin-hypercube sampling with a warning. \
%[text] - **`'Bandwidth'`** (default `[]`) — kernel bandwidth for `'SmoothBootstrap'`; empty uses Silverman's multivariate rule. \
%[text] - **`'Clip'`** (default `[]`, no clipping) — `Np x 2` bounds to clip joint draws back into. Left empty on purpose: `gsua_ia` narrows `T.Range` to a confidence interval of the *median*, so clipping joint draws to it would re-impose a narrowness that understates uncertainty. Pass the original model's physical bounds if clipping is wanted. \
M = gsua_dmatrix(T, 500, 'Show', 'on');
%[text] **Marginal vs joint.** When factors are correlated — which `gsua_ia`'s correlation heatmap shows, and which a `cost_band` far above `cost_data` in `gsua_covmetric` signals — independent (marginal) draws leave the identified manifold, so the band built from them is both too wide and badly centered. Joint sampling preserves the correlation structure:
%%
%[text]
%   Tia = gsua_ia(T, T.Estlsqc);                 % Tia.Est holds the accepted ensemble
%   Mjoint = gsua_dmatrix(Tia, 2000, 'Method','Joint', 'JointType','Bootstrap');
%
%%
%[text] ## 3. Uncertainty analysis — `gsua_ua`
%[text] `gsua_ua` propagates the sampled factors through the model and summarizes the spread of outputs. It automatically runs Monte-Carlo filtering (`gsua_MCF`) and draws it: `figure(1)` holds the output band, and one further figure per output shows, factor by factor, the empirical CDFs of the runs that landed below vs above the reference output.
%[text] Call syntax: `Y = gsua_ua(M, T, ...)`.
%[text] **Name–value arguments**
%[text] - **`'xdata'`** (default `[]`) — time points to evaluate/interpolate onto; empty uses the model's own domain grid. \
%[text] - **`'ynom'`** (default `[]`) — reference output the filtering splits runs against; empty uses the output at `T.Nominal`. Its length must match `xdata`. \
%[text] - **`'parallel'`** (default `true`) — evaluate the design matrix in parallel. \
%[text] - **`'verbose'`** (default `true`) — print the `Progress: N %` line. \
%[text] - **`'plot'`** (default `true`) — draw the figures. Set `false` in batch runs that call `gsua_ua` many times, where the figures would otherwise pile up; the returned `Y` is identical either way. \
Y = gsua_ua(M, T, 'xdata', xdata, 'parallel', false);
%[text] **Reading the filtering plots.** The larger the gap between the low-value (blue) and high-value (red) curves for a factor, the more that factor matters. If the blue curve is on the bottom, high factor values push the output down; if the red curve is on the bottom, high factor values push the output up. When a specific `'ynom'` is given the interpretation gets subtler — a sensitivity analysis (next section) is the cleaner tool.
%%
%[text] ## 4. Sensitivity analysis — `gsua_sa`
%[text] `gsua_sa` estimates variance-based sensitivity indices: first-order `Si` (a factor's own contribution to output variance) and total `STi` (its contribution including interactions). It can return two extra outputs for later analysis.
%[text] Call syntax: `[T,Y,J] = gsua_sa(M, T, ...)` — note the order: `Y` is the full output matrix (as from `gsua_ua`) and `J` is its scalar MSE-vs-reference representation.
%[text] **Name–value arguments**
%[text] - **`'SensMethod'`** (default `'Xiao'`) — estimator: `'Xiao'` (newest, good for time-dependent responses, `N*(Np/2+1)` simulations), `'Saltelli'` and `'Jansen'` (same cost), `'Sobol'` (classic, `N*(Np+1)`), `'OAT'` (one-at-a-time, non-global, `N*(Np+1)`), or `'brute-force'` (all combinations, `N+2*Np*N^2`). \
%[text] - **`'ynom'`** (default `[]`) — reference output for index estimation; empty uses the output at `T.Nominal`. Index values do not depend on this choice. Its length must match `xdata`. \
%[text] - **`'xdata'`** (default `[]`) — resolution `[init:step:final]` for the time-dependent indices; empty uses the model's domain grid. \
%[text] - **`'parallel'`** (default `true`) — run the simulations in parallel. \
%[text] - **`'verbose'`** (default `true`) — print progress. \
%[text] - **`'xAlpha'`** (default `1`) — exponent in the distance-component estimator, in `(0,2]`. \
%[text] - **`'bandwidth'`** (default `1`) — bandwidth for the distance-component estimator. \
[T,Y,J] = gsua_sa(M, T, 'SensMethod', 'Xiao', 'parallel', false);
%[text] **Visualizing the indices.** `gsua_plot` draws them in several forms (see also the scatter options for `gsua_eval`/UA results):
%[text] - `gsua_plot('Pie',T,T.STi)` — pie of total indices. \
%[text] - `gsua_plot('Bar',T,T.STi)` — bar chart of total indices. \
%[text] - `gsua_plot('TotalSensitivityArea',T,T.STi_vec,xdata)` — normalized area of the time-dependent total indices. \
%[text] - `gsua_plot('FractionalSensitivityArea',T,T.Si_vec,xdata)` — area of the time-dependent first-order indices. \
%[text] - `gsua_plot('ScatterOutput',T,J,M,xdata)` — scatter of the scalar output against each factor. \
clf
gsua_plot('Bar', T, T.STi)
%[text] **Interpreting the indices.** `STi` is the key measure of factor relevance. The sum of `STi` is always `>= 1` and the sum of `Si` always `<= 1`; if either equals 1 so does the other, and the model is purely additive. As a rule of thumb, if `sum(Si) > 0.65` the factor interactions are weak enough to treat as negligible, and local or analytical SA becomes appropriate. For the Saltelli and distance-based estimators, negative `Si` values are the signal that `N` was too small — raise it.
%%
%[text] ## 5. Parameter estimation — `gsua_pe`
%[text] `gsua_pe` fits the model to data by multistart optimization, returning the table with the estimates attached (`T.Est<solver>`, e.g. `T.Estlsqc`) and the per-run cost `res`.
%[text] Call syntax: `[T,res] = gsua_pe(T, xdata, ydata, ...)`.
%[text] **Name–value arguments**
%[text] - **`'N'`** (default `1`) — number of estimation runs (multistart restarts from different points). \
%[text] - **`'solver'`** (default `'lsqc'`) — optimizer: `'lsqc'` (lsqcurvefit), `'lsqn'` (lsqnonlin), `'fmincon'`, `'ga'`, `'particle'` (particleswarm), `'psearch'` (patternsearch), `'surrogate'`, or `'annealing'` (simulannealbnd). \
%[text] - **`'opt'`** (default `[]`) — an `optimoptions` object for the chosen optimizer (the optimizer name must match). \
%[text] - **`'Multistart'`** (default `1`) — number of estimations per cycle using MATLAB's MultiStart (lsqcurvefit only). \
%[text] - **`'ipoint'`** (default `[]`) — explicit initial points, with `N` rows and one column per free factor, instead of random starts. \
%[text] - **`'margin'`** (default `0`) and **`'alpha'`** (default `2`) — switch the objective from plain least squares to the correlation-penalized multi-objective cost (`gsua_costf`), which is what you want when one model must match several signals at once. `gsua_pe` records both on the table so `gsua_noisefloor` (section 7) recovers them automatically. \
%[text] - **`'A'`, `'B'`, `'Aeq'`, `'Beq'`, `'nonlcon'`** (default `[]`) — linear (`A*x<=B`, `Aeq*x=Beq`) and nonlinear constraints passed through to the optimizer. \
%[text] - **`'Show'`** (default `'off'`) — `'on'` plots the optimizer results; `'on2'` also calls `gsua_ia`. \
%[text] - **`'save'`** (default `false`) — if `true`, writes `Estimations.mat` to the current folder after each run as a crash-resume checkpoint, overwriting any existing file of that name. (It defaulted to `true` before v1.15.) \
%[text] - **`'timer'`** (default `true`) — print per-estimation progress. \
opt = optimoptions('lsqcurvefit', 'Display', 'off');
[T,res] = gsua_pe(T, xdata, yexp, 'solver', 'lsqc', 'N', 10, 'opt', opt, 'timer', false);
T
%%
%[text] ## 5.1 Practical identifiability — `gsua_ia` and `gsua_dia`
%[text] `gsua_ia` turns a pool of repeated estimates into identifiability diagnostics: per-factor CDFs and boxplots, a factor-correlation heatmap, and an identifiability index. It replaces `T.Range` with a distribution-free confidence interval of the median (`gsua_medianCI`), sets `T.Nominal` to the median estimate, `T.index` to the per-factor index (0 = well identified, 1 = poorly identified), and `T.Est` to the estimates actually used. `gsua_dia` computes the same statistics without drawing any figures — use it in scripted pipelines.
%[text] Call syntax: `[T,clusterInfo] = gsua_ia(T, T_est, correction, outlier, show, cluster, maxK, silThreshold, ...)`. The first arguments are **positional**:
%[text] - **`correction`** (default `false`) — clip the new range back to the original `T.Range` when it would otherwise widen. \
%[text] - **`outlier`** (default `false`) — remove multivariate outliers before computing statistics. Note the bundled detector is aggressive on clean, tightly clustered data; inspect `T.Est` if results surprise you. \
%[text] - **`show`** (default `false`) — plot the identifiability-index scatter, correlation graph, and (if `cluster` is true) the clustering diagnostics. \
%[text] - **`cluster`** (default `false`) — run spectral clustering to detect multiple global minima (distinct basins the optimizer landed in) rather than assuming one. Requires Statistics and Machine Learning Toolbox. \
%[text] - **`maxK`** (default `5`) — largest cluster count tested. \
%[text] - **`silThreshold`** (default `0.5`) — minimum mean silhouette score to accept a multi-cluster split over a single basin. \
%[text] **Name–value arguments** (fit-quality filtering, shared with `gsua_costcutoff`):
%[text] - **`'cost'`** (default `[]`) — the `res` output of `gsua_pe`, one entry per column of `T_est`; screens out runs that converged badly before any statistic is computed. \
%[text] - **`'CostMethod'`** (default `'rtol'`) — `'rtol'` (keep `cost <= min*(1+CostRtol)+CostAtol`) or `'gap'` (cut at the largest jump in the sorted costs). \
%[text] - **`'CostRtol'`** (default `0.1`), **`'CostAtol'`** (default `1e-8`), **`'CostGapRatio'`** (default `3`), **`'MinKeep'`** (default `3`) — tolerances for the two methods and a floor on how many runs survive. \
T_newRange = gsua_ia(T, T.Estlsqc)
%[text] When `cluster` is true and a genuine split is found, every statistic reflects the dominant (largest) cluster only — pooling separated basins into one interval is not meaningful — and `clusterInfo` reports each basin's candidate global minimum:
%%
%[text]
%   [Tia,clusterInfo] = gsua_ia(T, T.Estlsqc, false, false, true, true, 'cost', res);
%   clusterInfo.NumClusters   % > 1 means distinct basins were found
%   clusterInfo.Centers       % one candidate global minimum per basin
%
%%
%[text] ## 6. Uncertainty-interval calculus — `gsua_oatr` and `gsua_csb`
%[text] These size a confidence interval around an estimate by sampling. `gsua_oatr` expands each factor's range one at a time around the objective output stored in `T.Nominal`; `gsua_csb` then estimates the confidence subcontour box by iterated sampling. `gsua_oatr2` is a second-generation one-at-a-time refinement.
%[text] `gsua_oatr` call syntax: `T = gsua_oatr(T, ...)`.
%[text] - **`'lim'`** (default `0.3`) — distance limit for range expansion. \
%[text] - **`'parallel'`** (default `false`) — parallelize the range calculus. \
%[text] - **`'correct'`** (default `false`) — clip expansions to the original range. \
%[text] - **`'show'`** (default `false`) — plot the iterations. \
%[text] - **`'titerlimit'`** (default `inf`) — cap on total iterations (`gsua_oatr` only). \
%[text] `gsua_csb` call syntax: `[T,New_range,J_test,Y_test,sup] = gsua_csb(T, N, ...)`, where `N` is samples per cycle. **Do not run `gsua_csb` inside a Live Script** — run it in the command window.
%[text] - **`'reps'`** (default `100`) — number of cycles. \
%[text] - **`'stop'`** (default `0.95`) — confidence level of the interval. \
%[text] - **`'lim'`** (default `0.3`) — distance criterion. \
%[text] - **`'select'`** (default `0.5`) — selection criterion. \
%[text] - **`'switch'`** (default `0.75`) — cut criterion. (Older help called this `'recort'`; the actual argument is `'switch'`.) \
%[text] - **`'ynom'`** (default `[]`) — a new objective output. \
%[text] - **`'correct'`** (default `true`), **`'parallel'`** (default `false`), **`'show'`** (default `false`) — clip to bounds, parallelize, plot. \
T.Nominal = [1000 200 0.1 0.06]';
T_oat = gsua_oatr(T);
%%
%[text]
%   T_csb = gsua_csb(T_oat, 100)   % run in the command window, not here
%
%%
%[text] ## 7. Fit quality, acceptance, and confidence intervals
%[text] Sections 5 and 6 give an estimate and a range. The functions here answer what comes next: was the fit good enough to believe, which multistart runs to keep, how wide the acceptance band should be, and what the confidence interval of each factor is once the others are free to move.
%%
%[text] ## 7.1 Band coverage and tightness — `gsua_covmetric`
%[text] `gsua_covmetric` scores a Monte-Carlo band against the data with two numbers on the same normalized scale (below 1 = within tolerance): `cost_data`, the distance from the band median to the data (accuracy), and `cost_band`, the distance between the band edges (precision, independent of the data). A `cost_band` far larger than `cost_data` is the signature of factor confounding — the band is wide because the sampled factors are inconsistent with each other, which is exactly what the joint sampling in section 2 fixes.
%[text] Call syntax: `[cost_data,cost_band,P5,P50,P95] = gsua_covmetric(Y, ydata, ...)`.
%[text] - **`'margin'`** (default `0.1`) — relative tolerance for the shared normalization (values `<1` are read as a fraction and converted to `1+margin`). \
%[text] - **`'alpha'`** (default `2`) — exponent on each per-signal cost term. \
%[text] - **`'plevels'`** (default `[5 95]`) — `[low high]` percentile levels defining the band. \
%%
%[text]
%   M = gsua_dmatrix(T,2000);
%   Y = gsua_ua(M,T,'xdata',xdata);
%   [cost_data,cost_band] = gsua_covmetric(Y,yexp,'margin',0.1);
%
%%
%[text] ## 7.2 Which runs to keep — `gsua_costcutoff`
%[text] A multistart run that converged to a much worse cost than the best one contributes essentially arbitrary factor values; left in the pool it can make a well-identified factor look poorly identified. `gsua_costcutoff` decides which runs to keep. It is the filter `gsua_ia`/`gsua_dia` apply through their `'cost'` argument, exposed on its own.
%[text] Call syntax: `keep = gsua_costcutoff(cost, method, rtol, atol, gapRatio, minKeep)` — all six arguments are positional and required (the `gsua_ia` defaults are `'rtol', 0.1, 1e-8, 3, 3`). `method` is `'rtol'` (keep `cost <= min*(1+rtol)+atol`) or `'gap'` (cut at the largest jump in the sorted log-costs, only if it is at least `gapRatio` times the median gap). `minKeep` guarantees the pool is never emptied.
%%
%[text]
%   keep = gsua_costcutoff(res,'rtol',0.1,1e-8,3,3);
%   Est_kept = T.Estlsqc(:,keep);
%
%%
%[text] ## 7.3 A data-calibrated acceptance threshold — `gsua_noisefloor`
%[text] The common rule `res < 1.5*res(1)` is scale-dependent: it asks how close a run is to the best run, not whether it is statistically distinguishable from the truth, so it tightens as the fit improves and can reject curves that describe the data perfectly well. `gsua_noisefloor` replaces it with a parametric bootstrap: it simulates the best fit, estimates an observation-noise model from that fit's residuals, generates synthetic datasets under that noise, and scores the true model against each with the same cost `gsua_pe` minimized. The upper quantile of that distribution is a data-calibrated acceptance threshold. Because `gsua_pe` records the `margin` and `alpha` it used on the table, `gsua_noisefloor` recovers them automatically.
%[text] Call syntax: `out = gsua_noisefloor(T, xdata, ydata, Est, res, ...)`.
%[text] - **`'B'`** (default `2000`) — bootstrap sample size. \
%[text] - **`'noise'`** (default `'quasipoisson'`) — noise model: `'poisson'` (usually far too tight for real data), `'quasipoisson'` (recommended; variance `phi*mu`), or `'nbglobal'` (single pooled negative-binomial `k`, report as a bound only). `out.byModel` reports all three regardless, as a sanity check. \
%[text] - **`'quantile'`** (default `0.95`) — upper quantile of the bootstrap cost distribution that defines the threshold. \
%[text] - **`'cumulative'`** (default all `false`) — one logical flag per output row; a cumulative series has its noise estimated on the first-differenced (incident) series and re-accumulated. \
%[text] - **`'seed'`** (default none) — seed for reproducible draws; the caller's global RNG state is always restored on return. \
%[text] - **`'margin'`**, **`'alpha'`** (default recovered from the table) — override the values `gsua_pe` scored with, in `gsua_pe`'s own convention. \
%[text] Key outputs: `out.threshold` (the cutoff on the `res` scale), `out.accepted` (logical `res < threshold`), and `out.byModel` (all three noise models). If the best fit itself exceeds the floor, a warning fires — the model is rejected as a description of the data at that quantile.
%%
%[text]
%   [T,res] = gsua_pe(T,xdata,yexp,'solver','lsqc','N',50,'margin',0.1);
%   out = gsua_noisefloor(T,xdata,yexp,T.Estlsqc,res);
%   Tia = gsua_ia(T, T.Estlsqc(:,out.accepted));
%
%%
%[text] ## 7.4 Profile-likelihood confidence intervals — `gsua_likelihood`
%[text] `gsua_likelihood` scans each factor across a range while re-estimating all the others at every step, and reports the interval where the profiled cost stays below a likelihood-ratio threshold. It is slower than the distribution-free interval `gsua_ia` reports, but because it re-fits the remaining factors it accounts for their trade-offs — the better choice precisely when factors are correlated.
%[text] Call syntax (all arguments **positional**): `Tfinal = gsua_likelihood(T, xdata, ydata, alpha, step, margin, tolerance1, tolerance2, limit, reps, show, parallel, saver, pars, verbose)`. Here `alpha` is the confidence level (e.g. `0.95`), `step` the initial relative step, `margin` the relative standard deviation in the Gaussian likelihood weighting, `tolerance1`/`tolerance2` the bisection stopping tolerances, `limit` the iteration cap per bound, `reps` the estimations per inner `gsua_pe` call, and `pars` the factor indices to profile (`[]` = all). `verbose` (default `true`) is optional.
%%
%[text]
%   Tci = gsua_likelihood(T,xdata,yexp,0.95,0.05,1.1,0.01,0.01,15,10,false,false,false,[]);
%
%%
%[text] ## 8. Supporting functions
%[text] The functions above are the main workflow. These fill it in:
%[text] - **`gsua_dia`** — the identifiability statistics of `gsua_ia` without the plots (positional args `outlier, cluster, maxK, silThreshold`, same `'cost'`/`'CostMethod'`/... name–value set). \
%[text] - **`gsua_oatr2`** — second-generation one-at-a-time range refinement (`'lim' 0.3`, `'parallel'`, `'correct'`, `'show'`). \
%[text] - **`gsua_medianCI`** — distribution-free confidence interval of the median, `[lo,hi] = gsua_medianCI(X, alpha)` with `alpha` default `0.05`; this is what `gsua_ia` uses to build `T.Range`. \
%[text] - **`gsua_depth`** — generalized band depth, `[depth,idx] = gsua_depth(x, parallel)`, for ranking curves by centrality; `gsua_sa` uses it to pick a reference output when no `'ynom'` is given. \
%[text] - **`gsua_MCF`** — Monte-Carlo filtering on its own, `[Jsup,Jinf] = gsua_MCF(T, M, Y, y_exp, tIn, tOut)`, outside of `gsua_ua`. \
%[text] - **`gsua_userdefined`** — build `T` straight from a user function, bypassing `gsua_dataprep` (same name–value arguments; see the appendix). \
%[text] - **`gsua_deval`** and **`gsua_pardeval`** — single and batched model evaluation used internally; `gsua_pardeval(M,T,xdata,parallel,init,verbose)`. \
%[text] - **`gsua_costf`, `gsua_rcostf`, `gsua_costfMulti`, `gsua_likecost`** — the cost functions `gsua_pe` minimizes. \
%%
%[text] ## 9. Python port
%[text] A Python port is maintained in the `python` folder (package `gsua_csb`). It mirrors these functions under idiomatic names, keeping `gsua_`-prefixed aliases for cross-reference: `design_matrix`/`gsua_dmatrix`, `parameter_estimation`/`gsua_pe`, `identifiability_analysis`/`gsua_ia`, `noise_floor`/`gsua_noisefloor`, and so on. Simulink models are MATLAB-only; symbolic ODE and user-defined models are both supported there. See `python/USERGUIDE.md` for its guide and `python/README.md` for what is and is not ported.
%%
%[text] ## Appendix A — Simulink models
%[text] Prepare a Simulink model as described in the Simulink setup guide, then:
%[text] `[T,tout] = gsua_dataprep(model, Ranges, Parameter_names)` — `model` is the model file name, `Ranges` an `Np x 2` array (or a cell array, see `'rMethod'`), and `Parameter_names` a cell/string array of factor names. Each factor must be assigned in the model workspace as a *reference to a named variable* (the toolbox matches factors by name); give a factor a constant value instead to exclude it. At present a Simulink model exposes a single output — connect an output block to the objective signal and terminators elsewhere.
%[text] **Name–value arguments**
%[text] - **`'rMethod'`** (default `'range'`) — how `Ranges` is read: `'range'` `{lower, upper}`, `'std'` `{nominal, distance}`, `'percent'` `{nominal, percent}`, or `'normal'` `{nominal, std}`. A third column in a cell `Ranges` sets a per-factor method (do not mix in `'normal'`). The normal method also shapes how `gsua_dmatrix` draws samples. \
%[text] - **`'nominal'`** (default midpoint of `Ranges`) — an explicit nominal factor family. \
%[text] - **`'modelkind'`** (default `'Dynamic'`) — set `'Static'` for static models (initial time = final time = 0). \
%[text] - **`'Step'`** (default `[]`) — a fixed solver step magnitude. \
%%
%[text]
%   model = 'predatorprey';
%   ParIn = {'beta','delta','ro','P0'};
%   Ranges = [0.4 0.6; 0.4 0.6; 2.2 2.8; 2.7 3.2];
%   [T,time] = gsua_dataprep(model, Ranges, ParIn);
%
%%
%[text] ## Appendix B — Symbolic ODE models
%[text] Build the system with the Symbolic Math Toolbox, choose an order for the state variables, then: `[T,solver] = gsua_dataprep(odes, vars, domain, modelName)`. `odes` is the equation array, `vars` the state-variable array, `domain` a `1 x 2` time span, and `modelName` the name of the model file the function writes for later simulations. If your system is in `y^(n) = f(...)` form, reduce it first with `reduceDifferentialOrder`.
%[text] **Name–value arguments**
%[text] - **`'range'`** / **`'rMethod'`** (default `'range'`) — provide the factor ranges and how they are read (`'range'`, `'std'`, `'percent'`, `'normal'`), as for Simulink. Fix a factor with a degenerate range. \
%[text] - **`'nominal'`** (default midpoint) — explicit nominal factor family. \
%[text] - **`'output'`** (default `1:model_order`) — which state(s) to expose (one output for SA and CSB). \
%[text] - **`'solver'`** (default `'ode45'`) — `'ode45'` (adaptive) or `'ode4'` (fixed step; add a third column to `domain` for the step). \
%[text] - **`'opt'`** (default none) — an `odeset` options struct (for `ode4`, only `NonNegative` is supported). \
%%
%[text]
%   syms S(t) I(t) R(t) beta gamma
%   P = S+I+R;
%   odes = [diff(S)==-S*I*beta/P; diff(I)==S*I*beta/P-gamma*I; diff(R)==gamma*I];
%   vars = [R S I];
%   Ranges = [0 0; 100 3000; 1 500; 0 1; 0 1];
%   [T,solver] = gsua_dataprep(odes, vars, [0 264], 'SIR', 'range', Ranges);
%
%%
%[text] ## Appendix C — User-defined functions
%[text] Store the model as a `.m` file. A **time-dependent** function has the form `sol = user_function(factors, domain)` (or `(factors, domain, opt)`), where `factors` is `Np x N` and `domain` is `[t0 tf]`, and returns a struct with `sol.x` (time, `1 x N`) and `sol.y` (outputs, `nOutputs x N`). The native output of MATLAB's `ode*`/`dde*` solvers already meets this; the toolbox accepts either `ode45` output orientation and normalizes it internally. A **time-free** function has the form `out = user_function(factors)` and returns a numeric array (vectorized) or a scalar; if its output is a vector, treat it as time-dependent and ignore the `domain` argument inside the function.
%[text] Prepare it with `[T,sol] = gsua_dataprep(func, Ranges, 'domain', domain)` (`func` is the file name as a char). `gsua_userdefined` does the same directly.
%[text] **Name–value arguments**
%[text] - **`'domain'`** (default `1`) — time span passed to the function when it takes 2+ inputs. \
%[text] - **`'rMethod'`** (default `'range'`) — read `Ranges` as `'range'`, `'std'`, `'percent'`, or `'normal'`. Fix a factor with a degenerate range. \
%[text] - **`'nominal'`** (default midpoint of `Ranges`) — explicit nominal factor family. \
%[text] - **`'output'`** (default all) — indices of the function's outputs to expose. \
%[text] - **`'names'`** (default auto `1..Np`) — factor names (cellstr, string array, or char), shown on plots; TeX and LaTeX are both rendered. \
%[text] - **`'out_names'`** (default `{'out'}`) — output names, shown on plots. \
%[text] - **`'vectorized'`** (default `false`) — set `true` if the function accepts a batch of factor sets at once (single-input functions only). \
%[text] - **`'opt'`** (default `struct('action',1)`) — anything passed through to a three-input function's third argument. \
%%
%[text]
%   func = 'user_dependent';
%   ranges = [0.1 0.3; 10 14];
%   [T,sol] = gsua_dataprep(func, ranges, 'domain', [0 25], 'names', {'\tau','\lambda'});
%

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
