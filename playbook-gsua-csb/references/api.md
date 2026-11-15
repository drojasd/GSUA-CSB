# API card: Functions

Generated mechanically from the help text of 43 `gsua_*.m` files (source hash `da4251f1123ece02`). It repeats what the source says; regenerate it when the toolbox changes instead of editing it by hand.

### `gsua_bootstrap`

`[Est,costs,out] = gsua_bootstrap(T,xdata,ydata,Est0,res,varargin)`  
Source: `gsua_bootstrap.m`

```text
GSUA_BOOTSTRAP Parametric-bootstrap ensemble of parameter estimates for identifiability analysis

  [Est,costs] = gsua_bootstrap(T,xdata,ydata,Est0,res)
  [Est,costs,out] = gsua_bootstrap(T,xdata,ydata,Est0,res,'B',200,'solver','fmincon', ...
          'N',20,'margin',0.3,'alpha',2,'noise','quasipoisson','cumulative',[false], ...
          'seed',1,'parallel',false,'opt',optimoptions('fmincon','Display','off'))

  Re-fitting the SAME data from many starting points explores the cost function (defends against
  local minima), but it collapses to the same optimum: it does NOT measure how much the DATA
  NOISE moves the estimate. In the discrete-SIR test case, multistart drives T.Range of S0 to a
  near-degenerate [2982.7, 2982.7] when the truth is 5000. A parametric bootstrap measures that
  missing quantity: it perturbs the data (not the start point) with the model's own observation
  noise and re-estimates, so the spread of the resulting ensemble is the parameter uncertainty.

  Method:
    1. Take the best pool member (by res) as theta*; simulate it to get the fitted curves f.
    2. Estimate a pooled observation-noise model from theta*'s residuals (GSUA_NOISEMODEL) -- the
       SAME dispersion estimate GSUA_NOISEFLOOR uses.
    3. Draw B synthetic datasets from f under that noise model (GSUA_NOISEDRAW): 'poisson',
       'quasipoisson' (default) or 'nbglobal', with cumulative-series support.
    4. Re-fit each synthetic dataset with GSUA_PE (same solver/margin/alpha, N starts per
       replica) and keep each replica's best estimate.
    5. Return the ensemble Est (Np x B) and its costs, ready to pass to GSUA_IA as T_est.

  Multistart vs bootstrap (see also the GSUA_IA help): the multistart ensemble answers "did the
  optimizer find the global optimum?"; the bootstrap ensemble answers "how much would the
  estimate change under a fresh draw of the data's noise?". Pass the bootstrap ensemble to
  GSUA_IA for honest parameter-uncertainty intervals and correlations.

  Inputs:
    T     <-- summary table from GSUA_PE (must carry the Margin/Alpha CustomProperties GSUA_PE
              sets when it scored with GSUA_COSTF, unless 'margin'/'alpha' are supplied).
    xdata <-- domain the pool was fit against
    ydata <-- inputs x len reference data the pool was fit against (may contain NaN)
    Est0  <-- Np x N pool of estimated parameter vectors (e.g. T.Estfmincon)
    res   <-- N x 1 (or 1 x N) fit cost per pool member (same scale GSUA_PE returns)

  Name-value arguments:
    'B'          <-- number of bootstrap datasets / refits. Default 200.
    'solver'     <-- GSUA_PE solver used for the refits. Default 'fmincon'.
    'N'          <-- GSUA_PE starts per replica (multistart within each refit). Default 20.
    'margin'     <-- raw margin for the refits (GSUA_PE's own convention). Default: recovered
                     from T.Properties.CustomProperties.Margin (stored already offset; converted
                     back to raw here).
    'alpha'      <-- cost exponent for the refits. Default: recovered from T's Alpha property.
    'noise'      <-- noise model for the synthetic datasets. A single name applies to EVERY row:
                     'poisson', 'quasipoisson' (default) or 'nbglobal' -- these are COUNT models
                     (they assume integer-count data). For a joint fit mixing count rows and
                     fraction/continuous rows, pass a 1 x inputs cellstr with a PER-ROW model,
                     adding 'gaussian' (additive noise, std from that row's residuals) or 'none'
                     (no perturbation): e.g. {'quasipoisson','none'}. This avoids drawing count
                     noise (and zeros) on a 0-1 fraction row, which otherwise zeroes that row's
                     cost scale and makes the refit fail. See GSUA_NOISEFLOOR.
    'cumulative' <-- 1 x inputs logical; a true row is drawn on the incident series and
                     re-accumulated. Default: all false.
    'seed'       <-- seed for the synthetic-data draws (and, for parity, the GSUA_PE starts via
                     the global stream). The caller's global RNG state is restored on return.
    'parallel'   <-- if true, run the B refits with parfor. Default false.
    'opt'        <-- optimoptions passed through to GSUA_PE. Default: a Display-off options set
                     for the chosen solver, so the refits stay quiet (passing Display-on options
                     slows a desktop run dramatically without changing the estimates).
    'timer'      <-- if true (default), print progress.

  Outputs:
    Est   <-- Np x B ensemble of best-per-replica estimates (pass directly to GSUA_IA)
    costs <-- 1 x B best cost per replica
    out   <-- struct with thetaStar, f, phi, k, noise, quasiFallback, exampleSynthetic
              (first synthetic draw), and the B/solver/N/margin/alpha actually used.

  Example:
    [T,res] = gsua_pe(T,xdata,ydata,'solver','fmincon','margin',0.3,'N',200);
    Eb = gsua_bootstrap(T,xdata,ydata,T.Estfmincon,res,'B',200,'solver','fmincon','N',20,'seed',1);
    [T,ci] = gsua_ia(T,Eb,false,false,true,true);   % honest uncertainty from the bootstrap

  See also GSUA_PE, GSUA_NOISEFLOOR, GSUA_NOISEMODEL, GSUA_NOISEDRAW, GSUA_IA.
```

### `gsua_clustersummary`

`S = gsua_clustersummary(T,clusterInfo)`  
Source: `gsua_clustersummary.m`

```text
GSUA_CLUSTERSUMMARY Per-cluster, per-parameter summary with an identifiability verdict

  S = gsua_clustersummary(T,clusterInfo)

  Builds a table with one row per (cluster, parameter) from the clusterInfo returned by GSUA_IA
  or GSUA_DIA. For each parameter, in each cluster, it reports where the group CENTER sits and how
  WIDE the plausible range is, and a suggested identifiability verdict -- the summary used in the
  epidemiology hands-on of the Semillero course.

  Two kinds of interval are reported, and they answer different questions:
    * MedianCI / MeanCI describe the CENTER of the group (how well pinned the central value is):
      the distribution-free confidence interval of the median (GSUA_MEDIANCI) and the Student-t
      confidence interval of the mean. They shrink as the cluster grows.
    * The 2.5-97.5 PERCENTILE interval [P2_5, P97_5] describes the SPREAD of plausible parameter
      values across the ensemble. Its width relative to the parameter's ORIGINAL range
      (PctWidthRel) is what the verdict is based on -- it does NOT shrink just because more
      estimates were collected.

  Verdict (on PctWidthRel): < 0.25 -> 'identifiable'; > 0.50 -> 'not identifiable'; in between
  -> 'dubious'.

  Inputs:
    T           <-- summary table (supplies the parameter names via its RowNames)
    clusterInfo <-- second output of GSUA_IA/GSUA_DIA (uses .NumClusters, .Est, .OrigRange)

  Output:
    S <-- table with columns Cluster, Parameter, n, Median, MedianCI_low, MedianCI_high, Mean,
          MeanCI_low, MeanCI_high, P2_5, P97_5, PctWidthRel, Verdict

  Example:
    [T,ci] = gsua_ia(T,Est,false,false,true,true);
    S = gsua_clustersummary(T,ci);
    disp(S)

  See also GSUA_IA, GSUA_DIA, GSUA_MEDIANCI, GSUA_BOOTSTRAP.
```

### `gsua_corr2`

`r = gsua_corr2(a,b)`  
Source: `gsua_corr2.m`

```text
GSUA_CORR2 2-D (Pearson) correlation coefficient without Image Processing Toolbox.

  r = gsua_corr2(a,b)

  Returns the correlation coefficient between two equally sized arrays a
  and b, computed over all their elements. This reproduces the Image
  Processing Toolbox function CORR2 exactly -- including returning NaN
  (0/0) when either input is constant -- so the toolbox's cost and
  sensitivity functions no longer depend on that product.

  Inputs:
    a,b <-- equal-size vectors or matrices to correlate

  Outputs:
    r <-- correlation coefficient in [-1,1], or NaN if an input is constant

  See also GSUA_CORR2OMITNAN.
```

### `gsua_corr2omitnan`

`r = gsua_corr2omitnan(a,b)`  
Source: `gsua_corr2omitnan.m`

```text
GSUA_CORR2OMITNAN 2-D correlation coefficient over jointly non-NaN entries

  r = gsua_corr2omitnan(a,b)

  NaN-tolerant replacement for corr2(a,b): computes the standard 2-D
  correlation coefficient using only the positions where both a and b
  are non-NaN. Identical to corr2(a,b) when neither input contains
  NaN. Returns 0 (no correlation credit, rather than erroring or
  propagating NaN) if fewer than 2 jointly-valid entries remain.

  Exists so GSUA_COSTF's correlation penalty term degrades gracefully
  when real-world reference data has missing (NaN) periods, instead of
  poisoning the whole objective with a single NaN.

  Inputs:
    a,b <-- equal-size vectors or matrices to correlate

  Outputs:
    r <-- 2-D correlation coefficient in [-1,1], or 0 if fewer than 2
          jointly non-NaN entries remain

  See also GSUA_COSTF, CORR2.
```

### `gsua_costcutoff`

`keep = gsua_costcutoff(cost,method,rtol,atol,gapRatio,minKeep)`  
Source: `gsua_costcutoff.m`

```text
GSUA_COSTCUTOFF Decide which repeated-estimation runs to keep based on fit cost

  keep = gsua_costcutoff(cost,method,rtol,atol,gapRatio,minKeep)

  Shared fit-quality filter used by GSUA_IA/GSUA_DIA to drop repeated
  estimation runs that converged to a much worse cost than the best
  one, before any identifiability statistic is computed from their
  parameter values. A run that converged to a bad local optimum
  contributes essentially arbitrary parameter values -- left in, it can
  make a genuinely well-identified parameter look poorly identified
  purely because an optimizer run failed, not because of real
  non-identifiability.

  Inputs:
    cost     <-- 1xN (or Nx1) fit cost, one per repeated estimation run
                 (e.g. the res output of gsua_pe)
    method   <-- 'rtol' -- keep cost <= min(cost)*(1+rtol)+atol.
                 Simple and predictable, but a fixed tolerance is a
                 somewhat arbitrary answer to a genuinely data-dependent
                 question.
                 'gap' -- keep runs up to the largest jump in the
                 *sorted* costs (in log space), the automated version
                 of reading a waterfall plot by eye: a cluster of
                 converged runs, a jump, then a scattered tail of
                 failures. Only cuts if that jump is at least gapRatio
                 times the *median* gap between sorted costs, so a
                 smooth continuum of similar costs is left untouched
                 rather than sliced arbitrarily.
    rtol     <-- relative tolerance for method='rtol'
    atol     <-- absolute tolerance added to the 'rtol' threshold, so
                 the filter doesn't become vacuously strict when the
                 best cost is at or near zero (e.g. a synthetic,
                 noiseless fit)
    gapRatio <-- salience ratio for method='gap'
    minKeep  <-- minimum number of runs to keep regardless of method or
                 how aggressive its cutoff would otherwise be -- a
                 degenerate 0- or 1-point result would otherwise
                 silently propagate into meaningless/NaN correlation,
                 confidence-interval, and clustering results downstream
                 instead of a clear signal that the filter was too
                 strict

  Output:
    keep <-- 1xN logical, true for runs to keep

  See also GSUA_IA, GSUA_DIA, GSUA_PE.
```

### `gsua_costf`

`cos=gsua_costf(inputs,regulator,len,ydata,yfunction,alpha)`  
Source: `gsua_costf.m`

```text
GSUA_COSTF Scalar correlation-penalized MSE cost function

  cos = gsua_costf(inputs,regulator,len,ydata,yfunction,alpha)

  Internal objective function used by GSUA_PE (for the 'ga',
  'particle', 'psearch', 'surrogate', 'annealing' and 'fmincon'
  solvers) to score a candidate parameter set. For each output row,
  the mean squared error between ydata and yfunction is normalized by
  regulator and penalized by a correlation factor (2-corr2), then
  raised to alpha; the per-output costs are averaged into a single
  scalar so any MATLAB scalar optimizer can be used.

  NaN-tolerant per output row: a NaN in ydata (e.g. a gap in one
  signal's real-world sampling, for a multi-output joint fit where
  different signals have independently missing periods) is excluded
  from that row's MSE and correlation terms rather than propagating
  through and poisoning the whole scalar cost. The per-row MSE is
  normalized by that row's own count of non-NaN samples, not the
  fixed len, so a gappy row isn't silently deflated relative to a
  fully-populated one. Identical to the pre-NaN-tolerant behavior when
  ydata has no NaN (the per-row non-NaN count then equals len exactly).

  Inputs:
    inputs    <-- number of output rows (signals) being compared
    regulator <-- 1 x inputs normalization vector, typically the MSE
                  between ydata and a margin-perturbed ydata (see
                  gsua_pe)
    len       <-- number of samples per output row (retained for
                  signature compatibility; the per-row MSE divisor is
                  now each row's own non-NaN count, see above)
    ydata     <-- inputs x len matrix of reference/experimental data,
                  may contain NaN for missing samples
    yfunction <-- inputs x len matrix of model output for the
                  candidate parameter set
    alpha     <-- exponent applied to each per-output cost term

  Outputs:
    cos <-- scalar cost (mean of the per-output penalized costs)

  Margin convention: the regulator GSUA_PE builds is mean((ydata - ydata*(m+1)).^2), which equals
  (m*rms(ydata)).^2 per row, where m is the RAW margin you passed GSUA_PE ('margin',m). So with a
  good-shape fit (correlation ~ 1) this cost falls below 1 exactly when the residual RMS is within
  a fraction m of the data RMS (e.g. m = 0.3 -> "within ~30%"). Because the regulator does not
  depend on the parameters, m only RESCALES the cost -- it does not move the optimum, but it sets
  what "cost < 1" and the GSUA_NOISEFLOOR threshold mean. See GSUA_PE's 'margin' help.

  See also GSUA_PE, GSUA_RCOSTF, GSUA_COSTFMULTI, GSUA_CORR2OMITNAN.
```

### `gsua_costfMulti`

`cos=gsua_costfMulti(ydata,yfunction,margin,parallel,alpha)`  
Source: `gsua_costfMulti.m`

```text
GSUA_COSTFMULTI Vectorized correlation-penalized MSE cost over many simulations

  cos = gsua_costfMulti(ydata,yfunction)
  cos = gsua_costfMulti(ydata,yfunction,margin,parallel,alpha)

  Batched version of GSUA_COSTF: scores every simulation ("run") in a
  3-D array of model outputs against a single reference yfunction in
  one call. Used by GSUA_SA (multi-objective sensitivity analysis) and
  GSUA_CSB to reduce a whole population of Monte-Carlo/design-matrix
  simulations to one scalar cost per run.

  Inputs:
    ydata     <-- reps x len x inputs array, one page of len-sample
                  output per input signal, one row per simulation run
                  (typically the output of gsua_pardeval)
    yfunction <-- inputs x len reference/nominal output (or a
                  1 x len x inputs array reshaped internally)
    margin    <-- (optional) relative perturbation used to build the
                  internal per-signal regulator (values <=1 are
                  treated as 1+abs(margin)). Default: 0.1
    parallel  <-- (optional, logical) evaluate runs with parfor.
                  Default: false
    alpha     <-- (optional) exponent applied to each per-signal cost
                  term. Default: 1

  Outputs:
    cos <-- 1 x reps row vector of scalar costs, one per run (Inf for
            runs where the cost computation failed)

  See also GSUA_COSTF, GSUA_SA, GSUA_CSB, GSUA_PARDEVAL.
```

### `gsua_covmetric`

`[cost_data,cost_band,P5,P50,P95] = gsua_covmetric(Y,ydata,varargin)`  
Source: `gsua_covmetric.m`

```text
GSUA_COVMETRIC Percentile-band coverage/tightness metrics for uncertainty-analysis results

  [cost_data,cost_band,P5,P50,P95] = gsua_covmetric(Y,ydata)
  [cost_data,cost_band,P5,P50,P95] = gsua_covmetric(Y,ydata,'margin',m,'alpha',a,'plevels',[lo hi])

  Reduces a Monte-Carlo uncertainty-analysis output Y (from GSUA_UA / GSUA_PARDEVAL) to a
  low/median/high percentile band per output signal and time step, then scores that band with
  two costs built on the same margin-normalized, correlation-penalized MSE used throughout the
  toolbox (see GSUA_COSTF): a value below 1 means "within the acceptable margin", a value at or
  above 1 means "not yet within tolerance". This gives the semi-automation identification
  routine a single, directly interpretable stopping signal instead of ad hoc envelope checks.

  cost_data (accuracy)  <-- distance from the band's median curve to ydata. High means the
      identified region's central tendency does not actually track the observed data.
  cost_band (precision) <-- distance between the low and high percentile curves themselves,
      independent of ydata. Low means repeated estimations/samples landed in a tight region,
      i.e. the identification has converged, REGARDLESS of how well that region fits the data.
      This distinction matters: a model can have a genuinely thin, converged identifiable
      region that still does not perfectly track noisy or structurally-imperfect data. That is
      a model-fit limitation, not an identifiability problem -- a semi-automation routine
      should recognize low cost_band as "stop refining, this has converged" even when
      cost_data never reaches 0.

  Inputs:
    Y       <-- reps x len (x nout) array of simulated trajectories, one row per Monte-Carlo
                run/estimation (output of gsua_ua/gsua_pardeval). Runs containing any
                non-finite value are excluded from the percentile band, with a warning.
    ydata   <-- nout x len (or 1 x len for a single output) reference/experimental data
    margin  <-- (optional) relative tolerance used to build the shared normalization
                regulator, same dual-mode convention as gsua_rcostf/gsua_costfMulti (values
                <1 are treated as a fraction and converted to 1+margin). Default: 0.1
    alpha   <-- (optional) exponent applied to each per-signal cost term, same convention as
                gsua_costf. Default: 2
    plevels <-- (optional) [low,high] percentile levels defining the band. Default: [5 95]

  Outputs:
    cost_data <-- scalar normalized cost of the band median vs ydata
    cost_band <-- scalar normalized cost of the high percentile curve vs the low percentile
                  curve (band tightness)
    P5,P50,P95 <-- nout x len low/median/high percentile curves actually used (named for the
                   default plevels but reflecting whatever plevels was passed; returned for
                   plotting/reuse, e.g. with gsua_plot)

  Example (semi-automation routine, Phase 1 reachability / Phase 4 stopping check):
    Y = gsua_ua(M,T,'xdata',xdata,'ynom',ydata,'parallel',true);
    [cost_data,cost_band] = gsua_covmetric(Y,ydata,'margin',0.1);
    reachable = cost_data < 1;               % Phase 1
    converged = cost_data < 1 && cost_band < 1;  % Phase 4

  See also GSUA_COSTF, GSUA_UA, GSUA_IA, GSUA_MEDIANCI, GSUA_PE.
```

### `gsua_csb`

`[ParU,New_range,J_test,Y_test,sup,condition]= gsua_csb(Par,N,varargin)`  
Source: `gsua_csb.m`

```text
Function for uncertainty-based confidence intervals

T= gsua_csb(T,N)
Parameters:
T <-- summary table from gsua_oatr function
N <-- Number of samples per cycle
Outputs:
T <-- summary table with new confidence intervals in T.Range
Additional features:
It is possible to request the following additional positional outputs
New_range <-- array with range modification record.
J_test <-- MSE output of last iteration
Y_test <-- output of last iteration
sup    <-- record of good scalars
Also, it is possible to apply the next paired features
'ynom',ynom <-- to especify a new output objective
'reps',k <--to perform k cycles. Default: 100.
'switch',s <-- to set another cut criteria (s). Default: 0.75.
'select',s <-- to set another select criteria (s). Default: 0.5.
'lim',l <-- to set another distance criteria (l). Default: 0.3.
'stop',p <-- to set the confidence of the interval in (p). Default: 0.95.
'parallel',false <-- to avoid the speed up of parallel computing toolbox.
[T,New_range,J_test,Y_test,sup]= gsua_csb(T,N,'lim',0.5,'select',0.4)

See also GSUA_OATR, GSUA_OATR2, GSUA_DMATRIX, GSUA_PARDEVAL.
```

### `gsua_dataprep`

`[T,out1]=gsua_dataprep(p1,p2,p3,p4,varargin)`  
Source: `gsua_dataprep.m`

```text
Consolidation of data and environment for (c) GSUA-RC Toolbox

Working with simulink models
[T,tout]=gsua_dataprep(model_name,Ranges,Parameter_names)
Parameters:
model_name      <-- file name of the simulink model
Ranges          <-- double array of Npx2
Parameter_names <-- parameter names of the model, as a cellstr, string array, or char
Outputs:
T    <-- table array with all required system information
tout <-- array of output time for the model

Working with symbolic toolbox, matlab
[T,solver]=gsua_dataprep(odes,vars,domain,model_name)
Parameters:
odes       <-- array with the system of differential equations
vars       <-- array of state variables for the system of equations
domain     <-- array with the domain of the system of equations
model_name <-- name for the output file of the system
Outputs:
T      <-- table array with all required system information
solver <-- handle function for ode solvers

Working with User-defined functions, matlab
[T,solver]=gsua_dataprep(Function,Ranges)
Parameters:
Function       <-- Is the name of the user-defined function (char)
Ranges       <-- double array of Npx2
Outputs:
T      <-- table array with all required system information
solver <-- handle function to solve the system

To explore additional configuration, open the user guide:
open gsua_userguide <--- User-defined functions userguide

See also GSUA_DPMAT, SENS_DATAPREP, GSUA_USERDEFINED, GSUA_DMATRIX.
```

### `gsua_depth`

`[depth,idx]=gsua_depth(x,parallel)`  
Source: `gsua_depth.m`

```text
GSUA_DEPTH Generalized band depth of a set of curves, ranked deepest-first

  [depth,idx] = gsua_depth(x)
  [depth,idx] = gsua_depth(x,parallel)

  Computes the generalized band depth of every curve in x relative to
  the full set x (Lopez-Pintado & Romo style band depth), a measure of
  how "central" each curve is among the others. Used by GSUA_SA in the
  multi-objective sensitivity analysis branch to automatically pick a
  representative reference output (the deepest curve) when no
  experimental/nominal output ('ynom') is supplied.

  Inputs:
    x        <-- n x d x inputs array: n curves, d sample points per
                 curve, inputs output signals (depth is averaged
                 across signals)
    parallel <-- (optional, logical) use parfor for the depth
                 computation. Default: false

  Outputs:
    depth <-- n x 1 vector of average band depth per curve
    idx   <-- indices of the curves sorted from deepest (most central,
               idx(1)) to shallowest

  See also GSUA_SA.
```

### `gsua_deval`

`y = gsua_deval(par,Table2,xdata)`  
Source: `gsua_deval.m`

```text
GSUA_DEVAL Evaluate the model behind a GSUA table for one parameter set

  y = gsua_deval(par,T,xdata)

  Internal dispatcher used throughout the toolbox (GSUA_PE, GSUA_OATR,
  GSUA_LIKELIHOOD, ...) to evaluate a single candidate parameter
  vector against whatever model kind T was built from (Simulink,
  symbolic ODE, or user-defined function) and return output
  interpolated onto xdata. Dispatch is driven by
  T.Properties.CustomProperties.Kind; on any evaluation error a matrix
  of Inf is returned so calling optimizers can reject the point
  instead of erroring out.

  Inputs:
    par   <-- row vector (or N x Np matrix for Kind==1, Simulink) of
              parameter values to evaluate
    T     <-- summary table from gsua_dataprep (or its
              CustomProperties struct directly)
    xdata <-- domain at which to evaluate/interpolate the model output

  Outputs:
    y <-- model output evaluated at xdata (Inf-filled if evaluation
          failed)

  See also GSUA_PARDEVAL, GSUA_EVAL, GSUA_DATAPREP.
```

### `gsua_dia`

`[T,clusterInfo] = gsua_dia(T,T_est,varargin)`  
Source: `gsua_dia.m`

```text
GSUA_DIA Practical identifiability diagnostics (headless, no plots)

  T = gsua_dia(T,T_est)
  T = gsua_dia(T,T_est,outlier)
  [T,clusterInfo] = gsua_dia(T,T_est,outlier,cluster)
  [T,clusterInfo] = gsua_dia(T,T_est,outlier,cluster,maxK,silThreshold)
  [T,clusterInfo] = gsua_dia(...,'cost',cost,'CostMethod','gap','IndexMethod','spread')

  Computes an updated confidence range and an identifiability index
  for each parameter from a set of repeated estimations, without
  producing any figures. This is the non-plotting counterpart of
  GSUA_IA: same range/index computation (parametric normal CI on the
  median) and same clusterInfo, useful for scripted/batch pipelines
  where the diagnostic plots of GSUA_IA are not needed.

  Multiple global minima (cluster = true):
  Same spectral-clustering detection as GSUA_IA, without plots: the
  normalized estimates are clustered for candidate counts 2:maxK using
  spectralcluster, the split with the best mean silhouette score
  (evalclusters) is kept only if that score is >= silThreshold,
  otherwise the estimates are treated as one basin. Requires Statistics
  and Machine Learning Toolbox. See GSUA_IA for the full explanation.

  Fit-quality filtering ('cost'):
  A repeated-estimation run that converged to a much worse cost than
  the best one contributes essentially arbitrary parameter values --
  left in, it can make a genuinely well-identified parameter look
  poorly identified purely because an optimizer run failed. Passing
  'cost' (e.g. the res output of gsua_pe, one entry per column of
  T_est) screens those runs out before anything else is computed. See
  GSUA_COSTCUTOFF for the two available methods ('rtol'/'gap') and the
  minKeep safety floor. The surviving columns are reported in
  clusterInfo.Kept.

  Clustering runs BEFORE outlier removal, on purpose: global
  Mahalanobis-distance outlier detection assumes one unimodal
  population, so running it first would treat a genuine second basin
  as "outliers" relative to the pooled mean/covariance across both
  basins and could delete it before clustering ever finds it. When a
  genuine multi-cluster split is found, T.Range/T.index and the
  correlation matrix are computed from the dominant (largest)
  cluster's points only -- pooling separated basins into one interval
  or correlation isn't a meaningful summary -- and outlier removal, if
  requested, is scoped to that single cluster, where the unimodal
  assumption behind Mahalanobis distance is actually valid.

  Note on outlier=true: the bundled DetectMultVarOutliers utility
  (Add Funcs/AntonSemechko-Multivariate-Outliers-.../) defaults to
  assuming up to roughly half the sample could be contaminated, which
  makes it considerably more aggressive than a plain Mahalanobis cutoff
  on data that is actually clean and tightly clustered -- confirmed to
  remove the large majority of a low-noise synthetic single cluster in
  testing, independent of this function's cluster/outlier ordering.
  This is a pre-existing characteristic of that third-party utility,
  not something introduced here; inspect T.Est after the call if
  outlier=true results look surprising.

  Inputs:
    T            <-- summary table from gsua_dataprep
    T_est        <-- Np x N matrix, one column per repeated parameter
                     estimation (e.g. Table.Estlsqc from gsua_pe)
    outlier      <-- (optional, logical) remove multivariate outliers
                     before computing the index. Default: false
    cluster      <-- (optional, logical) run spectral clustering to
                     detect multiple global minima. Default: false
    maxK         <-- (optional) largest cluster count to test. Default: 5
    silThreshold <-- (optional) minimum mean silhouette score required
                     to accept a multi-cluster split. Default: 0.5
    'cost'         <-- (optional, paired) 1xN fit cost, one per column
                       of T_est. Default: [] (no filtering)
    'CostMethod'   <-- (optional, paired) 'rtol' (default) or 'gap'
    'CostRtol'     <-- (optional, paired) Default: 0.1
    'CostAtol'     <-- (optional, paired) Default: 1e-8
    'CostGapRatio' <-- (optional, paired) Default: 3
    'MinKeep'      <-- (optional, paired) Default: 3
    'IndexMethod'  <-- (optional, paired) 'spread' (default) or 'ci'. See
                       GSUA_IA: 'spread' scores interval width by the 2.5-97.5
                       percentile dispersion of the estimates over the ORIGINAL
                       range (a flat parameter is not called identified); 'ci' is
                       the legacy median-CI width over the original range.

  Outputs:
    T           <-- same table with T.Range replaced by the new
                    confidence interval, T.Nominal set to the median
                    estimate, T.index set to the per-parameter
                    identifiability index (0 = well identified, 1 =
                    poorly identified; see IndexMethod), and T.Est set
                    to the estimates actually used. Reflects the dominant
                    cluster when NumClusters > 1.
    clusterInfo <-- struct with the same fields GSUA_IA returns:
                    NumClusters, Idx (1xN over the ORIGINAL columns, NaN
                    for cost-dropped), Kept (1xN logical over the original
                    columns), Centers, ClusterSize, Silhouette,
                    SilhouetteByK, SilhouetteK, Dominant, and the
                    per-cluster cells Range, Index, Est, Median. See GSUA_IA.

  See also GSUA_IA, GSUA_PE, GSUA_BOOTSTRAP, GSUA_MEDIANCI, GSUA_COSTCUTOFF, GSUA_CLUSTERSUMMARY.
```

### `gsua_dmatrix`

`[M,T] = gsua_dmatrix(Table,N,varargin)`  
Source: `gsua_dmatrix.m`

```text
Function for design of experiments (factor space sampling)

[M,T2]=gsua_dmatrix(T,N)
Parameters:
T <-- summary table from gsua_dataprep
N <-- number of samples
Outpus:
M  <-- design matrix of NxNp for later routines
T2 <-- summary table with fixed parameters actualized
Additional features:
You can choose a method for factor space sampling between uniform
distribution design and latin hypercube design (default). To switch
between methods use the paired feature 'Method' and 'Uniform',
'Sobol', or 'LatinHypercube'.
Also you can visualize the sampling result using the paired feature
'Show', 'on'
M=gsua_dmatrix(T,N,'Method','Sobol','Show','on')

MARGINAL VS JOINT SAMPLING
The methods above are MARGINAL: each parameter is drawn independently
inside its own range. For a model with parameter confounding (correlated,
marginally-non-identifiable parameters) that is wrong in a way that
matters -- independent draws leave the identified manifold, so the
resulting uncertainty band is both too wide and badly centered, with a
median trajectory resembling no actual good fit. GSUA_COVMETRIC reports
this condition as cost_band far above cost_data.

'Method','Joint' instead draws whole parameter VECTORS from an ensemble of
accepted estimates, so the correlation structure among parameters
survives:
M=gsua_dmatrix(Tia,N,'Method','Joint')

The ensemble is taken from Tia.Est (set by GSUA_IA, already filtered by
fit quality, restricted to the dominant cluster and outlier-cleaned), or
supplied explicitly with 'Pool' as an Np x nPool matrix.

Joint sampling also fixes a second, separate problem. GSUA_IA replaces
T.Range with a distribution-free confidence interval OF THE MEDIAN of the
pool (see GSUA_MEDIANCI) -- a statement about where the pool's centre
lies, which NARROWS as the pool grows. It is not the spread of parameter
values consistent with the data, so sampling it understates uncertainty.
Drawing pool vectors inherits the pool's own empirical spread instead.
For that reason 'Clip' defaults to no clipping: clipping joint draws back
to Tia.Range would reimpose exactly the narrowness this avoids. Pass the
ORIGINAL model's Range (physical/feasible bounds) if clipping is wanted.

'JointType' selects how the ensemble is turned into samples:
  'Bootstrap'       (default) resample whole pool columns with
                    replacement. Exact correlations, no distributional
                    assumption, cannot leave the manifold, and unaffected
                    by a rank-deficient covariance. Draws are limited to
                    the nPool distinct vectors actually in the pool.
  'SmoothBootstrap' as above plus a variance-corrected Gaussian kernel, so
                    draws are new points rather than repeats while staying
                    near the manifold. Bandwidth defaults to Silverman's
                    multivariate rule; override with 'Bandwidth'.
  'Gaussian'        ensemble mean and covariance, sampled by
                    eigendecomposition so a rank-deficient covariance
                    (nPool < number of free parameters) is handled by
                    sampling only the non-null eigendirections. Assumes
                    the manifold is linear -- on a CURVED ridge a Gaussian
                    places mass beside it, off the manifold.

CAVEATS
When GSUA_IA detects multiple global minima, Tia.Est holds the DOMINANT
cluster only (pooling separated basins into one interval is not
meaningful), so a joint band built from it is conditional on that basin.
Use clusterInfo.Idx from GSUA_IA to build per-basin bands instead.
A pool of fewer than 'MinPoolN' (default 5) vectors is refused: correlation
from very few points is not merely noisy but degenerate -- two points
always correlate at exactly +-1 -- so the result would be confidently
wrong rather than approximately right. Such a call warns and falls back to
marginal Latin hypercube sampling.

See also GSUA_DATAPREP, GSUA_SA, GSUA_UA, GSUA_PARDEVAL, GSUA_IA,
GSUA_MEDIANCI, GSUA_COVMETRIC.
```

### `gsua_dpmat`

`[Table,solver] = gsua_dpmat(odes,vars,domain,modelName,varargin)`  
Source: `gsua_dpmat.m`

```text
Consolidation of data and symbolic matlab environment for (c) GSUA-RC Toolbox

[T,solver]=gsua_dpmat(odes,vars,domain,model_name)
Parameters:
odes       <-- array with the system of differential equations
vars       <-- array of state variables for the system of equations
domain     <-- array with the domain of the system of equations
model_name <-- name for the output file of the system
Outputs:
T      <-- table array with all required system information
solver <-- handle function for ode solvers

To explore additional configuration, open the respective user guide:
open gsua_matuserguide <--- symbolic matlab userguide

See also GSUA_DATAPREP, GSUA_ODEFUN, GSUA_DEVAL.
```

### `gsua_eval`

`[y,xdata] = gsua_eval(par,Table,varargin)`  
Source: `gsua_eval.m`

```text
Function for few model evaluations

Usage:
  Y = gsua_eval(values,Table)
  Y = gsua_eval(values,Table,xdata)
  Y = gsua_eval(values,Table,xdata,ydata)
  Y = gsua_eval(values,Table,xdata,ydata,parallel)
  Y = gsua_eval(values,Table,xdata,ydata,parallel,show)
  Y = gsua_eval(values,Table,xdata,ydata,parallel,show,verbose)

Inputs:
  par      <-- array of Np x N (number of factors x number of simulations)
  Table    <-- summary table from gsua_dataprep function
  xdata    <-- (optional) domain to which model output is interpolated
  ydata    <-- (optional) experimental / previous model data matching xdata
  parallel <-- (optional, logical) use parallel evaluation (default: false)
  show     <-- (optional, logical) if true, calls gsua_plot; if false, no plot
  verbose  <-- (optional, logical) if true, prints "Progress: N %" during the
               batch loop (default: true, preserving prior behavior). Pass false
               for a quiet run -- e.g. inside a script whose output is captured.

Outputs:
  y     <-- array with model output
  xdata <-- domain used for evaluation/interpolation

See also GSUA_DEVAL, GSUA_PARDEVAL, GSUA_PLOT.
```

### `gsua_fullpars`

`[full, names] = gsua_fullpars(T, freeValues)`  
Source: `gsua_fullpars.m`

```text
GSUA_FULLPARS Reinsert fixed parameters into a free-parameter vector.

  [full, names] = gsua_fullpars(T, freeValues)

  A user-defined model pins a factor by giving it a degenerate range ([v v]); GSUA_DATAPREP then
  drops that row from T entirely, so T (and every estimate produced from it, e.g. T.Estfmincon)
  carries only the FREE factors. To evaluate a derived quantity, or to compare against the full
  ground truth, you need the complete ordered vector with the fixed values put back. This utility
  does exactly that: it places freeValues at the free positions and the pinned values at the
  fixed positions, in the original factor order, and returns the matching names.

  Inputs:
    T          <-- summary table from GSUA_DATAPREP for a USER-DEFINED model (Kind 5/6). The fixed
                   bookkeeping is read from its CustomProperties (FullNames/FixedMask/FixedValues),
                   with a fallback to the stored `creator` for tables built before v1.18.
    freeValues <-- the free-factor values, as a vector (row or column) of length equal to the
                   number of free factors (height(T)), or an Np_free x K matrix whose COLUMNS are
                   each a free vector (e.g. a whole estimate ensemble T.Estfmincon). Each column is
                   expanded independently.

  Outputs:
    full  <-- full ordered parameter vector (Ntot x 1), or Ntot x K if freeValues had K columns.
    names <-- Ntot x 1 cellstr of the full ordered factor names.

  Example:
    [T,~] = gsua_dataprep('myModel',[0 10; 5 5; 0 1]);   % factor 2 pinned at 5 -> dropped from T
    [T,res] = gsua_pe(T,xdata,ydata,'solver','fmincon');
    [theta, names] = gsua_fullpars(T, T.Estfmincon(:,1));  % best fit, with the fixed 5 restored

  See also GSUA_DATAPREP, GSUA_USERDEFINED, GSUA_DEVAL.
```

### `gsua_ia`

`[T,clusterInfo] = gsua_ia(T,T_est,varargin)`  
Source: `gsua_ia.m`

```text
GSUA_IA Practical identifiability analysis with diagnostic plots

  T = gsua_ia(T,T_est)
  T = gsua_ia(T,T_est,correction,outlier,show)
  [T,clusterInfo] = gsua_ia(T,T_est,correction,outlier,show,cluster)
  [T,clusterInfo] = gsua_ia(T,T_est,correction,outlier,show,cluster,maxK,silThreshold)
  [T,clusterInfo] = gsua_ia(...,'cost',cost,'CostMethod','gap')
  [T,clusterInfo] = gsua_ia(...,'IndexMethod','spread','truth',truth,'xdata',xdata,'ydata',ydata)

  Analyzes practical parameter identifiability from a set of repeated
  estimations: normalized CDFs and boxplots per parameter, a
  parameter-correlation heatmap, and (optionally) an identifiability
  index scatter/graph highlighting strongly correlated parameter
  pairs. Updates T.Range to a distribution-free confidence interval
  of the median (see GSUA_MEDIANCI) and reports an identifiability
  index per parameter combining interval width, mean correlation, and
  strong-correlation count. This is the plotting counterpart of
  GSUA_DIA.

  Multiple global minima (cluster = true):
  Repeated estimations from a non-convex problem can land in distinct
  basins of attraction ("multiple global minima") instead of scattering
  around one point. When cluster is true, gsua_ia runs spectral
  clustering (spectralcluster) on the normalized estimates for a range
  of candidate cluster counts 2:maxK, and picks the count with the best
  mean silhouette score (evalclusters). If the best mean silhouette is
  below silThreshold, the estimates are treated as a single basin
  (NumClusters = 1) -- a low silhouette score means the "best" split
  found is not actually well separated, so the sample is more likely
  unimodal than confidently multimodal. Requires Statistics and Machine
  Learning Toolbox.

  Multistart vs bootstrap ensemble:
  Feeding gsua_ia a MULTISTART ensemble (many fits of the SAME data)
  measures whether the optimizer found the global optimum, not parameter
  uncertainty -- it collapses to the best optimum. For honest
  uncertainty intervals and correlations, feed a PARAMETRIC-BOOTSTRAP
  ensemble instead (GSUA_BOOTSTRAP: refits of synthetic datasets drawn
  from the best fit's own observation noise), whose spread reflects how
  much the data noise moves each estimate.

  Fit-quality filtering ('cost'):
  A repeated-estimation run that converged to a much worse cost than
  the best one contributes essentially arbitrary parameter values --
  left in, it can make a genuinely well-identified parameter look
  poorly identified purely because an optimizer run failed. Passing
  'cost' (e.g. the res output of gsua_pe, one entry per column of
  T_est) screens those runs out before anything else is computed. See
  GSUA_COSTCUTOFF for the two available methods ('rtol'/'gap') and the
  minKeep safety floor. The surviving columns are reported in
  clusterInfo.Kept (see below).

  Clustering runs BEFORE outlier removal, on purpose: global
  Mahalanobis-distance outlier detection assumes one unimodal
  population, so running it first would treat a genuine second basin
  as "outliers" relative to the pooled mean/covariance across both
  basins and could delete it before clustering ever finds it. When a
  genuine multi-cluster split is found, every plot and statistic below
  (ECDFs, boxplots, the correlation heatmap, T.Range/T.index) reflects
  the dominant (largest) cluster's points only -- pooling separated
  basins into one interval or correlation isn't meaningful -- except
  the "Candidate global minima" plot and the per-cluster fields of
  clusterInfo, which cover every cluster. Outlier removal, if
  requested, is scoped to the dominant cluster only, where the
  unimodal assumption behind Mahalanobis distance is actually valid.

  Note on outlier=true: the bundled DetectMultVarOutliers utility
  (Add Funcs/AntonSemechko-Multivariate-Outliers-.../) defaults to
  assuming up to roughly half the sample could be contaminated, which
  makes it considerably more aggressive than a plain Mahalanobis cutoff
  on data that is actually clean and tightly clustered -- confirmed to
  remove the large majority of a low-noise synthetic single cluster in
  testing, independent of this function's cluster/outlier ordering.
  This is a pre-existing characteristic of that third-party utility,
  not something introduced here; inspect T.Est after the call if
  outlier=true results look surprising.

  Inputs:
    T            <-- summary table from gsua_dataprep
    T_est        <-- Np x N matrix, one column per repeated parameter
                     estimation (e.g. Table.Estlsqc from gsua_pe)
    correction   <-- (optional, logical) clip the new range to the
                     original T.Range when it would otherwise expand
                     beyond it. Default: false
    outlier      <-- (optional, logical) remove multivariate outliers
                     before computing statistics. Default: false
    show         <-- (optional, logical) plot the identifiability index
                     scatter, correlation graph, and (if cluster is
                     true) the clustering diagnostic plots. Default: false
    cluster      <-- (optional, logical) run spectral clustering to
                     detect multiple global minima among the repeated
                     estimations. Default: false
    maxK         <-- (optional) largest cluster count to test. Default: 5
    silThreshold <-- (optional) minimum mean silhouette score required
                     to accept a multi-cluster split over a single
                     basin. Default: 0.5
    'cost'         <-- (optional, paired) 1xN fit cost, one per column
                       of T_est. Default: [] (no filtering)
    'CostMethod'   <-- (optional, paired) 'rtol' (default) or 'gap'
    'CostRtol'     <-- (optional, paired) Default: 0.1
    'CostAtol'     <-- (optional, paired) Default: 1e-8
    'CostGapRatio' <-- (optional, paired) Default: 3
    'MinKeep'      <-- (optional, paired) Default: 3
    'IndexMethod'  <-- (optional, paired) 'spread' (default) or 'ci'. How the
                       identifiability index measures interval width. 'spread'
                       uses the 2.5-97.5 percentile dispersion of the estimates
                       relative to the ORIGINAL range, so a parameter scattered
                       across its whole prior is not reported as identified even
                       when its median is well pinned (e.g. an initial-condition
                       factor the data barely constrains). 'ci' is the legacy
                       behaviour: median-CI width over the original range, which
                       shrinks as the sample grows and can rate a flat parameter
                       as identified. The correlation terms are unchanged.
    'truth'        <-- (optional, paired) Np x 1 reference ("true") parameter
                       values (e.g. for synthetic teaching data). When given and
                       show is true, they are marked on the CDFs, the boxplots and
                       the per-cluster scatter. Default: [] (none)
    'xdata'        <-- (optional, paired) domain for the per-cluster curve plot
                       (#7); required together with 'ydata' to draw it.
    'ydata'        <-- (optional, paired) observed data for the per-cluster curve
                       plot: inputs x numel(xdata). With show and >1 cluster, each
                       cluster's center curve is overlaid on the data.

  Outputs:
    T           <-- same table with T.Range replaced by the new
                    confidence interval, T.Nominal set to the median
                    estimate, T.index set to the per-parameter
                    identifiability index (0 = well identified, 1 =
                    poorly identified; see IndexMethod), and T.Est set
                    to the estimates actually used (after cost filtering,
                    cluster restriction and outlier removal -- may be a
                    subset of the input T_est). Reflects the dominant
                    cluster when NumClusters > 1.
    clusterInfo <-- struct describing the clustering and per-cluster
                    statistics:
                      NumClusters   number of candidate global minima
                      Idx           1xN cluster label per ORIGINAL T_est
                                    column; NaN for columns dropped by
                                    the cost filter
                      Kept          1xN logical over the ORIGINAL T_est
                                    columns: which survived the cost
                                    filter (all true if no 'cost' given)
                      Centers       NumClusters x Np per-cluster medians
                                    (candidate global minima)
                      ClusterSize   estimation count per cluster
                      Silhouette    mean silhouette of the chosen split
                      SilhouetteByK mean silhouette for each tested K
                      SilhouetteK   the K values SilhouetteByK refers to
                      Dominant      index of the largest cluster (the one
                                    T.Range/T.index reflect)
                      Range         1xNumClusters cell; Range{c} is the
                                    Np x 2 median CI of cluster c
                      Index         1xNumClusters cell; Index{c} is the
                                    Np x 1 identifiability index of cluster c
                      Est           1xNumClusters cell; Est{c} is the
                                    Np x n_c estimates of cluster c
                      Median        1xNumClusters cell; Median{c} is the
                                    Np x 1 median of cluster c
                    Single-basin default (NumClusters = 1) when cluster
                    is false or no split beat silThreshold; the per-cluster
                    cells then hold one entry equal to the dominant stats.

  Example:
    [T,ci] = gsua_ia(T,Table.Estlsqc,false,false,true,true);
    if ci.NumClusters > 1
        disp('Candidate global minima (one row per cluster):')
        disp(ci.Centers)
        disp('Per-cluster identifiability index of cluster 2:')
        disp(ci.Index{2})
    end

  See also GSUA_DIA, GSUA_PE, GSUA_BOOTSTRAP, GSUA_MEDIANCI, GSUA_COSTCUTOFF, GSUA_CLUSTERSUMMARY.
```

### `gsua_iastats`

`[lb,up,ind,med,RHO,boxin,corrin,extrin] = gsua_iastats(statT_est,origRange,indexMethod,ciMethod,correction)`  
Source: `gsua_iastats.m`

```text
GSUA_IASTATS Shared range + identifiability-index computation for GSUA_IA/GSUA_DIA.

  [lb,up,ind,med,RHO,boxin,corrin,extrin] = ...
      gsua_iastats(statT_est,origRange,indexMethod,ciMethod,correction)

  Computes, for a single (presumed-unimodal) point set, the confidence interval of the median,
  the per-parameter identifiability index, and the correlation matrix. Factored out so GSUA_IA,
  GSUA_DIA and the per-cluster statistics all use one implementation.

  Inputs:
    statT_est   <-- Np x n point set (one column per estimation)
    origRange   <-- Np x 2 ORIGINAL parameter bounds (the index denominator)
    indexMethod <-- 'spread' (estimate 2.5-97.5 percentile dispersion over origRange) or 'ci'
                    (legacy: median-CI width over origRange)
    ciMethod    <-- 'orderstat' (distribution-free GSUA_MEDIANCI, GSUA_IA's) or 'normal'
                    (parametric med +/- 1.96*sqrt(pi/2)*std/sqrt(n), GSUA_DIA's)
    correction  <-- (optional, logical) clip the CI to origRange. Default: false

  Outputs:
    lb,up   <-- Np x 1 lower/upper confidence interval of the median
    ind     <-- 1 x Np identifiability index (0 = well identified, 1 = poorly identified)
    med     <-- Np x 1 median estimate
    RHO     <-- Np x Np parameter correlation matrix
    boxin   <-- Np x 1 interval-width term of the index
    corrin  <-- 1 x Np mean-correlation term
    extrin  <-- 1 x Np strong-correlation-count term

  See also GSUA_IA, GSUA_DIA, GSUA_MEDIANCI.
```

### `gsua_intrp`

`Sxint = gsua_intrp(sol,xint,idx,extra)`  
Source: `gsua_intrp.m`

```text
GSUA_INTRP Interpolate ODE solver output at requested time points

  Sxint = gsua_intrp(sol,xint,idx)
  Sxint = gsua_intrp(sol,xint,idx,extra)

  Internal helper (DEVAL-style) used to sample the time-domain output
  of an ODE solution structure (from ode45 or the bundled ODE4
  integrator) at the requested time points xint, selecting the state
  components given by idx. Used by GSUA_DEVAL and GSUA_PARDEVAL for
  Kind==4/5 (custom fixed-step or user-defined-with-domain models).

  Inputs:
    sol   <-- ODE solution struct with fields x (time) and y (states,
              possibly with a leading simulation-run dimension)
    xint  <-- time points to interpolate at
    idx   <-- indices of the state/output components to extract
    extra <-- (optional) 1 to use the struct/sol-array form (default),
              2 to use the simpler per-index interp1 form

  Outputs:
    Sxint <-- interpolated output at xint for the requested components

  See also GSUA_DEVAL, GSUA_PARDEVAL.
```

### `gsua_labels`

`[labels, interp] = gsua_labels(names)`  
Source: `gsua_labels.m`

```text
GSUA_LABELS Prepare parameter or output names for display and pick the text interpreter

  [labels, interp] = gsua_labels(names)

  Figures show names with MATLAB's default 'tex' interpreter, which understands TeX
  markup such as '\beta_h' or '\mu_{w,l}' but not LaTeX: a name written as '$\phi$'
  keeps its dollar signs, and commands such as \hat, \mathrm or \frac are printed
  literally (with a "valid interpreter syntax" warning). Names copied from a LaTeX
  document therefore rendered badly.

  If any name looks like LaTeX (contains '$' or a LaTeX-only command), the whole set
  switches to the 'latex' interpreter -- a legend or a tick-label set has a single
  interpreter, so the names must agree -- and each name is converted: already in
  $...$ is kept, TeX-style math (backslash, '_' or '^') is wrapped in $...$, and a
  plain name is left as text. Otherwise the names come back unchanged with 'tex',
  which is exactly the previous behavior.

  Inputs:
    names  <-- char, string array, or cell array of names

  Outputs:
    labels <-- the names ready to display: a char vector for a char input, a cell
               array of the same size otherwise
    interp <-- 'tex' or 'latex', to pass as the 'Interpreter' (or
               'TickLabelInterpreter') of whatever displays the labels

  See also GSUA_PLOT, GSUA_MCF, GSUA_IA.
```

### `gsua_likecost`

`cost=gsua_likecost(ydata,ynom,margin)`  
Source: `gsua_likecost.m`

```text
GSUA_LIKECOST Gaussian negative-log-likelihood cost relative to a margin baseline

  cost = gsua_likecost(ydata,ynom,margin)

  Standalone likelihood-style objective function: computes the
  Gaussian negative log-likelihood of each run in ydata against ynom
  (assuming a standard deviation proportional to ynom*margin), then
  subtracts the negative log-likelihood of the margin-perturbed
  baseline itself. The result is a relative cost that is zero when a
  run matches ynom as closely as the margin baseline does, and grows
  as the run deviates further. This mirrors the margin<0 objective
  built inline in GSUA_PE and can be used as a standalone objective
  function for custom optimization setups.

  Inputs:
    ydata  <-- reps x len x 1 array, one row per simulation run
    ynom   <-- 1 x len reference/nominal output
    margin <-- relative standard deviation used for both the
               likelihood weighting and the baseline (absolute value
               is used internally)

  Outputs:
    cost <-- 1 x reps row vector of relative negative log-likelihood
             costs, one per run

  See also GSUA_PE, GSUA_LIKELIHOOD.
```

### `gsua_likelihood`

`Tfinal = gsua_likelihood(T,xdata,ydata,alpha,step,margin,tolerance1,tolerance2,limit,reps,show,parallel,saver,pars,verbose)`  
Source: `gsua_likelihood.m`

```text
GSUA_LIKELIHOOD Profile-likelihood confidence intervals for each parameter

  Tfinal = gsua_likelihood(T,xdata,ydata,alpha,step,margin,tolerance1,tolerance2,limit,reps,show,parallel,saver,pars)

  Computes a profile-likelihood-based confidence interval for each
  parameter in T: starting from the nominal estimate, each parameter
  is stepped away (independently, up and down) while T is
  re-estimated (via GSUA_PE) with that parameter fixed at each trial
  value, until the resulting likelihood-ratio statistic crosses the
  chi-squared threshold implied by alpha. This is an alternative to
  the sampling-based ranges used by GSUA_OATR/GSUA_CSB, based directly
  on the likelihood-ratio test.

  Inputs:
    T          <-- summary table from gsua_dataprep (T.Nominal is the
                   reference parameter set)
    xdata      <-- domain at which the model is evaluated
    ydata      <-- experimental/reference output at xdata (if empty,
                   the nominal model output is used instead)
    alpha      <-- confidence level for the chi-squared threshold
                   (e.g. 0.95)
    step       <-- initial relative step size used to perturb each
                   parameter away from its nominal value
    margin     <-- RAW relative standard deviation used in the Gaussian
                   likelihood weighting (see gsua_likecost), e.g. 0.08 = 8%% noise. This is the
                   same raw convention as gsua_pe/gsua_costf (unified in v1.18). NOTE: before v1.18
                   this was offset by one (1.08 meant 8%%); pass the raw fraction now. A value >= 1
                   triggers a MarginConvention warning.
    tolerance1 <-- convergence tolerance on the likelihood-ratio
                   statistic (bisection stopping criterion)
    tolerance2 <-- convergence tolerance on the step size (relative to
                   the initial nominal-to-bound distance)
    limit      <-- maximum number of bisection iterations per
                   parameter bound
    reps       <-- number of estimation runs per GSUA_PE call while
                   profiling
    show       <-- (logical) plot each parameter's likelihood profile
                   and corresponding model fit as it is computed
    parallel   <-- (logical) let the inner GSUA_PE estimations use
                   parallel (fmincon UseParallel) optimization
    saver      <-- (logical) append each parameter's result to
                   'likelihoodCI.mat' as it is computed
    pars       <-- indices of the parameters to profile (empty = all)
    verbose    <-- (optional, logical) print "Progress: N %" timing during
                   profiling. Default: true (preserving prior behavior); pass
                   false for a quiet run, e.g. inside a captured script.

  Outputs:
    Tfinal <-- copy of T with Tfinal.Range replaced by the
               profile-likelihood interval bounds [lower,upper] per
               parameter

  See also GSUA_PE, GSUA_LIKECOST, GSUA_OATR, GSUA_CSB.
```

### `gsua_load`

`T=gsua_load(T)`  
Source: `gsua_load.m`

```text
GSUA_LOAD Rebuild a GSUA table's Solver handle after loading from disk

  T = gsua_load(T)

  gsua_save strips the CustomProperties.Solver function handle before
  saving a GSUA table (function handles to model/ODE definitions are
  not always safe or compact to serialize). gsua_load reverses that:
  it reads the table's stored 'creator' metadata (the original
  arguments used to build the table) and re-runs the matching
  gsua_dataprep/gsua_dpmat/gsua_userdefined call to reconstruct the
  Solver handle and any other derived CustomProperties.

  Inputs:
    T <-- GSUA table previously loaded from a .mat file saved with
          gsua_save (must still carry its CustomProperties.Kind and
          CustomProperties.creator fields)

  Outputs:
    T <-- same table with CustomProperties.Solver (and, for Simulink-
          or symbolic-derived tables, Domain) restored

  Example:
    load('portable.mat','T')
    T = gsua_load(T);

  See also GSUA_SAVE, GSUA_DATAPREP.
```

### `gsua_MCF`

`[Jsup,Jinf] = gsua_MCF(T,M,Y,y_exp,tIn,tOut)`  
Source: `gsua_MCF.m`

```text
GSUA_MCF Monte-Carlo filtering plots (behavioral vs non-behavioral parameter sets)

  [Jsup,Jinf] = gsua_MCF(T,M,Y,y_exp)
  [Jsup,Jinf] = gsua_MCF(T,M,Y,y_exp,tIn,tOut)

  Splits a design/sample matrix M into "low" and "high" subsets based
  on whether each run's scalar output is below/above the experimental
  reference y_exp (or, if tIn/tOut are given, based on the output
  value at a single time instant tIn), then plots the empirical CDF of
  each parameter for the full sample and both subsets. Called by
  GSUA_UA after an uncertainty analysis to visualize which parameter
  regions drive high/low model output.

  Inputs:
    T     <-- summary table from gsua_dataprep (used for
              CustomProperties.Fixed and parameter names)
    M     <-- N x Np design/sample matrix (from gsua_dmatrix)
    Y     <-- model output for each row of M (scalar per row, or
              N x Nd if dimension>1 and tIn/tOut are not given)
    y_exp <-- experimental/reference output to compare against
    tIn   <-- (optional) single time instant to filter on
    tOut  <-- (optional) time vector matching the columns of Y, used
              to locate tIn

  Outputs:
    Jsup <-- 1 x Np cell array, parameter values for runs above y_exp
    Jinf <-- 1 x Np cell array, parameter values for runs below y_exp

  See also GSUA_UA, GSUA_DMATRIX.
```

### `gsua_medianCI`

`[CI_lower, CI_upper] = gsua_medianCI(X, alpha)`  
Source: `gsua_medianCI.m`

```text
BINOMIALMEDIANCI  Compute a nonparametric (distribution-free) confidence
interval for the median of each row of X, using the order-statistics and
binomial approach.

  X is an M-by-N matrix of data. Each row is a separate variable
  (N samples per variable).

  ALPHA is the significance level (e.g., 0.05 for a 95% confidence interval).
  If ALPHA is omitted, it defaults to 0.05.

  [CI_LOWER, CI_UPPER] = BINOMIALMEDIANCI(...) returns two M-by-1 vectors
  containing the lower and upper confidence limits for each row's median.

  This method does NOT assume normality; it uses the distribution-free
  construction:
      CI = [ X_(r) ,  X_(s) ]
  where X_(k) denotes the k-th order statistic, and r and s are chosen
  by inverting a Binomial(n, 1/2) distribution to achieve coverage
  approximately 1 - ALPHA.

  See also GSUA_IA, GSUA_DIA.
```

### `gsua_mtest`

`[] = gsua_mtest(name,T,simvals,xdata,gname,pause)`  
Source: `gsua_mtest.m`

```text
GSUA_MTEST Animate model output as one parameter is swept, saved as a GIF

  gsua_mtest(name,T,simvals,xdata,gname,pause)

  Diagnostic/presentation utility: repeatedly overrides a single
  parameter (T{name,'gaming'}) with each value in simvals, evaluates
  the model against its nominal curve using GSUA_EVAL, captures each
  frame, and (when simvals has more than one row) assembles the frames
  into an animated GIF. Useful for visually exploring model
  sensitivity to one parameter at a time.

  Inputs:
    name    <-- row name (char) of the parameter to sweep
    T       <-- summary table from gsua_dataprep, with an
                Estlsqc/estimation column providing the nominal
                parameter set
    simvals <-- Nframes x Np matrix of parameter sets to evaluate, one
                row per animation frame
    xdata   <-- domain at which to evaluate the model
    gname   <-- base file name (without extension) for the output GIF
    pause   <-- delay (seconds) between GIF frames

  Outputs: none (writes '<gname>.gif' to the current folder when
  simvals has more than one row)

  See also GSUA_EVAL, GSUA_DEVAL.
```

### `gsua_noisedraw`

`yb = gsua_noisedraw(model,muClampedGrid,cumulative,ydata,phi,kGlobal,quasiFallback,B,sigmaRow)`  
Source: `gsua_noisedraw.m`

```text
GSUA_NOISEDRAW Draw B synthetic datasets under a count-noise model (optionally per output row).

  yb = gsua_noisedraw(model,muClampedGrid,cumulative,ydata,phi,kGlobal,quasiFallback,B)
  yb = gsua_noisedraw(model,...,B,sigmaRow)

  Shared synthetic-data generator used by GSUA_NOISEFLOOR (which scores the true model against
  the draws) and GSUA_BOOTSTRAP (which refits each draw). Factored out so both calibrations draw
  from exactly the same noise models rather than keeping two copies in sync.

  Inputs:
    model          <-- noise model. Either a single name applied to every row -- 'poisson',
                       'quasipoisson' or 'nbglobal' (count models; see GSUA_NOISEFLOOR) -- or a
                       1 x inputs cellstr giving a PER-ROW model, where each entry may also be
                       'gaussian' (additive N(0,sigmaRow(i)^2), for fraction/continuous rows) or
                       'none' (leave the row at its fitted mean, no perturbation). A per-row model
                       is what lets a joint count+fraction fit avoid drawing count noise (and
                       zeros) on a 0-1 fraction row.
    muClampedGrid  <-- inputs x len mean grid, clamped at 0 (from GSUA_NOISEMODEL)
    cumulative     <-- 1 x inputs logical; a true row is drawn on the incident series and
                       re-accumulated with cumsum
    ydata          <-- inputs x len reference data; its NaN mask is reapplied to every draw
    phi            <-- pooled Pearson dispersion (used by 'quasipoisson')
    kGlobal        <-- method-of-moments global NB dispersion (used by 'nbglobal'; Inf -> Poisson)
    quasiFallback  <-- true if quasi-Poisson must fall back to Poisson draws (phi<=1)
    B              <-- number of synthetic datasets to draw
    sigmaRow       <-- (optional) 1 x inputs additive-Gaussian std per row, used only by 'gaussian'
                       rows. Default: zeros (a 'gaussian' row with sigma 0 is then just its mean).

  Output:
    yb <-- inputs x len x B synthetic datasets

  See also GSUA_NOISEMODEL, GSUA_NOISEFLOOR, GSUA_BOOTSTRAP.
```

### `gsua_noisefloor`

`out = gsua_noisefloor(T,xdata,ydata,Est,res,varargin)`  
Source: `gsua_noisefloor.m`

```text
GSUA_NOISEFLOOR Noise-calibrated fit-acceptance threshold via parametric bootstrap

  out = gsua_noisefloor(T,xdata,ydata,Est,res)
  out = gsua_noisefloor(T,xdata,ydata,Est,res,'margin',m,'alpha',a,'B',2000, ...
          'noise','quasipoisson','quantile',0.95,'cumulative',[false true],'seed',1)

  Replaces the scale-dependent "lims = sum(res < 1.5*res(1))" fit-acceptance rule with a
  threshold calibrated against the data's own observation noise instead of against the best fit
  found. The 1.5x rule asks "how close is this to the best fit"; this asks "is this
  statistically distinguishable from a perfect model given how noisy the data is" -- the former
  gets narrower as the fit improves (exactly backwards), the latter does not.

  Method (parametric bootstrap of the noise floor):
    1. Take the best pool member (by res) as theta*, simulate it to get fitted curves f_i(t).
    2. Estimate a pooled observation-noise model from theta*'s own residuals.
    3. Generate B synthetic datasets from f under that noise model.
    4. Score the TRUE model f against each synthetic dataset with GSUA_COSTF, recomputing
       regulator from the synthetic data exactly as GSUA_PE does -- so the resulting cost
       distribution is directly comparable to res.
    5. That distribution is the cost a perfect model incurs from noise alone; its upper
       quantile (default 95%) is the largest cost still consistent with being correct.
    6. Accept every pool member with res < threshold.

  Positioning: this is a goodness-of-fit BAND-SIZING CALIBRATION, not a formal confidence
  region -- it does not replace profile likelihood (GSUA_LIKELIHOOD) or the Confidence
  Sub-contour Box (GSUA_CSB), it is complementary and cheap (no refitting). It is a PARAMETRIC
  bootstrap: the noise floor is conditional on theta* being approximately correct -- under model
  misspecification the floor is optimistic. The quantile is a user choice with real
  consequences; 95% is a convention, not a derivation. Report the threshold alongside the old
  1.5x rule so the change in band size is visible.

  Inputs:
    T     <-- summary table from GSUA_PE. Must carry the Margin/Alpha CustomProperties GSUA_PE
              sets automatically when it scored with GSUA_COSTF (i.e. was called with a
              'margin' name-value pair that is neither 0 nor -1, in GSUA_PE's raw convention --
              see 'margin' below), unless both are supplied explicitly.
    xdata <-- domain GSUA_PE was fit against
    ydata <-- inputs x len reference data GSUA_PE was fit against, same orientation as GSUA_PE.
              May contain NaN for missing samples.
    Est   <-- Np x N pool of estimated parameter vectors (e.g. T.Estfmincon)
    res   <-- N x 1 (or 1 x N) fit cost per pool member, same scale as GSUA_COSTF/GSUA_PE's res

  Name-value arguments:
    'margin'     <-- Override for the margin GSUA_PE was called with, in GSUA_PE's own raw
                     user-facing convention (i.e. pass the same value you gave GSUA_PE's
                     'margin', NOT the internally-offset one). The raw value is a relative
                     residual band (e.g. 0.3 ~ 30% of the data RMS); it is used internally as
                     m+1 and only rescales the cost without moving the optimum -- see GSUA_PE/
                     GSUA_COSTF. Default: recovered from T.Properties.CustomProperties.Margin
                     (already offset there).
    'alpha'      <-- Override for the cost exponent GSUA_PE was called with. Default: recovered
                     from T.Properties.CustomProperties.Alpha.
    'B'          <-- Bootstrap sample size. Default 2000.
    'noise'      <-- Noise model used for the primary out.threshold/out.cost: 'poisson',
                     'quasipoisson' (default, recommended -- see below), or 'nbglobal'.
    'quantile'   <-- Upper quantile of the bootstrap cost distribution defining the threshold.
                     Default 0.95.
    'cumulative' <-- 1 x inputs logical, one flag per fitted output row (default: all false).
                     For a row whose ydata is a CUMULATIVE series, noise is estimated and
                     generated on the INCIDENT (first-differenced) series and re-accumulated --
                     applying independent noise directly to a cumulative series would imply
                     independent errors on a running sum, which is wrong.
    'seed'       <-- Seed for the bootstrap draws, for reproducibility. Default: not set (draws
                     from the current global stream). The caller's global rng state is always
                     restored on return (or on error), regardless of whether 'seed' was given.
    'Show'       <-- {false,true}. If true, plot the sorted pool costs (log scale) with the
                     threshold and accepted/rejected runs marked, plus the bootstrap cost
                     distribution with the chosen quantile. Default false.

  Noise models (do not default to Poisson -- it is almost always far too tight for real
  surveillance-style data, since real observation processes are typically overdispersed
  relative to a Poisson count model):
    'poisson'      var = mu. Included for completeness/comparison, never the default.
    'quasipoisson' var = phi*mu, phi = pooled Pearson dispersion sum(r^2/mu)/(n-p). RECOMMENDED
                   DEFAULT. Implemented as NB2 with a per-point k = mu/(phi-1) (yields
                   var = mu + mu^2/k = phi*mu exactly). Falls back to Poisson draws with a
                   warning if phi <= 1 (data equi/under-dispersed relative to Poisson, for which
                   this reparameterization is degenerate).
    'nbglobal'     var = mu + mu^2/k, a single pooled k fitted by method of moments
                   (k = sum(mu^2)/(sum(r^2)-sum(mu))). Report as a bound only -- a single global
                   k is fitted almost entirely by the largest-mu points (e.g. epidemic peaks)
                   and can predict implausible spread elsewhere. Falls back to Poisson draws
                   with a warning if the method-of-moments estimate is degenerate (k<=0).
  out.byModel reports the threshold/distribution under all three regardless of 'noise', as a
  sanity check -- quasi-Poisson and global-NB agreeing closely is reassuring; a large
  disagreement is itself diagnostic.

  Outputs (out struct):
    out.threshold     <-- scalar cutoff on the GSUA_COSTF/res scale, for the chosen 'noise' model
    out.cost           <-- B x 1 bootstrap cost distribution, for the chosen 'noise' model
    out.accepted       <-- N x 1 logical, res < out.threshold -- replaces "res < 1.5*res(1)"
    out.phi, out.k     <-- pooled dispersion estimates (Pearson phi, method-of-moments global k)
    out.byModel        <-- struct with .poisson/.quasipoisson/.nbglobal, each a struct with
                            .threshold and .cost, for the sensitivity comparison above
    out.model          <-- which model out.threshold/out.cost came from
    out.quasiFallback  <-- true if quasi-Poisson degenerated to Poisson (phi<=1)
    out.exampleSynthetic <-- struct with .poisson/.quasipoisson/.nbglobal, each one example
                            inputs x len synthetic dataset (the first bootstrap draw) -- for
                            plotting a representative noise-floor sample against real data

  If the best pool member's cost exceeds the threshold, a warning fires: the model itself is
  rejected as a description of the data at the chosen quantile, which is diagnostically
  important and must not pass silently.

  Idiom replacing "lims = sum(res < 1.5*res(1)); Est = T.Estfmincon(:,1:lims);":
      out = gsua_noisefloor(T,xdata,ydata,T.Estfmincon,res);
      Tia = gsua_ia(T, T.Estfmincon(:,out.accepted), ...);

  See also GSUA_PE, GSUA_COSTF, GSUA_CORR2OMITNAN, GSUA_IA, GSUA_LIKELIHOOD, GSUA_CSB.
```

### `gsua_noisemodel`

`[f,muClampedGrid,phi,kGlobal,quasiFallback,nValid,pFree] = gsua_noisemodel(T,xdata,ydata,thetaStar,cumulative)`  
Source: `gsua_noisemodel.m`

```text
GSUA_NOISEMODEL Simulate theta* and estimate a pooled count-noise model from its residuals.

  [f,muClampedGrid,phi,kGlobal,quasiFallback,nValid,pFree] = ...
      gsua_noisemodel(T,xdata,ydata,thetaStar,cumulative)

  Shared dispersion-estimation step behind GSUA_NOISEFLOOR and GSUA_BOOTSTRAP: simulate the best
  fit, pool its residuals across every fitted output, and estimate the quasi-Poisson (Pearson
  phi) and global negative-binomial (method-of-moments k) dispersion. Deliberately emits NO
  warnings -- it returns the raw diagnostics (nValid, pFree, and the degeneracy flags implied by
  kGlobal/quasiFallback) so each caller can warn with its OWN message identifier
  (GSUA_NOISEFLOOR keeps its 'gsua_noisefloor:...' ids; GSUA_BOOTSTRAP uses its own).

  Inputs:
    T          <-- summary table (supplies T.Range for the free-parameter count)
    xdata      <-- domain the pool was fit against
    ydata      <-- inputs x len reference data (may contain NaN)
    thetaStar  <-- 1 x Np best-fit parameter vector
    cumulative <-- (optional) 1 x inputs logical; a true row is differenced before dispersion
                   estimation. Default: all false.

  Outputs:
    f             <-- inputs x len simulated curves at thetaStar
    muClampedGrid <-- inputs x len mean grid (incident for cumulative rows), clamped at 0
    phi           <-- pooled Pearson dispersion sum(r^2/mu)/(nValid-pFree)
    kGlobal       <-- method-of-moments global NB k (Inf if degenerate / underdispersed)
    quasiFallback <-- true if phi <= 1 (quasi-Poisson must fall back to Poisson)
    nValid        <-- number of valid (finite, mu>0) pooled residuals used
    pFree         <-- number of free parameters (for the dispersion degrees of freedom)

  See also GSUA_NOISEDRAW, GSUA_NOISEFLOOR, GSUA_BOOTSTRAP.
```

### `gsua_oatr`

`ParT=gsua_oatr(ParT,varargin)`  
Source: `gsua_oatr.m`

```text
Function for one-at-a-time range expansion

T=gsua_oatr(T)
Parameters:
T <-- summary table from gsua_dataprep function. The objective output is
given by T.Nominal.
Outputs:
T <-- summary table with new factor intervals in T.Range
Paired features:
'lim',l        <-- distance limit for range expansion. Default: 0.3.
'parallel',tf  <-- parallelize the range calculus. Default: false.
'correct',tf   <-- clip expansions back to the original range. Default: false.
'show',tf      <-- plot the iterations. Default: false.
'titerlimit',n <-- cap on the total number of iterations. Default: inf.
T=gsua_oatr(T,'lim',0.2,'parallel',true)

See also GSUA_OATR2, GSUA_CSB, GSUA_DEVAL.
```

### `gsua_oatr2`

`ParT=gsua_oatr2(ParT,varargin)`  
Source: `gsua_oatr2.m`

```text
Function for one-at-a-time range expansion (second-generation)

T=gsua_oatr2(T)
Parameters:
T <-- summary table from gsua_dataprep function. The objective output is
given by T.Nominal.
Outputs:
T <-- summary table with new factor intervals in T.Range
Paired features:
'lim',l       <-- distance limit for range expansion. Default: 0.3.
'parallel',tf <-- parallelize the range calculus. Default: false.
'correct',tf  <-- clip expansions back to the original range. Default: false.
'show',tf     <-- plot the iterations. Default: false.
T=gsua_oatr2(T,'lim',0.2,'parallel',true)

See also GSUA_OATR, GSUA_CSB, GSUA_DEVAL.
```

### `gsua_odefun`

`f = gsua_odefun(expr, vars, varargin)`  
Source: `gsua_odefun.m`

```text
ODEFUNCTION   Convert a system (vector or matrix) of symbolic algebraic
  expressions to a MATLAB function handle suitable for ODE45, ODE15s,
  and other ODE solvers.

  f = odeFunction(expr, vars, options)
  f = odeFunction(expr, vars, p1, p2, ..., options)

  expr is a vector or matrix of symbolic expressions in the variables
  specified by the vector vars, which can consist of symbolic variables,
  univariate symbolic functions, or univariate symbolic function calls.

  If expr contains symbolic parameters other than the variables given
  in vars, you must specify these parameters as additional arguments
  p1,p2,... .
  These parameters can be symbolic variables, univariate or multivariate
  symbolic functions, or univariate or multivariate symbolic function
  calls.

  You can also use the name-value pair arguments acceptable by
  MATLABFUNCTION, except that you cannot use 'vars' to control the
  order of the input variables. ODEFUNCTION passes these name-value
  pair arguments to an internal MATLABFUNCTION call.

  The output f is a function handle to a function f(t,Y,p1,p2,...).
  The scalar input parameter t represents the 'time' variable. The
  (univariate) symbolic functions or symbolic function calls listed
  in vars depend on this variable.
  The column vector Y represents the state variables (functions of
  'time') given by the input argument vars.
  The arguments p1,p2,... represent the symbolic input parameters
  (any symbolic variables that do not appear in vars).

  Given numerical values v1,v2,... for the symbolic parameters
  p1,p2,..., the reduced function handle
    F = @(t,Y) f(t,Y,v1,v2,...)
  represents the function F(t,Y) = expr. The reduced function handle
  F is suitable as the first argument to various numeric ODE solvers
  such as ODE45, ODE15s and so on.

  Example:
    Consider the following system of differential algebraic equations
    (DAEs):

    >> syms x1(t) x2(t) a b r(t);
       expr = [diff(x1(t),t) == a*x1(t) + b*x2(t)^2,...
               x1(t)^2 + x2(t)^2 == r(t)^2];
       vars = [x1(t),x2(t)]

    x1(t) and x2(t) are the state variables of the system. The system
    also involves the parameters a and b, and the parameter function
    r(t).
    Extract the mass matrix of this DAE system:

    >>  [m,f] = massMatrixForm(expr,vars)

    Convert the expressions m and f to function handles depending on
    the variables x1, x2 and the parameters a, b, and r(t):

    >> m = odeFunction(m,vars,a,b,r(t))
    >> f = odeFunction(f,vars,a,b,r(t))

    Before using the numerical MATLAB ODE solvers, you must assign
    numerical values to the parameters of the system.

    >> a = -0.6;
       b = -0.1;
       r = @(t) cos(t)/(1+t^2);

    Find the reduced function handles M and F.

    >> M = @(t,Y) m(t,Y,a,b,r(t));
    >> F = @(t,Y) f(t,Y,a,b,r(t));

    Find consistent initial conditions for this DAE system:

    >> t0 = 0;
       y0 = [-r(t0)*sin(0.1); r(t0)*cos(0.1)];
       yp0= [a*y0(1) + b*y0(2)^2; 1.234];

    Set the mass matrix and call ODE15S with these initial conditions.

    >> opt = odeset('mass',M,'InitialSlope',yp0);
    >> ode15s(F,[t0,1],y0,opt)

See also: MASSMATRIXFORM, ODE45, ODE15S.
```

### `gsua_pardeval`

`y = gsua_pardeval(M,Table2,xdata,parallel,init,verbose)`  
Source: `gsua_pardeval.m`

```text
GSUA_PARDEVAL Evaluate the model behind a GSUA table for a batch of parameter sets

  y = gsua_pardeval(M,T,xdata,parallel)
  y = gsua_pardeval(M,T,xdata,parallel,init)

  Batched counterpart of GSUA_DEVAL: evaluates the model behind T once
  per row of the design/sample matrix M (optionally with parfor), and
  assembles the results into an output array. Used internally by
  GSUA_SA, GSUA_UA and GSUA_CSB to run the many simulations required
  for sensitivity/uncertainty analysis.

  Inputs:
    M        <-- N x Np design matrix, one parameter set per row
                 (from gsua_dmatrix)
    T        <-- summary table from gsua_dataprep (or its
                 CustomProperties struct directly)
    xdata    <-- domain at which to evaluate/interpolate the model
                 output
    parallel <-- (logical) evaluate rows of M with parfor
    init     <-- (optional) [i0,Nsim] progress-tracking info passed to
                 gsua_timer-style reporting: i0 is the starting index
                 of this call within a larger run, Nsim is the total
                 number of simulations across that run

  Outputs:
    y <-- N x numel(xdata) (or N x numel(xdata) x nout for
          multi-output models) array of model outputs

  See also GSUA_DEVAL, GSUA_SA, GSUA_UA, GSUA_CSB, GSUA_DMATRIX.
```

### `gsua_pe`

`[Table,res,IP] = gsua_pe(Table,xdata,ydata,varargin)`  
Source: `gsua_pe.m`

```text
Parameter estimation function
[T,res] = gsua_pe(T,xdata,ydata,N)
Parameters:
T     <-- summary table from gsua_dataprep function
xdata <-- array of points where the model will be evaluated
ydata <-- array with expected model output
Outputs:
T   <-- summary table with parameter estimation results. Besides T.Est<solver> (Np x N, columns
        sorted best-first), it now carries T.IP<solver> (Np x N, the initial point actually used
        for each run, same column order) and the CustomProperty Retries<solver> (1 x N, how many
        silent restarts from a random start point each estimate cost), so the retry behaviour is
        traceable instead of silent.
res <-- cost functions for each estimation
Also, for the gradient solvers ('fmincon','lsqc','lsqn'), if you do not pass an 'opt' with a
TypicalX set, gsua_pe fills TypicalX with the nominal factor magnitudes (max(abs(T.Nominal),eps)),
which greatly improves finite-difference scaling when factors span many orders of magnitude.
Additional paired features:
'N',N     <-- number of parameter estimations
'Multistart',k <-- activate multistart feature for lsqcurvefit optimizer,
perform k parameter estimations in each cycle
'solver',{'lsqc','lsqn','ga','particle','psearch','surrogate','annealing'}
<-- allows to choose among several matlab optimizers. Default: lsqc
'opt',optimoptions(optimizer,...) <-- allows to configure the respective
matlab optimizer, optimizer must match with the real matlab optimizer
name.
'ipoint',point <-- allows to give to the optimizer the initial points for
estimations. point must have the same number of columns as factors to
estimate and the same number of rows as N.
'Show',{'off','on'}. If activated, the function plot optimizer results
Default: off.
'save',{false,true} <-- if true, writes Estimations.mat to the CURRENT folder after
each estimation (a crash-resume checkpoint), OVERWRITING any existing
Estimations.mat there. Default: false (it was true before v1.15).
'margin',m <-- relative tolerance band for the correlation-penalized cost (GSUA_COSTF), used by
the 'ga','particle','psearch','surrogate','annealing' and 'fmincon' solvers. The value you pass
is a RELATIVE band on the residual: with a good-shape fit (correlation ~ 1) the cost drops below
1 exactly when the residual RMS is within a fraction m of the data's own RMS -- so 'margin',0.3
means "accept, as cost < 1, fits whose residual is within ~30% of the data scale". Internally the
regulator is built with m+1 (so 0.3 is used as 1.3); you ALWAYS pass the raw band m, and
GSUA_NOISEFLOOR/GSUA_BOOTSTRAP take the same raw m. The regulator does not depend on the
parameters, so margin only RESCALES the cost: it does NOT move the optimum, only the
interpretation of cost thresholds (cost < 1, the noise floor, the fit-quality cutoff). m = 0
selects plain MSE; m < 0 selects the Gaussian negative-log-likelihood path (GSUA_LIKELIHOOD).
Default: 0.
'alpha',a <-- exponent applied to each per-output term in GSUA_COSTF. Default: 2.
[T,res] = gsua_pe(T,xdata,ydata,'Show','on','solver','particle','N',3)
[T,res] = gsua_pe(T,xdata,ydata,'solver','fmincon','margin',0.3)  % cost<1 ~ residual within 30%

See also GSUA_DMATRIX, GSUA_DEVAL, GSUA_COSTF, GSUA_IA, GSUA_DIA.
```

### `gsua_plot`

`[]=gsua_plot(plot_type,p1,p2,p3,p4,p5)`  
Source: `gsua_plot.m`

```text
Visualization function for GSUA toolbox

[]=gsua_plot(plot_type,p1,p2,p3,p4,p5)
Check gsua_userguide for function features (section 4.2)
open gsua_userguide

See also GSUA_EVAL, GSUA_UA, GSUA_SA, GSUA_MCF, GSUA_DMATRIX.
```

### `gsua_rcostf`

`cos=gsua_rcostf(ydata,yfunction,margin,alpha)`  
Source: `gsua_rcostf.m`

```text
GSUA_RCOSTF Scalar correlation-penalized MSE cost with a self-derived regulator

  cos = gsua_rcostf(ydata,yfunction)
  cos = gsua_rcostf(ydata,yfunction,margin,alpha)

  Standalone alternative to GSUA_COSTF for use as a scalar objective
  function with MATLAB optimizers (fmincon, ga, particleswarm, ...).
  Unlike GSUA_COSTF, the normalization ("regulator") is derived
  internally from ydata itself (the MSE between ydata and a
  margin-perturbed copy of ydata), so no external regulator needs to
  be precomputed by the caller.

  Inputs:
    ydata     <-- inputs x len matrix of reference/experimental data
    yfunction <-- inputs x len matrix of model output for the
                  candidate parameter set
    margin    <-- (optional) relative perturbation used to build the
                  internal regulator (values <1 are treated as
                  1+abs(margin)). Default: 1.1
    alpha     <-- (optional) exponent applied to each per-output cost
                  term. Default: 2

  Outputs:
    cos <-- scalar cost (mean of the per-output penalized costs)

  See also GSUA_COSTF, GSUA_COSTFMULTI, GSUA_PE.
```

### `gsua_ref`

`optim=gsua_ref(est,table,xdata,ydata,exe)`  
Source: `gsua_ref.m`

```text
GSUA_REF Iteratively re-estimate a parameter set while nudging ranges away from their bounds

  T = gsua_ref(est,T,xdata,ydata,exe)

  Repeatedly re-runs GSUA_PE with fmincon starting from est: whenever
  a (non-excluded) parameter's current estimate sits within 10% of
  its range bound, that bound is stretched (widened by 50% on the
  binding side) and the parameter set is re-estimated. This continues
  until every non-excluded parameter estimate settles comfortably
  inside its range. Useful as a quick range-adequacy check/refinement
  after an initial GSUA_PE call, when you suspect the parameter
  bounds are actively constraining the fit rather than the model
  itself.

  Inputs:
    est   <-- initial parameter estimate (row vector, one value per
              table row)
    table <-- summary table from gsua_dataprep (or a prior gsua_pe
              result); table.Range is expanded in place as needed
    xdata <-- domain at which the model is evaluated
    ydata <-- experimental/reference output at xdata
    exe   <-- indices of parameters to exclude from the bound check
              (kept fixed at their current range)

  Outputs:
    optim <-- updated summary table after convergence, with Range
              expanded as needed and a final GSUA_PE estimation column

  See also GSUA_PE, GSUA_OATR, GSUA_OATR2.
```

### `gsua_sa`

`[ParT,Y,J] = gsua_sa(M,ParT,varargin)`  
Source: `gsua_sa.m`

```text
Function for sensitivity analysis

T=gsua_sa(M,T)
Parameters:
M <-- design matrix from gsua_dmatrix function
T <-- summary table from gsua_dataprep function
Output:
T <-- summary table with sensitivity indices.
Additional features:
It is possible to request 2 extra positional outputs, in this order:
Y <-- model output matrix (same as gsua_ua output), 2nd output
J <-- vector of model outputs in scalar representation (MSE), 3rd output
It is possible to choose between 7 SA methods using the paired feature
'SensMethod',method. Where method is one of the following:
'Saltelli'    <-- requires N*(Np/2+1) simulations
'Jansen'      <-- requires N*(Np/2+1) simulations
'Xiao'        <-- newest method (Default), requires N*(Np/2+1)
                  simulations
'Sobol'       <-- classic method, requires N*(Np+1) simulations
'Shapley'     <-- Shapley effects (see SHAPLEY EFFECTS below); writes
                  T.Sh / T.Sh_vec instead of T.Si / T.STi
'brute-force' <-- explores all possible combinations, N + 2Np*N^2
                  simulations
'OAT'         <-- non-global method, requires N*(Np+1) simulations.
                  Varies one factor at a time over its N sampled values
                  with the others at T.Nominal; Si is each factor's share
                  of the summed OAT variance of the cost J and STi that
                  variance itself. No interactions, no Si_vec/STi_vec.

WHAT EACH METHOD ESTIMATES
'Sobol' and 'Jansen' estimate the variance-based (Sobol') first-order
and total indices of the output Y at each time (Si_vec, STi_vec).
'Xiao' (Xiao, Lu & Wang 2018, distance components / energy distance)
estimates STi = E|YA-YABi|^a / E|YA-YB|^a with a = 'xAlpha' (default 1).
With 'xAlpha',2 it coincides with the Jansen (variance-based) index; with
the default 1 it is a moment-independent, distribution-based index, which
is more robust but can rank the factors differently from Sobol'.
'Saltelli' scores the squared-deviation field (Y-ynom).^2 rather than Y
itself, i.e. it measures which factors drive the squared departure from
the nominal trajectory -- a different quantity from Sobol'/Jansen.
For every method the scalar Si/STi are indices of the cost
J = sum((Y-ynom).^2) (time and series summed).

INDICES OUTSIDE [0,1]
Sensitivity indices lie in [0,1] by definition, but their Monte Carlo
ESTIMATES do not: with a small N an estimated STi can exceed 1 and an Si
can be negative. That is sampling noise of the estimator, not an error,
and it is deliberately not clipped (clipping would hide that N is too
small). It is most frequent for 'Saltelli', occasional for 'Sobol' and
'Jansen', and rare for 'Xiao' with xAlpha=1. Increase N until
sum(Si)/sum(abs(Si)) is stable near 1 (a ratio well below 1 means
spurious negative Si, i.e. N is too small).

SHAPLEY EFFECTS ('SensMethod','Shapley')
Sh_i shares the output variance among the factors, splitting every
interaction fairly between the factors involved (Owen 2014; Song, Nelson
& Staum 2016). They always sum to 1 and, for independent factors,
Si <= Sh_i <= STi. Unlike Sobol' indices they remain meaningful when the
factors are DEPENDENT (correlated or structurally coupled, e.g.
R0 = beta*S0/gamma). Choose the estimator with 'ShapleyType':
  'Independent' (default) <-- factors sampled independently from T.Range
      (as gsua_dmatrix does). Nested Monte Carlo over random permutations
      (Song, Nelson & Staum 2016). Cost: N + m*(Np-1)*No*Ni simulations,
      with N = size(M,1) used for Var(Y). Options: 'ShapleyPerm' m
      (default N; when Np! <= m every permutation is used equally often),
      'ShapleyNo' outer draws (default 1), 'ShapleyNi' inner draws
      (default 3).
  'Given' <-- a JOINT sample whose dependence is kept: the rows of M
      (e.g. from gsua_dmatrix(Tia,N,'Method','Joint','JointType',
      'SmoothBootstrap'); the default 'Bootstrap' only repeats pool
      vectors) or 'Pool', an Np x nPool ensemble such as Tia.Est.
      Nearest-neighbour estimator
      of Broto, Bachoc & Depecker (2020); costs one simulation per
      distinct sample vector (duplicates are dropped). Option
      'ShapleyK' neighbours (default 3). It answers "how is the output
      variance shared WITHIN this joint distribution": with an accepted-
      estimate ensemble (Tia.Est) that is sensitivity inside the
      identified region, not over the prior ranges. Use a few hundred
      distinct vectors or more; small samples give noisy estimates.
Like any Monte Carlo index, an estimated Sh_i can come out slightly
negative at small N; the scalar T.Sh, built on the heavy-tailed cost J,
needs a larger N than T.Sh_vec to settle.
Results: T.Sh (Np x 1, indices of the cost J) and T.Sh_vec
(Np x Nd x nSeries, indices of Y at each time). Plot them like the other
indices, e.g. gsua_plot('Pie',T,T.Sh) or
gsua_plot('FractionalSensitivityArea',T,T.Sh_vec,xdata).

Also, it is possible to specify an output for indices estimation using
the paired feature 'ynom'. If you want to run the analysis with a
different resolution for the time-dependent sensitivity indices, you can
provide the paired feature 'xdata', [init:step:final].
Finally, you can avoid parallel speed up
for simulation process using the paired feature 'parallel', false.
[ParT,Y,J] = gsua_sa(M,T,'SensMethod','Xiao','ynom',ynom)
Note: ynom size must coincide with xdata size, where xdata is the next
xdata=linspace(Domain(1),Domain(2),Domain(2)-Domain(1)+1) and Domain is
the array of the model domain. If you are using the paired feature
'xdata', ignore the last advice.

See also GSUA_UA, GSUA_DMATRIX, GSUA_PARDEVAL, GSUA_COSTFMULTI, GSUA_DEPTH.
```

### `gsua_save`

`gsua_save(T,filename)`  
Source: `gsua_save.m`

```text
GSUA_SAVE Save a GSUA table to disk in a portable form

  gsua_save(T)
  gsua_save(T,filename)

  Removes the CustomProperties.Solver function handle from T before
  saving, since it can reference local functions, symbolic engine
  state, or Simulink model handles that do not serialize reliably
  across sessions or MATLAB versions. Use gsua_load after loading the
  file back in to reconstruct the Solver handle.

  Inputs:
    T        <-- GSUA table (from gsua_dataprep or a downstream step)
    filename <-- (optional) output .mat file name. Default: 'portable.mat'

  Example:
    gsua_save(T,'myModel.mat')

  See also GSUA_LOAD, GSUA_DATAPREP.
```

### `gsua_timer`

`t1 = gsua_timer(rem_compute,init,verbose)`  
Source: `gsua_timer.m`

```text
GSUA_TIMER Print elapsed/remaining-time progress for long simulation loops

  gsua_timer(rem_compute,init)
  t1 = gsua_timer(rem_compute,init)

  Internal progress-reporting helper used by GSUA_PE and
  GSUA_LIKELIHOOD. Call once with rem_compute=0 to start the shared
  timer (kept in a persistent variable), then repeatedly with rem_compute=1
  to print progress percentage, elapsed time, and estimated remaining
  time to the command window.

  Inputs:
    rem_compute <-- 0 to (re)start the timer, 1 to print a progress
                    update
    init        <-- [i,N] current iteration index and total number of
                    iterations

  Outputs:
    t1 <-- current time reading in seconds (only set when
           rem_compute is 1)

  See also GSUA_PE, GSUA_LIKELIHOOD.
```

### `gsua_ua`

`y = gsua_ua(M,T,varargin)`  
Source: `gsua_ua.m`

```text
Function for uncertainty analysis

Y=gsua_ua(M,T)
Parameters:
M <-- design matrix from gsua_dmatrix function
T <-- summary table from gsua_dataprep function
Output:
Y <-- result of Monte-Carlo simulation
Additional features:
The function automatically apply a Monte-Carlo filtering over results and
present it as figures (gsua_MCF function). It is possible to perform UA
over multiple model outputs. To obtain extrapolated results, use the
paired feature 'xdata'. To compare results with a specific output use
the paired feature 'ynom' (ynom length must coindice with xdata length).
To avoid parallel computing (no speed up), use the paired feature
'parallel',false. To silence the "Progress: N %" print during evaluation,
use the paired feature 'verbose',false (default true). To skip the
figures altogether -- e.g. in batch runs that call gsua_ua many times,
where one figure(1) plus one Monte-Carlo-filtering figure per output
would otherwise pile up -- use the paired feature 'plot',false
(default true). The returned Y is the same either way.
Y=gsua_ua(M,T,'xdata',xdata,'ynom',ynom,'parallel',false)
Y=gsua_ua(M,T,'plot',false,'verbose',false)   % quiet, no figures

See also GSUA_SA, GSUA_DMATRIX, GSUA_MCF, GSUA_PARDEVAL.
```

### `gsua_userdefined`

`[Table,solver] = gsua_userdefined(func,Range,varargin)`  
Source: `gsua_userdefined.m`

```text
GSUA_USERDEFINED Build a GSUA table around a user-defined MATLAB function

  [T,solver] = gsua_userdefined(func,Range)
  [T,solver] = gsua_userdefined(func,Range,'names',names,'domain',domain,...)

  Called by GSUA_DATAPREP when the model is neither a Simulink model
  nor a symbolic ODE system: func is any MATLAB function (by name)
  that maps a parameter vector to model output. This builds the
  summary table T (Range/Nominal plus the CustomProperties consumed
  by the rest of the toolbox: Solver, Kind, Domain, Fixed, ...) and a
  solver handle wrapping func with fixed parameters substituted in.

  Inputs:
    func  <-- name (char) of the user-defined function. Its number of
              declared inputs (1, 2, or 3) selects how it is called:
              solver(pars), solver(pars,domain), or
              solver(pars,domain,opt)
    Range <-- Np x 2 array of parameter ranges (or perturbation values,
              depending on 'rMethod')

  Paired features (name-value):
    'names',names       <-- parameter names (row names), as a cellstr, a string
                            array ["a" "b"], or a char. Omit to auto-number 1..Np
    'domain',domain     <-- domain passed through to func when it
                             accepts 2+ inputs. Default: 1
    'rMethod',method     <-- how to interpret Range: 'range' (Np x 2
                             bounds, default), 'percent', 'normal', or
                             'std'
    'nominal',nominal    <-- explicit nominal parameter values.
                             Default: midpoint of Range
    'output',output      <-- indices of func's output(s) to expose.
                             Default: all
    'out_names',names    <-- output names, as a cellstr, string array, or char.
                            Default: {'out'}
    'vectorized',tf      <-- true if func accepts a batch of parameter
                             sets at once (only used when func has a
                             single input). Default: false
    'opt',opt            <-- options struct passed through to func when
                             it accepts 3 inputs

  Outputs:
    Table  <-- summary table with all required system information
    solver <-- handle function wrapping func for direct evaluation

  See also GSUA_DATAPREP, GSUA_DPMAT, GSUA_DEVAL.
```
