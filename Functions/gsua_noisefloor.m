function out = gsua_noisefloor(T,xdata,ydata,Est,res,varargin)
%GSUA_NOISEFLOOR Noise-calibrated fit-acceptance threshold via parametric bootstrap
%
%   out = gsua_noisefloor(T,xdata,ydata,Est,res)
%   out = gsua_noisefloor(T,xdata,ydata,Est,res,'margin',m,'alpha',a,'B',2000, ...
%           'noise','quasipoisson','quantile',0.95,'cumulative',[false true],'seed',1)
%
%   Replaces the scale-dependent "lims = sum(res < 1.5*res(1))" fit-acceptance rule with a
%   threshold calibrated against the data's own observation noise instead of against the best fit
%   found. The 1.5x rule asks "how close is this to the best fit"; this asks "is this
%   statistically distinguishable from a perfect model given how noisy the data is" -- the former
%   gets narrower as the fit improves (exactly backwards), the latter does not.
%
%   Method (parametric bootstrap of the noise floor):
%     1. Take the best pool member (by res) as theta*, simulate it to get fitted curves f_i(t).
%     2. Estimate a pooled observation-noise model from theta*'s own residuals.
%     3. Generate B synthetic datasets from f under that noise model.
%     4. Score the TRUE model f against each synthetic dataset with GSUA_COSTF, recomputing
%        regulator from the synthetic data exactly as GSUA_PE does -- so the resulting cost
%        distribution is directly comparable to res.
%     5. That distribution is the cost a perfect model incurs from noise alone; its upper
%        quantile (default 95%) is the largest cost still consistent with being correct.
%     6. Accept every pool member with res < threshold.
%
%   Positioning: this is a goodness-of-fit BAND-SIZING CALIBRATION, not a formal confidence
%   region -- it does not replace profile likelihood (GSUA_LIKELIHOOD) or the Confidence
%   Sub-contour Box (GSUA_CSB), it is complementary and cheap (no refitting). It is a PARAMETRIC
%   bootstrap: the noise floor is conditional on theta* being approximately correct -- under model
%   misspecification the floor is optimistic. The quantile is a user choice with real
%   consequences; 95% is a convention, not a derivation. Report the threshold alongside the old
%   1.5x rule so the change in band size is visible.
%
%   Inputs:
%     T     <-- summary table from GSUA_PE. Must carry the Margin/Alpha CustomProperties GSUA_PE
%               sets automatically when it scored with GSUA_COSTF (i.e. was called with a
%               'margin' name-value pair that is neither 0 nor -1, in GSUA_PE's raw convention --
%               see 'margin' below), unless both are supplied explicitly.
%     xdata <-- domain GSUA_PE was fit against
%     ydata <-- inputs x len reference data GSUA_PE was fit against, same orientation as GSUA_PE.
%               May contain NaN for missing samples.
%     Est   <-- Np x N pool of estimated parameter vectors (e.g. T.Estfmincon)
%     res   <-- N x 1 (or 1 x N) fit cost per pool member, same scale as GSUA_COSTF/GSUA_PE's res
%
%   Name-value arguments:
%     'margin'     <-- Override for the margin GSUA_PE was called with, in GSUA_PE's own raw
%                      user-facing convention (i.e. pass the same value you gave GSUA_PE's
%                      'margin', NOT the internally-offset one). Default: recovered from
%                      T.Properties.CustomProperties.Margin (already offset there).
%     'alpha'      <-- Override for the cost exponent GSUA_PE was called with. Default: recovered
%                      from T.Properties.CustomProperties.Alpha.
%     'B'          <-- Bootstrap sample size. Default 2000.
%     'noise'      <-- Noise model used for the primary out.threshold/out.cost: 'poisson',
%                      'quasipoisson' (default, recommended -- see below), or 'nbglobal'.
%     'quantile'   <-- Upper quantile of the bootstrap cost distribution defining the threshold.
%                      Default 0.95.
%     'cumulative' <-- 1 x inputs logical, one flag per fitted output row (default: all false).
%                      For a row whose ydata is a CUMULATIVE series, noise is estimated and
%                      generated on the INCIDENT (first-differenced) series and re-accumulated --
%                      applying independent noise directly to a cumulative series would imply
%                      independent errors on a running sum, which is wrong.
%     'seed'       <-- Seed for the bootstrap draws, for reproducibility. Default: not set (draws
%                      from the current global stream). The caller's global rng state is always
%                      restored on return (or on error), regardless of whether 'seed' was given.
%
%   Noise models (do not default to Poisson -- it is almost always far too tight for real
%   surveillance-style data, since real observation processes are typically overdispersed
%   relative to a Poisson count model):
%     'poisson'      var = mu. Included for completeness/comparison, never the default.
%     'quasipoisson' var = phi*mu, phi = pooled Pearson dispersion sum(r^2/mu)/(n-p). RECOMMENDED
%                    DEFAULT. Implemented as NB2 with a per-point k = mu/(phi-1) (yields
%                    var = mu + mu^2/k = phi*mu exactly). Falls back to Poisson draws with a
%                    warning if phi <= 1 (data equi/under-dispersed relative to Poisson, for which
%                    this reparameterization is degenerate).
%     'nbglobal'     var = mu + mu^2/k, a single pooled k fitted by method of moments
%                    (k = sum(mu^2)/(sum(r^2)-sum(mu))). Report as a bound only -- a single global
%                    k is fitted almost entirely by the largest-mu points (e.g. epidemic peaks)
%                    and can predict implausible spread elsewhere. Falls back to Poisson draws
%                    with a warning if the method-of-moments estimate is degenerate (k<=0).
%   out.byModel reports the threshold/distribution under all three regardless of 'noise', as a
%   sanity check -- quasi-Poisson and global-NB agreeing closely is reassuring; a large
%   disagreement is itself diagnostic.
%
%   Outputs (out struct):
%     out.threshold     <-- scalar cutoff on the GSUA_COSTF/res scale, for the chosen 'noise' model
%     out.cost           <-- B x 1 bootstrap cost distribution, for the chosen 'noise' model
%     out.accepted       <-- N x 1 logical, res < out.threshold -- replaces "res < 1.5*res(1)"
%     out.phi, out.k     <-- pooled dispersion estimates (Pearson phi, method-of-moments global k)
%     out.byModel        <-- struct with .poisson/.quasipoisson/.nbglobal, each a struct with
%                             .threshold and .cost, for the sensitivity comparison above
%     out.model          <-- which model out.threshold/out.cost came from
%     out.quasiFallback  <-- true if quasi-Poisson degenerated to Poisson (phi<=1)
%     out.exampleSynthetic <-- struct with .poisson/.quasipoisson/.nbglobal, each one example
%                             inputs x len synthetic dataset (the first bootstrap draw) -- for
%                             plotting a representative noise-floor sample against real data
%
%   If the best pool member's cost exceeds the threshold, a warning fires: the model itself is
%   rejected as a description of the data at the chosen quantile, which is diagnostically
%   important and must not pass silently.
%
%   Idiom replacing "lims = sum(res < 1.5*res(1)); Est = T.Estfmincon(:,1:lims);":
%       out = gsua_noisefloor(T,xdata,ydata,T.Estfmincon,res);
%       Tia = gsua_ia(T, T.Estfmincon(:,out.accepted), ...);
%
%   See also GSUA_PE, GSUA_COSTF, GSUA_CORR2OMITNAN, GSUA_IA, GSUA_LIKELIHOOD, GSUA_CSB.

p = inputParser;
addRequired(p,'T',@istable);
addRequired(p,'xdata',@isnumeric);
addRequired(p,'ydata',@isnumeric);
addRequired(p,'Est',@isnumeric);
addRequired(p,'res',@isnumeric);
addParameter(p,'margin',[],@(x) isempty(x) || isnumeric(x));
addParameter(p,'alpha',[],@(x) isempty(x) || isnumeric(x));
addParameter(p,'B',2000,@(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p,'noise','quasipoisson',@(x) any(validatestring(x,{'poisson','quasipoisson','nbglobal'})));
addParameter(p,'quantile',0.95,@(x) isnumeric(x) && isscalar(x) && x>0 && x<1);
addParameter(p,'cumulative',[],@(x) isempty(x) || islogical(x));
addParameter(p,'seed',[],@(x) isempty(x) || isnumeric(x));
parse(p,T,xdata,ydata,Est,res,varargin{:});

xdata = p.Results.xdata;
ydata = p.Results.ydata;
Est = p.Results.Est;
res = p.Results.res(:);
B = round(p.Results.B);
noiseChoice = validatestring(p.Results.noise,{'poisson','quasipoisson','nbglobal'});
quantileLevel = p.Results.quantile;
seed = p.Results.seed;

inputs = size(ydata,1);
len = length(xdata);

cumulative = p.Results.cumulative;
if isempty(cumulative)
    cumulative = false(1,inputs);
end
if numel(cumulative) ~= inputs
    error('gsua_noisefloor:BadCumulative', ...
        'cumulative must have one entry per fitted output row (%d), got %d', inputs, numel(cumulative));
end

% margin/alpha: an explicit override is in GSUA_PE's own raw user-facing convention (the same
% number you would pass as GSUA_PE's 'margin'), and gets the identical +1 offset GSUA_PE applies
% before using it in the regulator formula; a table-recovered value is already offset. Both paths
% converge on marginInternal before it is ever used, so there is no way to apply the offset twice
% or not at all depending on which source was used.
if isempty(p.Results.margin)
    try
        marginInternal = T.Properties.CustomProperties.Margin;
    catch
        error('gsua_noisefloor:NoMarginAlpha', ...
            ['T has no Margin CustomProperty (only set by GSUA_PE when it scored with GSUA_COSTF) ' ...
             'and no ''margin'' override was given.']);
    end
else
    marginInternal = p.Results.margin + 1;
end
if isempty(p.Results.alpha)
    try
        costAlpha = T.Properties.CustomProperties.Alpha;
    catch
        error('gsua_noisefloor:NoMarginAlpha', ...
            ['T has no Alpha CustomProperty (only set by GSUA_PE when it scored with GSUA_COSTF) ' ...
             'and no ''alpha'' override was given.']);
    end
else
    costAlpha = p.Results.alpha;
end
if marginInternal == 1 || marginInternal <= 0
    error('gsua_noisefloor:InvalidMargin', ...
        ['margin (internal value %.6g) must not be 1 (degenerates the regulator) and must be ' ...
         'positive -- the noise floor is only meaningful for a T whose GSUA_PE run scored with ' ...
         'GSUA_COSTF.'], marginInternal);
end

if ~isempty(seed)
    priorState = rng;
    cleanupObj = onCleanup(@() rng(priorState));
    rng(seed);
end

% theta*: best pool member by res, not assumed to be column 1 -- robust to caller sort order.
[bestRes,ibest] = min(res);
thetaStar = Est(:,ibest).';
f = gsua_deval(thetaStar,T,xdata);   % inputs x len, assumed NaN-free (a clean simulation)

pFree = sum(T.Range(:,1) ~= T.Range(:,2));

% Per-output incident/cumulative transform, built once and reused for both dispersion estimation
% and synthetic generation. For a cumulative row, diff([0,row]) lets NaN propagate/smear naturally
% across any real gap in ydata -- mathematically correct (a cumulative gap genuinely does not tell
% you the incident breakdown across it), and only ever feeds a pooled scalar (phi/k), not a
% positional array, so no gap-bridging logic is needed.
muGrid = zeros(inputs,len);
rGrid = nan(inputs,len);
for i = 1:inputs
    if cumulative(i)
        fInc = diff([0, f(i,:)]);
        yInc = diff([0, ydata(i,:)]);
        muGrid(i,:) = fInc;
        rGrid(i,:) = yInc - fInc;
    else
        muGrid(i,:) = f(i,:);
        rGrid(i,:) = ydata(i,:) - f(i,:);
    end
end
muClampedGrid = max(muGrid,0);

% Pooled dispersion across every fitted output (not per-row): a single phi/k, matching how this
% calibration is reported and used downstream.
valid = isfinite(rGrid) & isfinite(muGrid) & muGrid > 0;
R = rGrid(valid);
MU = muGrid(valid);
n = numel(R);
if n - pFree <= 0
    warning('gsua_noisefloor:LowDF', ...
        'Only %d valid (non-NaN, mu>0) pooled residuals for %d free parameters; dispersion estimates may be unreliable.', ...
        n, pFree);
end
denom = max(n - pFree,1);
phi = sum(R.^2 ./ MU) / denom;

kDenom = sum(R.^2) - sum(MU);
if kDenom <= 0
    kGlobal = Inf;
    warning('gsua_noisefloor:UnderdispersedNB', ...
        'Method-of-moments global-NB k is degenerate (data at or below Poisson dispersion); nbglobal falls back to Poisson.');
else
    kGlobal = sum(MU.^2) / kDenom;
end

quasiFallback = phi <= 1;
if quasiFallback
    warning('gsua_noisefloor:QuasiPoissonDegenerate', ...
        'Pooled Pearson dispersion phi=%.6g <= 1; quasipoisson falls back to Poisson draws.', phi);
end

models = {'poisson','quasipoisson','nbglobal'};
byModel = struct();
exampleSynthetic = struct();
for m = 1:numel(models)
    model = models{m};
    yb = generate_bootstrap(model,muClampedGrid,cumulative,ydata,phi,kGlobal,quasiFallback,B,inputs,len);
    exampleSynthetic.(model) = yb(:,:,1);
    cost = zeros(B,1);
    for b = 1:B
        ybSlice = yb(:,:,b);
        regulator_b = sum((ybSlice - ybSlice*marginInternal).^2,2,'omitnan') ./ sum(~isnan(ybSlice),2);
        cost(b) = gsua_costf(inputs,regulator_b,len,ybSlice,f,costAlpha);
    end
    byModel.(model) = struct('threshold',quantile(cost,quantileLevel),'cost',cost);
end

out = struct();
out.threshold = byModel.(noiseChoice).threshold;
out.cost = byModel.(noiseChoice).cost;
out.accepted = res < out.threshold;
out.phi = phi;
out.k = kGlobal;
out.byModel = byModel;
out.model = noiseChoice;
out.quasiFallback = quasiFallback;
% One example synthetic dataset per model (inputs x len, first bootstrap draw) -- lets a caller
% plot a representative noise-floor sample against their real data/pool, and gives a concrete,
% inspectable artifact for verifying e.g. that a cumulative row's NaN mask matches ydata's own.
out.exampleSynthetic = exampleSynthetic;

if bestRes > out.threshold
    warning('gsua_noisefloor:BestFitExceedsNoiseFloor', ...
        ['Best pool cost (%.6g) exceeds the noise-floor threshold (%.6g) at the %.0f%% level: ' ...
         'the model is rejected as a description of the data.'], ...
        bestRes, out.threshold, quantileLevel*100);
end

end


function yb = generate_bootstrap(model,muClampedGrid,cumulative,ydata,phi,kGlobal,quasiFallback,B,inputs,len)
%GENERATE_BOOTSTRAP Draw B synthetic datasets (inputs x len x B) under the given noise model.
muRep = repmat(muClampedGrid,[1 1 B]);
switch model
    case 'poisson'
        S = poissrnd(muRep);
    case 'quasipoisson'
        if quasiFallback
            S = poissrnd(muRep);
        else
            kPoint = max(muClampedGrid ./ (phi-1), realmin);
            Rn = repmat(kPoint,[1 1 B]);
            Pn = Rn ./ (Rn + muRep);
            S = nbinrnd(Rn,Pn);
        end
    case 'nbglobal'
        if isinf(kGlobal)
            S = poissrnd(muRep);
        else
            Rn = max(kGlobal,realmin) * ones(size(muRep));
            Pn = Rn ./ (Rn + muRep);
            S = nbinrnd(Rn,Pn);
        end
end
S(isnan(S)) = 0;

yb = zeros(inputs,len,B);
for i = 1:inputs
    if cumulative(i)
        yb(i,:,:) = cumsum(S(i,:,:),2);
    else
        yb(i,:,:) = S(i,:,:);
    end
    yb(i,isnan(ydata(i,:)),:) = NaN;
end
end
