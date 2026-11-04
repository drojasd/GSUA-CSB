function [Est,costs,out] = gsua_bootstrap(T,xdata,ydata,Est0,res,varargin)
%GSUA_BOOTSTRAP Parametric-bootstrap ensemble of parameter estimates for identifiability analysis
%
%   [Est,costs] = gsua_bootstrap(T,xdata,ydata,Est0,res)
%   [Est,costs,out] = gsua_bootstrap(T,xdata,ydata,Est0,res,'B',200,'solver','fmincon', ...
%           'N',20,'margin',0.3,'alpha',2,'noise','quasipoisson','cumulative',[false], ...
%           'seed',1,'parallel',false,'opt',optimoptions('fmincon','Display','off'))
%
%   Re-fitting the SAME data from many starting points explores the cost function (defends against
%   local minima), but it collapses to the same optimum: it does NOT measure how much the DATA
%   NOISE moves the estimate. In the discrete-SIR test case, multistart drives T.Range of S0 to a
%   near-degenerate [2982.7, 2982.7] when the truth is 5000. A parametric bootstrap measures that
%   missing quantity: it perturbs the data (not the start point) with the model's own observation
%   noise and re-estimates, so the spread of the resulting ensemble is the parameter uncertainty.
%
%   Method:
%     1. Take the best pool member (by res) as theta*; simulate it to get the fitted curves f.
%     2. Estimate a pooled observation-noise model from theta*'s residuals (GSUA_NOISEMODEL) -- the
%        SAME dispersion estimate GSUA_NOISEFLOOR uses.
%     3. Draw B synthetic datasets from f under that noise model (GSUA_NOISEDRAW): 'poisson',
%        'quasipoisson' (default) or 'nbglobal', with cumulative-series support.
%     4. Re-fit each synthetic dataset with GSUA_PE (same solver/margin/alpha, N starts per
%        replica) and keep each replica's best estimate.
%     5. Return the ensemble Est (Np x B) and its costs, ready to pass to GSUA_IA as T_est.
%
%   Multistart vs bootstrap (see also the GSUA_IA help): the multistart ensemble answers "did the
%   optimizer find the global optimum?"; the bootstrap ensemble answers "how much would the
%   estimate change under a fresh draw of the data's noise?". Pass the bootstrap ensemble to
%   GSUA_IA for honest parameter-uncertainty intervals and correlations.
%
%   Inputs:
%     T     <-- summary table from GSUA_PE (must carry the Margin/Alpha CustomProperties GSUA_PE
%               sets when it scored with GSUA_COSTF, unless 'margin'/'alpha' are supplied).
%     xdata <-- domain the pool was fit against
%     ydata <-- inputs x len reference data the pool was fit against (may contain NaN)
%     Est0  <-- Np x N pool of estimated parameter vectors (e.g. T.Estfmincon)
%     res   <-- N x 1 (or 1 x N) fit cost per pool member (same scale GSUA_PE returns)
%
%   Name-value arguments:
%     'B'          <-- number of bootstrap datasets / refits. Default 200.
%     'solver'     <-- GSUA_PE solver used for the refits. Default 'fmincon'.
%     'N'          <-- GSUA_PE starts per replica (multistart within each refit). Default 20.
%     'margin'     <-- raw margin for the refits (GSUA_PE's own convention). Default: recovered
%                      from T.Properties.CustomProperties.Margin (stored already offset; converted
%                      back to raw here).
%     'alpha'      <-- cost exponent for the refits. Default: recovered from T's Alpha property.
%     'noise'      <-- 'poisson', 'quasipoisson' (default) or 'nbglobal'. See GSUA_NOISEFLOOR.
%     'cumulative' <-- 1 x inputs logical; a true row is drawn on the incident series and
%                      re-accumulated. Default: all false.
%     'seed'       <-- seed for the synthetic-data draws (and, for parity, the GSUA_PE starts via
%                      the global stream). The caller's global RNG state is restored on return.
%     'parallel'   <-- if true, run the B refits with parfor. Default false.
%     'opt'        <-- optimoptions passed through to GSUA_PE. Default: a Display-off options set
%                      for the chosen solver, so the refits stay quiet (passing Display-on options
%                      slows a desktop run dramatically without changing the estimates).
%     'timer'      <-- if true (default), print progress.
%
%   Outputs:
%     Est   <-- Np x B ensemble of best-per-replica estimates (pass directly to GSUA_IA)
%     costs <-- 1 x B best cost per replica
%     out   <-- struct with thetaStar, f, phi, k, noise, quasiFallback, exampleSynthetic
%               (first synthetic draw), and the B/solver/N/margin/alpha actually used.
%
%   Example:
%     [T,res] = gsua_pe(T,xdata,ydata,'solver','fmincon','margin',0.3,'N',200);
%     Eb = gsua_bootstrap(T,xdata,ydata,T.Estfmincon,res,'B',200,'solver','fmincon','N',20,'seed',1);
%     [T,ci] = gsua_ia(T,Eb,false,false,true,true);   % honest uncertainty from the bootstrap
%
%   See also GSUA_PE, GSUA_NOISEFLOOR, GSUA_NOISEMODEL, GSUA_NOISEDRAW, GSUA_IA.

p = inputParser;
addRequired(p,'T',@istable);
addRequired(p,'xdata',@isnumeric);
addRequired(p,'ydata',@isnumeric);
addRequired(p,'Est0',@isnumeric);
addRequired(p,'res',@isnumeric);
addParameter(p,'B',200,@(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p,'solver','fmincon',@(x) any(validatestring(x,{'lsqc','lsqn','ga','particle','psearch','surrogate','annealing','fmincon'})));
addParameter(p,'N',20,@(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p,'margin',[],@(x) isempty(x) || isnumeric(x));
addParameter(p,'alpha',[],@(x) isempty(x) || isnumeric(x));
addParameter(p,'noise','quasipoisson',@(x) any(validatestring(x,{'poisson','quasipoisson','nbglobal'})));
addParameter(p,'cumulative',[],@(x) isempty(x) || islogical(x));
addParameter(p,'seed',[],@(x) isempty(x) || isnumeric(x));
addParameter(p,'parallel',false,@islogical);
addParameter(p,'opt',[]);
addParameter(p,'timer',true,@islogical);
parse(p,T,xdata,ydata,Est0,res,varargin{:});

xdata = p.Results.xdata;
ydata = p.Results.ydata;
Est0 = p.Results.Est0;
res = p.Results.res(:);
B = round(p.Results.B);
solver = validatestring(p.Results.solver,{'lsqc','lsqn','ga','particle','psearch','surrogate','annealing','fmincon'});
Nstarts = round(p.Results.N);
noiseChoice = validatestring(p.Results.noise,{'poisson','quasipoisson','nbglobal'});
parallelRefit = p.Results.parallel;
Opt = p.Results.opt;
timer = p.Results.timer;

% The B internal refits must stay quiet. With no options, fmincon/ga/patternsearch/... print a
% result block per start, and on an interactive desktop that console rendering dominates the runtime
% -- a 40x8 fmincon bootstrap measured ~234 s this way versus ~54 s for the identical refits with
% Display off. Build a Display-off default for the chosen solver; a caller-supplied 'opt' is passed
% through untouched. Display is cosmetic, so this does not change the estimates for a given seed.
if isempty(Opt)
    Opt = local_quietOptions(solver);
end

% Orient ydata to outputs x samples (accept a single output passed as a column vector).
if size(ydata,2) ~= length(xdata) && size(ydata,1) == length(xdata)
    ydata = ydata.';
end
inputs = size(ydata,1);

cumulative = p.Results.cumulative;
if isempty(cumulative)
    cumulative = false(1,inputs);
end
if numel(cumulative) ~= inputs
    error('gsua_bootstrap:BadCumulative', ...
        'cumulative must have one entry per fitted output row (%d), got %d', inputs, numel(cumulative));
end

% margin/alpha for the refits: an explicit override is in GSUA_PE's raw convention; a table-
% recovered value is already offset by +1, so convert it back to raw before handing it to GSUA_PE.
if isempty(p.Results.margin)
    try
        marginRaw = T.Properties.CustomProperties.Margin - 1;
    catch
        marginRaw = 0;   % no recorded margin: refit with plain MSE/least-squares, like GSUA_PE's default
    end
else
    marginRaw = p.Results.margin;
end
if isempty(p.Results.alpha)
    try
        costAlpha = T.Properties.CustomProperties.Alpha;
    catch
        costAlpha = 2;
    end
else
    costAlpha = p.Results.alpha;
end

seed = p.Results.seed;
if ~isempty(seed)
    priorState = rng;
    cleanupObj = onCleanup(@() rng(priorState)); %#ok<NASGU>
    rng(seed);
end

% theta*: best pool member by res, and the pooled noise model from its residuals (shared with
% GSUA_NOISEFLOOR). Surface the degenerate-dispersion cases with this function's own identifiers.
[~,ibest] = min(res);
thetaStar = Est0(:,ibest).';
[f,muClampedGrid,phi,kGlobal,quasiFallback,nValid,pFree] = gsua_noisemodel(T,xdata,ydata,thetaStar,cumulative);
if nValid - pFree <= 0
    warning('gsua_bootstrap:LowDF', ...
        'Only %d valid (non-NaN, mu>0) pooled residuals for %d free parameters; the noise model may be unreliable.', ...
        nValid, pFree);
end
if quasiFallback && strcmp(noiseChoice,'quasipoisson')
    warning('gsua_bootstrap:QuasiPoissonDegenerate', ...
        'Pooled Pearson dispersion phi=%.6g <= 1; quasipoisson falls back to Poisson draws.', phi);
end
if isinf(kGlobal) && strcmp(noiseChoice,'nbglobal')
    warning('gsua_bootstrap:UnderdispersedNB', ...
        'Method-of-moments global-NB k is degenerate; nbglobal falls back to Poisson draws.');
end

% Draw the B synthetic datasets once (inputs x len x B), then refit each.
yb = gsua_noisedraw(noiseChoice,muClampedGrid,cumulative,ydata,phi,kGlobal,quasiFallback,B);

Np = size(Est0,1);
Est = zeros(Np,B);
costs = inf(1,B);
estCol = ['Est' solver];

if timer
    fprintf('gsua_bootstrap: %d parametric-bootstrap refits (%s, N=%d per replica, noise=%s)\n', ...
        B, solver, Nstarts, noiseChoice);
end

if parallelRefit
    parfor b = 1:B
        [Eb,cb] = local_refit(T,xdata,yb(:,:,b),solver,Nstarts,marginRaw,costAlpha,Opt,estCol);
        Est(:,b) = Eb;
        costs(b) = cb;
    end
    if timer, fprintf('gsua_bootstrap: done (%d refits, parallel)\n', B); end
else
    nextTick = 0;
    for b = 1:B
        [Eb,cb] = local_refit(T,xdata,yb(:,:,b),solver,Nstarts,marginRaw,costAlpha,Opt,estCol);
        Est(:,b) = Eb;
        costs(b) = cb;
        if timer && b/B*100 >= nextTick
            fprintf('  %3.0f%% (%d/%d)\n', b/B*100, b, B);
            nextTick = nextTick + 10;
        end
    end
end

out = struct();
out.thetaStar = thetaStar;
out.f = f;
out.phi = phi;
out.k = kGlobal;
out.noise = noiseChoice;
out.quasiFallback = quasiFallback;
out.exampleSynthetic = yb(:,:,1);
out.B = B;
out.solver = solver;
out.N = Nstarts;
out.margin = marginRaw;
out.alpha = costAlpha;
end


function [Eb,cb] = local_refit(T,xdata,yb_b,solver,Nstarts,marginRaw,costAlpha,Opt,estCol)
%LOCAL_REFIT Refit one synthetic dataset and return its best estimate + cost.
if isempty(Opt)
    [Tb,resb] = gsua_pe(T,xdata,yb_b,'solver',solver,'N',Nstarts, ...
        'margin',marginRaw,'alpha',costAlpha,'save',false,'timer',false);
else
    [Tb,resb] = gsua_pe(T,xdata,yb_b,'solver',solver,'N',Nstarts, ...
        'margin',marginRaw,'alpha',costAlpha,'opt',Opt,'save',false,'timer',false);
end
estimates = Tb.(estCol);     % Np x N, columns sorted best-first by gsua_pe
[cb,jx] = min(resb);
Eb = estimates(:,jx);
end


function opt = local_quietOptions(solver)
%LOCAL_QUIETOPTIONS Display-off options for the GSUA_PE solver, so the refits never flood the
%   console. Only the Display setting is touched; every other option stays at the solver default,
%   so the estimates are identical to the previous (option-less) path for a given seed.
switch solver
    case 'fmincon',   opt = optimoptions('fmincon','Display','off');
    case 'lsqc',      opt = optimoptions('lsqcurvefit','Display','off');
    case 'lsqn',      opt = optimoptions('lsqnonlin','Display','off');
    case 'ga',        opt = optimoptions('ga','Display','off');
    case 'particle',  opt = optimoptions('particleswarm','Display','off');
    case 'psearch',   opt = optimoptions('patternsearch','Display','off');
    case 'surrogate', opt = optimoptions('surrogateopt','Display','off');
    case 'annealing', opt = optimoptions('simulannealbnd','Display','off');
    otherwise,        opt = [];   % unknown solver: leave gsua_pe to use its own default
end
end
