function [M,T] = gsua_dmatrix(Table,N,varargin)
% Function for design of experiments (factor space sampling)
%
% [M,T2]=gsua_dmatrix(T,N)
% Parameters:
% T <-- summary table from gsua_dataprep
% N <-- number of samples
% Outpus:
% M  <-- design matrix of NxNp for later routines
% T2 <-- summary table with fixed parameters actualized
% Additional features:
% You can choose a method for factor space sampling between uniform
% distribution design and latin hypercube design (default). To switch
% between methods use the paired feature 'Method' and 'Uniform',
% 'Sobol', or 'LatinHypercube'.
% Also you can visualize the sampling result using the paired feature
% 'Show', 'on'
% M=gsua_dmatrix(T,N,'Method','Sobol','Show','on')
%
% MARGINAL VS JOINT SAMPLING
% The methods above are MARGINAL: each parameter is drawn independently
% inside its own range. For a model with parameter confounding (correlated,
% marginally-non-identifiable parameters) that is wrong in a way that
% matters -- independent draws leave the identified manifold, so the
% resulting uncertainty band is both too wide and badly centered, with a
% median trajectory resembling no actual good fit. GSUA_COVMETRIC reports
% this condition as cost_band far above cost_data.
%
% 'Method','Joint' instead draws whole parameter VECTORS from an ensemble of
% accepted estimates, so the correlation structure among parameters
% survives:
% M=gsua_dmatrix(Tia,N,'Method','Joint')
%
% The ensemble is taken from Tia.Est (set by GSUA_IA, already filtered by
% fit quality, restricted to the dominant cluster and outlier-cleaned), or
% supplied explicitly with 'Pool' as an Np x nPool matrix.
%
% Joint sampling also fixes a second, separate problem. GSUA_IA replaces
% T.Range with a distribution-free confidence interval OF THE MEDIAN of the
% pool (see GSUA_MEDIANCI) -- a statement about where the pool's centre
% lies, which NARROWS as the pool grows. It is not the spread of parameter
% values consistent with the data, so sampling it understates uncertainty.
% Drawing pool vectors inherits the pool's own empirical spread instead.
% For that reason 'Clip' defaults to no clipping: clipping joint draws back
% to Tia.Range would reimpose exactly the narrowness this avoids. Pass the
% ORIGINAL model's Range (physical/feasible bounds) if clipping is wanted.
%
% 'JointType' selects how the ensemble is turned into samples:
%   'Bootstrap'       (default) resample whole pool columns with
%                     replacement. Exact correlations, no distributional
%                     assumption, cannot leave the manifold, and unaffected
%                     by a rank-deficient covariance. Draws are limited to
%                     the nPool distinct vectors actually in the pool.
%   'SmoothBootstrap' as above plus a variance-corrected Gaussian kernel, so
%                     draws are new points rather than repeats while staying
%                     near the manifold. Bandwidth defaults to Silverman's
%                     multivariate rule; override with 'Bandwidth'.
%   'Gaussian'        ensemble mean and covariance, sampled by
%                     eigendecomposition so a rank-deficient covariance
%                     (nPool < number of free parameters) is handled by
%                     sampling only the non-null eigendirections. Assumes
%                     the manifold is linear -- on a CURVED ridge a Gaussian
%                     places mass beside it, off the manifold.
%
% CAVEATS
% When GSUA_IA detects multiple global minima, Tia.Est holds the DOMINANT
% cluster only (pooling separated basins into one interval is not
% meaningful), so a joint band built from it is conditional on that basin.
% Use clusterInfo.Idx from GSUA_IA to build per-basin bands instead.
% A pool of fewer than 'MinPoolN' (default 5) vectors is refused: correlation
% from very few points is not merely noisy but degenerate -- two points
% always correlate at exactly +-1 -- so the result would be confidently
% wrong rather than approximately right. Such a call warns and falls back to
% marginal Latin hypercube sampling.
%
% See also GSUA_DATAPREP, GSUA_SA, GSUA_UA, GSUA_PARDEVAL, GSUA_IA,
% GSUA_MEDIANCI, GSUA_COVMETRIC.
p=inputParser;
defaultShow='off';
validShow={'off' 'on'};
checkShow = @(x) any(validatestring(x,validShow));
defaultMethod='LatinHypercube';
validMethod={'LatinHypercube','Uniform','Sobol','Joint'};
checkMethod = @(x) any(validatestring(x,validMethod));
validJointType={'Bootstrap','SmoothBootstrap','Gaussian'};
checkJointType = @(x) any(validatestring(x,validJointType));

addRequired(p,'Table');
addRequired(p,'N',@isnumeric);
addParameter(p,'Method',defaultMethod,checkMethod);
addParameter(p,'Show',defaultShow,checkShow);
addParameter(p,'JointType','Bootstrap',checkJointType);
addParameter(p,'Pool',[],@isnumeric);
addParameter(p,'MinPoolN',5,@(x) isnumeric(x) && isscalar(x) && x>=2);
addParameter(p,'Bandwidth',[],@(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>=0));
addParameter(p,'Clip',[],@(x) isempty(x) || (isnumeric(x) && size(x,2)==2));
addParameter(p,'Seed',[],@(x) isempty(x) || isnumeric(x));

parse(p,Table,N,varargin{:})
T=p.Results.Table;
N=p.Results.N;
method=validatestring(p.Results.Method,validMethod);
show=p.Results.Show;
jointType=validatestring(p.Results.JointType,validJointType);
pool=p.Results.Pool;
minPoolN=p.Results.MinPoolN;
bandwidth=p.Results.Bandwidth;
clipRange=p.Results.Clip;
seed=p.Results.Seed;
Range=T.Range';
Np=size(T,1);

if ~strcmp(method,'Joint') && ~ismember('JointType',p.UsingDefaults)
    error('gsua_dmatrix:JointTypeWithoutJoint',...
        '''JointType'' only applies to ''Method'',''Joint'' (got ''%s'').',method)
end

if ~isempty(seed)
    priorRngState=rng;
    rngCleanup=onCleanup(@() rng(priorRngState));
    rng(seed);
end
try
    Table2=T.Properties.CustomProperties;
catch
    TP=load('ATable.mat');
    Table2=TP.Table2;
end

N=floor(N);

if strcmp(method,'Joint')
    [M,method]=local_jointmatrix(T,N,Np,Range,jointType,pool,minPoolN,bandwidth,clipRange);
elseif strcmp(Table2.rMethod,{'normal'})
    disp('Given ranges represent a normal distribution')
    M = zeros(N,Np);
    for k = 1:Np
        pdfun = makedist('normal','mu',Range(1,k),'sigma',Range(2,k));
        M(:,k) = random(pdfun,1,N); % Normal distribution between two values given by the range of parameters
    end
else
    switch method
        case 'LatinHypercube'
            M = lhsdesign(N,Np).*( ones(N,1)*(Range(2,:)-Range(1,:)) ) + ones(N,1)*Range(1,:);
        case 'Uniform'
            M = zeros(N,Np);
            for k = 1:Np
                pdfun = makedist('Uniform','lower',Range(1,k),'upper',Range(2,k));
                M(:,k) = random(pdfun,1,N); % uniform distribution between two values given by the range of parameters
            end
        case 'Sobol'
            p=sobolset(Np);
            p = scramble(p,'MatousekAffineOwen');
            M = net(p,N).*( ones(N,1)*(Range(2,:)-Range(1,:)) ) + ones(N,1)*Range(1,:);
     end
end
try
    T.Properties.CustomProperties.Fixed=Range(1,:)==Range(2,:);
catch 
    Table2.Fixed=Range(1,:)==Range(2,:);
    save('ATable','Table2');
end

if strcmp(show,'on')
    figure
    clf
    gsua_plot('ScatterParameter',T,method,M)
end
end

function [M,methodUsed]=local_jointmatrix(T,N,Np,Range,jointType,pool,minPoolN,bandwidth,clipRange)
%LOCAL_JOINTMATRIX Correlation-preserving sampling from an accepted-estimate ensemble.
methodUsed='Joint';

if isempty(pool)
    if ~ismember('Est',T.Properties.VariableNames)
        error('gsua_dmatrix:NoPool',...
            ['''Method'',''Joint'' needs an ensemble of accepted estimates. Pass one as '...
            '''Pool'' (Np x nPool), or use a table from gsua_ia, which stores it as T.Est.'])
    end
    pool=T.Est;
end

if size(pool,1)~=Np
    error('gsua_dmatrix:PoolSizeMismatch',...
        'Pool has %u rows but the table has %u parameters; Pool must be Np x nPool.',...
        size(pool,1),Np)
end

free=Range(1,:)~=Range(2,:);
% Start every sample at the fixed value and overwrite only the free columns, matching
% how the marginal branches leave fixed parameters constant.
M=ones(N,1)*Range(1,:);
if ~any(free)
    return
end

poolFree=pool(free,:);
keep=all(isfinite(poolFree),1);
if ~any(keep)
    error('gsua_dmatrix:PoolAllNonFinite',...
        'Every pool column holds a non-finite value in at least one free parameter.')
end
if ~all(keep)
    warning('gsua_dmatrix:PoolNonFinite',...
        '%u of %u pool columns hold non-finite values and were dropped.',...
        sum(~keep),numel(keep))
    poolFree=poolFree(:,keep);
end

nPool=size(poolFree,2);
if nPool<minPoolN
    warning('gsua_dmatrix:PoolTooSmall',...
        ['Pool has %u vectors (fewer than MinPoolN = %u). Correlation estimated from so few '...
        'points is degenerate rather than merely noisy, so a joint sample would be confidently '...
        'wrong; falling back to marginal Latin hypercube sampling.'],nPool,minPoolN)
    M=lhsdesign(N,Np).*( ones(N,1)*(Range(2,:)-Range(1,:)) ) + ones(N,1)*Range(1,:);
    methodUsed='LatinHypercube';
    return
end

dFree=size(poolFree,1);
switch jointType
    case 'Bootstrap'
        X=poolFree(:,randi(nPool,1,N));
    case 'SmoothBootstrap'
        if isempty(bandwidth)
            % Silverman's multivariate rule, so the default needs no arbitrary constant.
            bandwidth=(4/((dFree+2)*nPool))^(1/(dFree+4));
        end
        L=local_covfactor(poolFree);
        mu=mean(poolFree,2);
        base=poolFree(:,randi(nPool,1,N));
        % The 1/sqrt(1+h^2) is what keeps this a smoothing of the ensemble rather than an
        % inflation of it: without it the added kernel would scale the covariance by (1+h^2).
        X=mu+(base-mu+bandwidth*(L*randn(size(L,2),N)))/sqrt(1+bandwidth^2);
    case 'Gaussian'
        L=local_covfactor(poolFree);
        mu=mean(poolFree,2);
        X=mu+L*randn(size(L,2),N);
end

M(:,free)=X';

if ~isempty(clipRange)
    if size(clipRange,1)~=Np
        error('gsua_dmatrix:ClipSizeMismatch',...
            'Clip has %u rows but the table has %u parameters; Clip must be Np x 2.',...
            size(clipRange,1),Np)
    end
    M=min(max(M,ones(N,1)*clipRange(:,1)'),ones(N,1)*clipRange(:,2)');
end
end

function L=local_covfactor(poolFree)
%LOCAL_COVFACTOR Square-root factor of the ensemble covariance, restricted to its non-null
%eigendirections so a rank-deficient covariance (nPool <= number of free parameters) stays usable.
C=cov(poolFree');
C=(C+C')/2;
[V,D]=eig(C);
d=real(diag(D));
tol=max(d)*numel(d)*eps;
keep=d>max(tol,0);
if ~any(keep)
    L=zeros(size(C,1),1);
    return
end
L=V(:,keep)*diag(sqrt(d(keep)));
end