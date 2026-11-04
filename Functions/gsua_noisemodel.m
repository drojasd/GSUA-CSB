function [f,muClampedGrid,phi,kGlobal,quasiFallback,nValid,pFree] = gsua_noisemodel(T,xdata,ydata,thetaStar,cumulative)
%GSUA_NOISEMODEL Simulate theta* and estimate a pooled count-noise model from its residuals.
%
%   [f,muClampedGrid,phi,kGlobal,quasiFallback,nValid,pFree] = ...
%       gsua_noisemodel(T,xdata,ydata,thetaStar,cumulative)
%
%   Shared dispersion-estimation step behind GSUA_NOISEFLOOR and GSUA_BOOTSTRAP: simulate the best
%   fit, pool its residuals across every fitted output, and estimate the quasi-Poisson (Pearson
%   phi) and global negative-binomial (method-of-moments k) dispersion. Deliberately emits NO
%   warnings -- it returns the raw diagnostics (nValid, pFree, and the degeneracy flags implied by
%   kGlobal/quasiFallback) so each caller can warn with its OWN message identifier
%   (GSUA_NOISEFLOOR keeps its 'gsua_noisefloor:...' ids; GSUA_BOOTSTRAP uses its own).
%
%   Inputs:
%     T          <-- summary table (supplies T.Range for the free-parameter count)
%     xdata      <-- domain the pool was fit against
%     ydata      <-- inputs x len reference data (may contain NaN)
%     thetaStar  <-- 1 x Np best-fit parameter vector
%     cumulative <-- (optional) 1 x inputs logical; a true row is differenced before dispersion
%                    estimation. Default: all false.
%
%   Outputs:
%     f             <-- inputs x len simulated curves at thetaStar
%     muClampedGrid <-- inputs x len mean grid (incident for cumulative rows), clamped at 0
%     phi           <-- pooled Pearson dispersion sum(r^2/mu)/(nValid-pFree)
%     kGlobal       <-- method-of-moments global NB k (Inf if degenerate / underdispersed)
%     quasiFallback <-- true if phi <= 1 (quasi-Poisson must fall back to Poisson)
%     nValid        <-- number of valid (finite, mu>0) pooled residuals used
%     pFree         <-- number of free parameters (for the dispersion degrees of freedom)
%
%   See also GSUA_NOISEDRAW, GSUA_NOISEFLOOR, GSUA_BOOTSTRAP.

inputs = size(ydata,1);
len = length(xdata);
if nargin < 5 || isempty(cumulative)
    cumulative = false(1,inputs);
end

f = gsua_deval(thetaStar,T,xdata);   % inputs x len, assumed NaN-free (a clean simulation)
pFree = sum(T.Range(:,1) ~= T.Range(:,2));

% Per-output incident/cumulative transform, built once and reused for both dispersion estimation
% and synthetic generation. For a cumulative row, diff([0,row]) lets NaN propagate/smear naturally
% across any real gap in ydata.
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

% Pooled dispersion across every fitted output (not per-row): a single phi/k.
valid = isfinite(rGrid) & isfinite(muGrid) & muGrid > 0;
R = rGrid(valid);
MU = muGrid(valid);
nValid = numel(R);
denom = max(nValid - pFree,1);
phi = sum(R.^2 ./ MU) / denom;

kDenom = sum(R.^2) - sum(MU);
if kDenom <= 0
    kGlobal = Inf;
else
    kGlobal = sum(MU.^2) / kDenom;
end

quasiFallback = phi <= 1;
end
