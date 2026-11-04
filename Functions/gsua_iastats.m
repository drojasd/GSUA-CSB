function [lb,up,ind,med,RHO,boxin,corrin,extrin] = gsua_iastats(statT_est,origRange,indexMethod,ciMethod,correction)
%GSUA_IASTATS Shared range + identifiability-index computation for GSUA_IA/GSUA_DIA.
%
%   [lb,up,ind,med,RHO,boxin,corrin,extrin] = ...
%       gsua_iastats(statT_est,origRange,indexMethod,ciMethod,correction)
%
%   Computes, for a single (presumed-unimodal) point set, the confidence interval of the median,
%   the per-parameter identifiability index, and the correlation matrix. Factored out so GSUA_IA,
%   GSUA_DIA and the per-cluster statistics all use one implementation.
%
%   Inputs:
%     statT_est   <-- Np x n point set (one column per estimation)
%     origRange   <-- Np x 2 ORIGINAL parameter bounds (the index denominator)
%     indexMethod <-- 'spread' (estimate 2.5-97.5 percentile dispersion over origRange) or 'ci'
%                     (legacy: median-CI width over origRange)
%     ciMethod    <-- 'orderstat' (distribution-free GSUA_MEDIANCI, GSUA_IA's) or 'normal'
%                     (parametric med +/- 1.96*sqrt(pi/2)*std/sqrt(n), GSUA_DIA's)
%     correction  <-- (optional, logical) clip the CI to origRange. Default: false
%
%   Outputs:
%     lb,up   <-- Np x 1 lower/upper confidence interval of the median
%     ind     <-- 1 x Np identifiability index (0 = well identified, 1 = poorly identified)
%     med     <-- Np x 1 median estimate
%     RHO     <-- Np x Np parameter correlation matrix
%     boxin   <-- Np x 1 interval-width term of the index
%     corrin  <-- 1 x Np mean-correlation term
%     extrin  <-- 1 x Np strong-correlation-count term
%
%   See also GSUA_IA, GSUA_DIA, GSUA_MEDIANCI.

if nargin<5 || isempty(correction), correction=false; end

med = median(statT_est,2);
switch ciMethod
    case 'orderstat'
        [lb,up] = gsua_medianCI(statT_est,0.05);
    otherwise   % 'normal' (parametric, GSUA_DIA's)
        n = size(statT_est,2);
        desv = std(statT_est,[],2);
        lb = med - 1.96*sqrt(pi/2)*desv/sqrt(n);
        up = med + 1.96*sqrt(pi/2)*desv/sqrt(n);
end

lo = origRange(:,1);
hi = origRange(:,2);
if correction
    lb(lb<lo) = lo(lb<lo);
    up(up>hi) = hi(up>hi);
end

RHO = corr(statT_est');
rangeW = hi - lo;
switch indexMethod
    case 'spread'
        spreadW = prctile(statT_est,97.5,2) - prctile(statT_est,2.5,2);
        boxin = min(spreadW ./ rangeW, 1);
    otherwise   % 'ci' (legacy)
        boxin = (up-lb) ./ rangeW;
end
len = length(boxin);
corrin = sum(abs(RHO))/len;
extrin = sum(abs(RHO)>0.5)/len;
ind = (2*boxin'+corrin+extrin)/4;
end
