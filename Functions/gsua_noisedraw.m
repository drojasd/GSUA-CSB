function yb = gsua_noisedraw(model,muClampedGrid,cumulative,ydata,phi,kGlobal,quasiFallback,B)
%GSUA_NOISEDRAW Draw B synthetic datasets under a count-noise model.
%
%   yb = gsua_noisedraw(model,muClampedGrid,cumulative,ydata,phi,kGlobal,quasiFallback,B)
%
%   Shared synthetic-data generator used by GSUA_NOISEFLOOR (which scores the true model against
%   the draws) and GSUA_BOOTSTRAP (which refits each draw). Factored out so both calibrations draw
%   from exactly the same noise models rather than keeping two copies in sync.
%
%   Inputs:
%     model          <-- 'poisson', 'quasipoisson' or 'nbglobal' (see GSUA_NOISEFLOOR for the
%                        variance definitions and when each degenerates to Poisson)
%     muClampedGrid  <-- inputs x len mean grid, clamped at 0 (from GSUA_NOISEMODEL)
%     cumulative     <-- 1 x inputs logical; a true row is drawn on the incident series and
%                        re-accumulated with cumsum
%     ydata          <-- inputs x len reference data; its NaN mask is reapplied to every draw
%     phi            <-- pooled Pearson dispersion (used by 'quasipoisson')
%     kGlobal        <-- method-of-moments global NB dispersion (used by 'nbglobal'; Inf -> Poisson)
%     quasiFallback  <-- true if quasi-Poisson must fall back to Poisson draws (phi<=1)
%     B              <-- number of synthetic datasets to draw
%
%   Output:
%     yb <-- inputs x len x B synthetic datasets
%
%   See also GSUA_NOISEMODEL, GSUA_NOISEFLOOR, GSUA_BOOTSTRAP.

[inputs,len] = size(muClampedGrid);
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
    otherwise
        error('gsua_noisedraw:UnknownModel','Unknown noise model: %s',model);
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
