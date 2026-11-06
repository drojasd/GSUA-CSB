function yb = gsua_noisedraw(model,muClampedGrid,cumulative,ydata,phi,kGlobal,quasiFallback,B,sigmaRow)
%GSUA_NOISEDRAW Draw B synthetic datasets under a count-noise model (optionally per output row).
%
%   yb = gsua_noisedraw(model,muClampedGrid,cumulative,ydata,phi,kGlobal,quasiFallback,B)
%   yb = gsua_noisedraw(model,...,B,sigmaRow)
%
%   Shared synthetic-data generator used by GSUA_NOISEFLOOR (which scores the true model against
%   the draws) and GSUA_BOOTSTRAP (which refits each draw). Factored out so both calibrations draw
%   from exactly the same noise models rather than keeping two copies in sync.
%
%   Inputs:
%     model          <-- noise model. Either a single name applied to every row -- 'poisson',
%                        'quasipoisson' or 'nbglobal' (count models; see GSUA_NOISEFLOOR) -- or a
%                        1 x inputs cellstr giving a PER-ROW model, where each entry may also be
%                        'gaussian' (additive N(0,sigmaRow(i)^2), for fraction/continuous rows) or
%                        'none' (leave the row at its fitted mean, no perturbation). A per-row model
%                        is what lets a joint count+fraction fit avoid drawing count noise (and
%                        zeros) on a 0-1 fraction row.
%     muClampedGrid  <-- inputs x len mean grid, clamped at 0 (from GSUA_NOISEMODEL)
%     cumulative     <-- 1 x inputs logical; a true row is drawn on the incident series and
%                        re-accumulated with cumsum
%     ydata          <-- inputs x len reference data; its NaN mask is reapplied to every draw
%     phi            <-- pooled Pearson dispersion (used by 'quasipoisson')
%     kGlobal        <-- method-of-moments global NB dispersion (used by 'nbglobal'; Inf -> Poisson)
%     quasiFallback  <-- true if quasi-Poisson must fall back to Poisson draws (phi<=1)
%     B              <-- number of synthetic datasets to draw
%     sigmaRow       <-- (optional) 1 x inputs additive-Gaussian std per row, used only by 'gaussian'
%                        rows. Default: zeros (a 'gaussian' row with sigma 0 is then just its mean).
%
%   Output:
%     yb <-- inputs x len x B synthetic datasets
%
%   See also GSUA_NOISEMODEL, GSUA_NOISEFLOOR, GSUA_BOOTSTRAP.

[inputs,len] = size(muClampedGrid);
if nargin < 9 || isempty(sigmaRow)
    sigmaRow = zeros(1,inputs);
end

% Homogeneous single count model -> keep the original vectorized draw byte-for-byte, so
% GSUA_NOISEFLOOR and the default GSUA_BOOTSTRAP reproduce exactly as before for a given seed.
homogeneous = (ischar(model) || isstring(model)) || (iscell(model) && numel(unique(model))==1);
if homogeneous
    if iscell(model), m1 = model{1}; else, m1 = char(model); end
    if any(strcmp(m1,{'poisson','quasipoisson','nbglobal'}))
        muRep = repmat(muClampedGrid,[1 1 B]);
        switch m1
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
        return
    end
end

% Per-row path: normalize model to a 1 x inputs cellstr and draw each row under its own model.
if ischar(model) || isstring(model)
    model = repmat({char(model)},1,inputs);
elseif iscell(model) && numel(model)==1
    model = repmat(model,1,inputs);
end
if numel(model) ~= inputs
    error('gsua_noisedraw:BadModel', ...
        'noise must be a single model or one model per output row (%d), got %d.', inputs, numel(model));
end

yb = zeros(inputs,len,B);
for i = 1:inputs
    mu_i = muClampedGrid(i,:);
    muRep = repmat(mu_i,[1 1 B]);
    switch lower(char(model{i}))
        case 'poisson'
            S = poissrnd(muRep);
        case 'quasipoisson'
            if quasiFallback
                S = poissrnd(muRep);
            else
                kPoint = max(mu_i ./ (phi-1), realmin);
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
        case 'gaussian'
            S = muRep + sigmaRow(i)*randn(size(muRep));   % additive, not clamped (fractions)
        case 'none'
            S = muRep;                                    % no perturbation: fitted mean
        otherwise
            error('gsua_noisedraw:UnknownModel','Unknown noise model: %s',char(model{i}));
    end
    S(isnan(S)) = 0;
    if cumulative(i)
        yb(i,:,:) = cumsum(S,2);
    else
        yb(i,:,:) = S;
    end
    yb(i,isnan(ydata(i,:)),:) = NaN;
end
end
