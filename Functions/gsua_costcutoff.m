function keep = gsua_costcutoff(cost,method,rtol,atol,gapRatio,minKeep)
%GSUA_COSTCUTOFF Decide which repeated-estimation runs to keep based on fit cost
%
%   keep = gsua_costcutoff(cost,method,rtol,atol,gapRatio,minKeep)
%
%   Shared fit-quality filter used by GSUA_IA/GSUA_DIA to drop repeated
%   estimation runs that converged to a much worse cost than the best
%   one, before any identifiability statistic is computed from their
%   parameter values. A run that converged to a bad local optimum
%   contributes essentially arbitrary parameter values -- left in, it can
%   make a genuinely well-identified parameter look poorly identified
%   purely because an optimizer run failed, not because of real
%   non-identifiability.
%
%   Inputs:
%     cost     <-- 1xN (or Nx1) fit cost, one per repeated estimation run
%                  (e.g. the res output of gsua_pe)
%     method   <-- 'rtol' -- keep cost <= min(cost)*(1+rtol)+atol.
%                  Simple and predictable, but a fixed tolerance is a
%                  somewhat arbitrary answer to a genuinely data-dependent
%                  question.
%                  'gap' -- keep runs up to the largest jump in the
%                  *sorted* costs (in log space), the automated version
%                  of reading a waterfall plot by eye: a cluster of
%                  converged runs, a jump, then a scattered tail of
%                  failures. Only cuts if that jump is at least gapRatio
%                  times the *median* gap between sorted costs, so a
%                  smooth continuum of similar costs is left untouched
%                  rather than sliced arbitrarily.
%     rtol     <-- relative tolerance for method='rtol'
%     atol     <-- absolute tolerance added to the 'rtol' threshold, so
%                  the filter doesn't become vacuously strict when the
%                  best cost is at or near zero (e.g. a synthetic,
%                  noiseless fit)
%     gapRatio <-- salience ratio for method='gap'
%     minKeep  <-- minimum number of runs to keep regardless of method or
%                  how aggressive its cutoff would otherwise be -- a
%                  degenerate 0- or 1-point result would otherwise
%                  silently propagate into meaningless/NaN correlation,
%                  confidence-interval, and clustering results downstream
%                  instead of a clear signal that the filter was too
%                  strict
%
%   Output:
%     keep <-- 1xN logical, true for runs to keep
%
%   See also GSUA_IA, GSUA_DIA, GSUA_PE.

cost = cost(:)';
n = numel(cost);

switch method
    case 'rtol'
        threshold = min(cost)*(1+rtol) + atol;
        keep = cost <= threshold;
    case 'gap'
        keep = true(1,n);
        if n >= 3
            [sortedCost,gapOrder] = sort(cost);
            logCost = log(max(sortedCost,realmin));
            gaps = diff(logCost);
            maxGap = max(gaps);
            medGap = median(gaps);
            if maxGap > 0 && medGap > 0 && maxGap >= gapRatio*medGap
                [~,cutRel] = max(gaps); % keep sortedCost(1:cutRel)
                keepSorted = false(1,n);
                keepSorted(1:cutRel) = true;
                keep(gapOrder) = keepSorted;
            end
        end
    otherwise
        error('gsua_costcutoff:UnknownMethod', ...
            'Unknown cost method: %s (expected ''rtol'' or ''gap'')',method)
end

minKeep = min(minKeep,n);
if sum(keep) < minKeep
    [~,costOrder] = sort(cost);
    rankOf = zeros(1,n);
    rankOf(costOrder) = 1:n;
    keep = rankOf <= minKeep;
end

end
