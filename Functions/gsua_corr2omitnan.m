function r = gsua_corr2omitnan(a,b)
%GSUA_CORR2OMITNAN 2-D correlation coefficient over jointly non-NaN entries
%
%   r = gsua_corr2omitnan(a,b)
%
%   NaN-tolerant replacement for corr2(a,b): computes the standard 2-D
%   correlation coefficient using only the positions where both a and b
%   are non-NaN. Identical to corr2(a,b) when neither input contains
%   NaN. Returns 0 (no correlation credit, rather than erroring or
%   propagating NaN) if fewer than 2 jointly-valid entries remain.
%
%   Exists so GSUA_COSTF's correlation penalty term degrades gracefully
%   when real-world reference data has missing (NaN) periods, instead of
%   poisoning the whole objective with a single NaN.
%
%   Inputs:
%     a,b <-- equal-size vectors or matrices to correlate
%
%   Outputs:
%     r <-- 2-D correlation coefficient in [-1,1], or 0 if fewer than 2
%           jointly non-NaN entries remain
%
%   See also GSUA_COSTF, CORR2.
valid = ~isnan(a) & ~isnan(b);
a = a(valid);
b = b(valid);
if numel(a) < 2
    r = 0;
else
    r = corr2(a,b);
end
end
