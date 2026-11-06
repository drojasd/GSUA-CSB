function y = twoScaleUserModel(pars)
%TWOSCALEUSERMODEL Two-series user model at very different scales (v1.18 multi-output example).
%
%   y = twoScaleUserModel(pars), with pars = [a b c d], returns a 2 x 40 matrix:
%     row 1 ~ O(1e3) :  a*exp(-b*t)        -- a decaying, count-like series
%     row 2 ~ O(1e-2):  c*(1-exp(-d*t))    -- a small, saturating, fraction-like series
%   over t = 1:40. The two rows differ in scale by ~5 orders of magnitude, and the factors span
%   1e-2..1e3, which is exactly the situation that motivates:
%     - per-row noise in gsua_bootstrap (count noise on row 1, gaussian/none on row 2),
%     - the automatic TypicalX scaling in gsua_pe (factors across many magnitudes), and
%     - per-series time-dependent sensitivity in gsua_sa.
%   Being a ONE-argument function, gsua_dataprep classifies it as a Kind-6 user model whose output
%   is a full matrix -- the case this release makes flow through gsua_ua/gsua_sa/gsua_pe/gsua_dia
%   without any scalar adapter.
%
%   See also GSUA_DATAPREP, GSUA_UA, GSUA_SA, GSUA_PE, GSUA_DIA, GSUA_BOOTSTRAP, GSUA_FULLPARS.
t = 1:40;
a = pars(1); b = pars(2); c = pars(3); d = pars(4);
y = [a*exp(-b*t); c*(1-exp(-d*t))];
end
