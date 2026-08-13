function y = noisefloorFixtureFunc(pars)
%NOISEFLOORFIXTUREFUNC Test fixture: 2-output count-like toy model for gsua_noisefloor tests.
%
%   y = noisefloorFixtureFunc(pars)
%
%   Domain-less (Kind 6) user model. Row 1 = pars(1)*exp(pars(2)*d) (incident-like, e.g. weekly
%   case counts), row 2 = cumsum of row 1 scaled by pars(3) (cumulative-like, e.g. running total).
%   Domain is hardcoded to d = linspace(0,4,15), matching the fixture tests' xdata exactly.
d = linspace(0,4,15);
row1 = pars(1) * exp(pars(2)*d);
row2 = pars(3) * cumsum(row1);
y = [row1; row2];
end
