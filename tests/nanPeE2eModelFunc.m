function y = nanPeE2eModelFunc(pars)
%NANPEE2EMODELFUNC Test fixture: 2-output toy model for gsuaCostfNanTolerantTest.
%
%   y = nanPeE2eModelFunc(pars)
%
%   Domain-less (Kind 6) user model, used only by
%   gsuaCostfNanTolerantTest's end-to-end gsua_pe check. Domain is
%   hardcoded to match that test's xdata = linspace(0,5,20) exactly.
d = linspace(0, 5, 20);
y = [pars(1)*exp(-pars(2)*d); pars(3)*d];
end
