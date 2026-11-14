function y = ishigamiUserModel(pars)
%ISHIGAMIUSERMODEL Ishigami test function (a=7, b=0.1) as a one-argument user model.
%   y = ishigamiUserModel(pars) with pars = [x1 x2 x3] (one row per sample) on [-pi,pi]^3.
%   Analytic first-order / total Sobol' indices: S = [0.3139 0.4424 0], ST = [0.5576 0.4424 0.2437];
%   Shapley effects: Sh = [0.4358 0.4424 0.1218] (the x1-x3 interaction split evenly).
%   Test fixture for gsuaShapleyTest.
a = 7; b = 0.1;
if isvector(pars)
    pars = reshape(pars,1,[]);      % gsua_deval passes one parameter set as a column
end
y = sin(pars(:,1)) + a*sin(pars(:,2)).^2 + b*(pars(:,3).^4).*sin(pars(:,1));
end
