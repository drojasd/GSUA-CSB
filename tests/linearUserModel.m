function y = linearUserModel(pars)
%LINEARUSERMODEL y = x1 + 2*x2 as a one-argument user model (one row per sample).
%   With (x1,x2) standard normal and corr(x1,x2) = rho, the Shapley effects are
%   Sh1 = (Var Y - 4(1-rho^2) + (1-rho^2)) / (2 Var Y), Sh2 = 1 - Sh1, Var Y = 5 + 4 rho;
%   for rho = 0.5: Sh = [0.3393 0.6607]. Test fixture for gsuaShapleyTest ('Given' estimator).
if isvector(pars)
    pars = reshape(pars,1,[]);
end
y = pars(:,1) + 2*pars(:,2);
end
