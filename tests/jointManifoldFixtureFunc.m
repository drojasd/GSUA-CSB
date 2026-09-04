function y = jointManifoldFixtureFunc(pars)
%JOINTMANIFOLDFIXTUREFUNC Test fixture: deliberately confounded 3-parameter model.
%
%   y = jointManifoldFixtureFunc(pars)
%
%   Domain-less (Kind 6) user model for the gsua_dmatrix joint-sampling tests.
%   y = a*b*exp(-k*t) with pars = [a; b; k]. The data constrains only the PRODUCT
%   a*b and the rate k, so the identified manifold is exactly the hyperbola
%   a*b = const -- a curved ridge, known in closed form. That makes it the analytic
%   control for whether a joint sampler stays on the manifold: bootstrap draws must,
%   and a Gaussian fitted to a curved ridge must not.
%
%   Domain is hardcoded to match the tests' xdata = linspace(0,5,20) exactly.
d = linspace(0, 5, 20);
y = pars(1)*pars(2)*exp(-pars(3)*d);
end
