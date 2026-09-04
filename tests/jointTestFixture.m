function [T,pool] = jointTestFixture(nPool)
%JOINTTESTFIXTURE Table plus an accepted-estimate pool on a known curved manifold.
%
%   [T,pool] = jointTestFixture()
%   [T,pool] = jointTestFixture(nPool)
%
%   Shared fixture for the gsua_dmatrix joint-sampling tests. Builds a 3-parameter
%   table over jointManifoldFixtureFunc (y = a*b*exp(-k*t)) and a pool of nPool
%   accepted estimates lying EXACTLY on the hyperbola a*b = 2, which is the manifold
%   the data of that model actually identifies. Because the manifold is known in
%   closed form, a sampler's fidelity to it is measurable rather than a matter of
%   opinion, and because it is curved, matching pairwise correlation is not enough to
%   stay on it.
%
%   The k row varies deterministically (sin, not randn) so the fixture never touches
%   the global RNG and cannot perturb the seed-reproducibility tests.
if nargin < 1
    nPool = 40;
end
T = gsua_userdefined('jointManifoldFixtureFunc', [0.3 3; 0.3 3; 0.1 1.2]);
a = linspace(0.5, 2.5, nPool);
pool = [a; 2./a; 0.5 + 0.02*sin(1:nPool)];
T.Est = pool;
end
