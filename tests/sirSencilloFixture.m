function [T, xdata, ydata, truth] = sirSencilloFixture()
%SIRSENCILLOFIXTURE Discrete SIR (SIRS) practical-identifiability test case.
%
%   [T, xdata, ydata, truth] = sirSencilloFixture()
%
%   Builds the summary table, true parameters and negative-binomial (k=10, rng(1)) observation of
%   I(t) used by the gsua_bootstrap / gsua_ia regression tests. This is the course "Modelacion
%   Experimental" (EAFIT) discrete-SIR case: multistart alone collapses S0 to ~2982.7 (truth
%   5000), which the parametric bootstrap and the flat-parameter index are designed to expose.
%
%   Requires sir_sencillo on the path (Examples/). Factors [S0 I0 R0 beta gamma alpha]; only I(t)
%   (output 2) is observed, on xdata = 1:50.

ranges = [1e3 1e4; 0 10; 10 100; 0 1; 0 1; 0 1];
[T,~] = gsua_dataprep('sir_sencillo', ranges, 'domain',[1 200], ...
    'names',{'S0','I0','R0','beta','gamma','alpha'}, 'out_names',{'S','I','R'});
T.Properties.CustomProperties.output = 2;   % observe only Infected

truth = [5000;5;50;0.5;0.2;0.01];
xdata = 1:50;
yI = squeeze(gsua_eval(truth,T,xdata,[],false,false,false));
yI = yI(:).';

rng(1);
k = 10;                                  % negative-binomial dispersion
ydata = nbinrnd(k, k./(k+yI));           % mean yI, var yI + yI.^2/k
end
