function [T, xdata, ydata, truePars] = noisefloorTestFixture()
%NOISEFLOORTESTFIXTURE Shared T/xdata/ydata fixture for gsua_noisefloor tests.
%
%   [T,xdata,ydata,truePars] = noisefloorTestFixture()
%
%   Builds a gsua_userdefined table around noisefloorFixtureFunc (2-output: row 1 incident-like,
%   row 2 cumulative-like), a fixed noiseless reference dataset at truePars, and runs gsua_pe once
%   with a nonzero margin so T carries the Margin/Alpha CustomProperties gsua_noisefloor needs.
xdata = linspace(0,4,15);
truePars = [4, 0.4, 0.8];
ydata = noisefloorFixtureFunc(truePars);

T = gsua_userdefined('noisefloorFixtureFunc', [0.1 20; 0.01 2; 0.01 5]);
[T,~] = gsua_pe(T, xdata, ydata, 'N', 1, 'solver', 'fmincon', 'margin', 0.1, ...
    'ipoint', truePars);
end
