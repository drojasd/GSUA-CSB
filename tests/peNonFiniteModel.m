function sol = peNonFiniteModel(pars, domain, ~)
%PENONFINITEMODEL Test fixture: a Kind-5 model that returns non-finite output.
%
%   Used by gsuaPeDiagnosticsTest to exercise gsua_pe's pre-flight check: a model
%   whose output is not finite at the nominal parameters must make gsua_pe fail early
%   with GSUA:gsua_pe:ModelNotFinite instead of silently returning an all-Inf res.
t = linspace(domain(1), domain(2), 50);
sol.x = t;
sol.y = nan(1, numel(t)) * pars(1);   % non-finite on purpose
sol.solver = '';
end
