function sol = logistica_dos_entradas(pars, domain)
% LOGISTICA_DOS_ENTRADAS La logística discreta declarada con DOS argumentos.
%
%   Solo existe para el módulo 0: gsua_dataprep la acepta, pero gsua_eval
%   llama al modelo como modelo(pars, domain, opt) y falla con
%   "Too many input arguments". La versión correcta es LOGISTICA_DISCRETA.

N0 = pars(1); r = pars(2); K = pars(3);
t = domain(1):domain(2);
N = zeros(1, numel(t));
N(1) = N0;
for k = 1:numel(t)-1
    N(k+1) = N(k) + r*N(k)*(1 - N(k)/K);
end
sol.x = t;
sol.y = N;
end
