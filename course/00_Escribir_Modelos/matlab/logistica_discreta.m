function sol = logistica_discreta(pars, domain, ~)
% LOGISTICA_DISCRETA Crecimiento logístico en tiempo discreto (modelo de usuario).
%
%   sol = logistica_discreta(pars, domain, opt)
%
%   Ecuación en diferencias:  N(t+1) = N(t) + r*N(t)*(1 - N(t)/K)
%
%   Factores, en este orden (primero parámetros, después condición inicial):
%       pars(1) -> r  : tasa de crecimiento por paso (adimensional)
%       pars(2) -> K  : capacidad de carga (individuos)
%       pars(3) -> N0 : población en t0 (individuos)
%
%   domain = [t0, tf]: lo fija 'domain' en gsua_dataprep.
%   El tercer argumento (opt) no se usa aquí, pero DEBE declararse:
%   gsua_eval siempre llama al modelo con tres argumentos.
%
%   Salida: sol.x (1 x Nt) con los pasos t0:tf
%           sol.y (1 x Nt) con N(t): una fila por salida.

r  = pars(1);
K  = pars(2);
N0 = pars(3);

t  = domain(1):domain(2);       % un punto por paso de tiempo
N  = zeros(1, numel(t));        % 1 salida -> 1 fila
N(1) = N0;
for k = 1:numel(t)-1
    N(k+1) = N(k) + r*N(k)*(1 - N(k)/K);
end

sol.x = t;
sol.y = N;
end
