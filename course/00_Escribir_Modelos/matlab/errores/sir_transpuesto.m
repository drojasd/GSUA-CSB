function sol = sir_transpuesto(pars, domain, opt)
% SIR_TRANSPUESTO El SIR continuo devolviendo sol.y transpuesta (tiempo x salidas).
%
%   Solo existe para el módulo 0: desde la versión 1.13 gsua_eval la
%   transpone sola (reconoce el eje del tiempo por numel(sol.x)), salvo
%   que haya tantas salidas como tiempos. No dependas de esa corrección.

sol = sir_continuo(pars, domain, opt);
sol.y = sol.y.';                 % Nt x 3: MAL
end
