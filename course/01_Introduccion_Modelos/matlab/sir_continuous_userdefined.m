function sol = sir_continuous_userdefined(params, domain, ~)
% SIR_CONTINUOUS_USERDEFINED Modelo epidémico SIR continuo resuelto con ode45.
%
%   sol = sir_continuous_userdefined(params, domain, opt)
%
%   Factores (Parámetros + Condiciones Iniciales):
%       params(1) -> beta  : Tasa de transmisión / contacto infectivo (1/(persona*día))
%       params(2) -> gamma : Tasa de recuperación (1/día)
%       params(3) -> S0    : Susceptibles iniciales en t0
%       params(4) -> I0    : Infectados iniciales en t0
%       params(5) -> R0_pop: Recuperados iniciales en t0 (opcional, por defecto 0)
%
%   Salidas (3 componentes):
%       Fila 1 -> S(t): Susceptibles
%       Fila 2 -> I(t): Infectados
%       Fila 3 -> R(t): Recuperados
%
%   Retorna:
%       Estructura sol compatible con GSUA_EVAL / GSUA_INTRP.

beta  = params(1);
gamma = params(2);
S0    = params(3);
I0    = params(4);

if length(params) >= 5
    R0_pop = params(5);
else
    R0_pop = 0;
end

y0 = [S0; I0; R0_pop];

% Opciones para precisión numérica del integrador ode45
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

% Integrar el sistema EDO continuo
sol = ode45(@(t, y) sir_rhs(t, y, beta, gamma), domain, y0, options);

% Asignar sol.solver = '' para compatibilidad con el interpolador gsua_intrp
sol.solver = '';
end

function dydt = sir_rhs(~, y, beta, gamma)
S = y(1);
I = y(2);

dSdt = -beta * S * I;
dIdt =  beta * S * I - gamma * I;
dRdt =  gamma * I;

dydt = [dSdt; dIdt; dRdt];
end
