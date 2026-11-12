function sol = sir_continuo_crudo(pars, domain, ~)
% SIR_CONTINUO_CRUDO El mismo SIR que SIR_CONTINUO, pero devuelve ode45 tal cual.
%
%   Solo existe para comparar en el módulo 0: sol.x son los pasos adaptativos
%   de ode45 (pocos y desiguales) y gsua_eval interpola linealmente entre ellos.

S0    = pars(1);
I0    = pars(2);
beta  = pars(3);
gamma = pars(4);
y0    = [S0; I0; 0];

rhs = @(t, y) [ -beta*y(1)*y(2);
                 beta*y(1)*y(2) - gamma*y(2);
                 gamma*y(2) ];

sol = ode45(rhs, domain, y0);   % tolerancias por defecto
end
