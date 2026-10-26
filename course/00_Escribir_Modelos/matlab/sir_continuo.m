function sol = sir_continuo(pars, domain, opt)
% SIR_CONTINUO Modelo epidémico SIR en tiempo continuo (modelo de usuario).
%
%   sol = sir_continuo(pars, domain, opt)
%
%   dS/dt = -beta*S*I
%   dI/dt =  beta*S*I - gamma*I
%   dR/dt =  gamma*I
%
%   Factores, en este orden (primero parámetros, después condiciones iniciales):
%       pars(1) -> beta  : tasa de transmisión (1/(persona*día))
%       pars(2) -> gamma : tasa de recuperación (1/día)
%       pars(3) -> S0    : susceptibles en t0 (personas)
%       pars(4) -> I0    : infectados en t0 (personas)
%   R(t0) = 0 no es un factor: es una constante del modelo.
%
%   domain = [t0, tf] en días.
%   opt.dt (opcional) = separación de la rejilla de salida, en días (0.1 por defecto).
%
%   Salida: sol.x (1 x Nt) rejilla t0:dt:tf
%           sol.y (3 x Nt) filas = salidas [S; I; R]

beta  = pars(1);
gamma = pars(2);
y0    = [pars(3); pars(4); 0];

dt = 0.1;
if isstruct(opt) && isfield(opt, 'dt')
    dt = opt.dt;
end

rhs = @(t, y) [ -beta*y(1)*y(2);
                 beta*y(1)*y(2) - gamma*y(2);
                 gamma*y(2) ];

ode = ode45(rhs, domain, y0, odeset('RelTol', 1e-6, 'AbsTol', 1e-8));

% ode45 avanza con pasos adaptativos y desiguales. deval evalúa su
% interpolante de alto orden en una rejilla fina y regular; así la
% interpolación lineal que hace gsua_eval después es casi exacta.
sol.x = domain(1):dt:domain(2);
sol.y = deval(ode, sol.x);      % 3 x Nt
end
