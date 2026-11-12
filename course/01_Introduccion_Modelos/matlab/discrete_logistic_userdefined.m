function sol = discrete_logistic_userdefined(params, domain, ~)
% DISCRETE_LOGISTIC_USERDEFINED Modelo de crecimiento logístico en tiempo discreto.
%
%   sol = discrete_logistic_userdefined(params, domain, opt)
%
%   Ecuación en diferencias:
%       N(t+1) = N(t) + r * N(t) * (1 - N(t) / K)
%
%   Factores (primero la condición inicial, después los parámetros):
%       params(1) -> N0: Población inicial en t = t0
%       params(2) -> r : Tasa intrínseca de crecimiento (1/tiempo)
%       params(3) -> K : Capacidad de carga poblacional
%
%   Dominio:
%       domain = [t0, tf] (pasos discretos de tiempo)
%
%   Retorna:
%       Estructura sol con campos sol.x (puntos de tiempo) y sol.y (población).

% gsua reinserta los factores fijos en su posición, así que el modelo
% siempre recibe el vector completo (N0 incluido) en el orden declarado.
N0 = params(1);
r  = params(2);
K  = params(3);

t0 = floor(domain(1));
tf = ceil(domain(2));
time_steps = t0:tf;
num_steps = length(time_steps);

N = zeros(1, num_steps);
N(1) = N0;

for t = 1:(num_steps - 1)
    N(t+1) = N(t) + r * N(t) * (1 - N(t) / K);
    if N(t+1) < 0
        N(t+1) = 0;
    end
end

sol.x = time_steps;
sol.y = N;
end
