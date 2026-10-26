function sol = sir_discreto_plantilla(pars, domain, ~)
% SIR_DISCRETO_PLANTILLA Ejercicio del módulo 0: un SIR en tiempo discreto.
%
%   Copia este archivo como mi_sir_discreto.m (y cambia también el nombre
%   de la función en la primera línea) y completa los TODO.
%
%   Con un paso de un día, el SIR discreto es:
%       S(t+1) = S(t) - beta*S(t)*I(t)
%       I(t+1) = I(t) + beta*S(t)*I(t) - gamma*I(t)
%       R(t+1) = R(t) + gamma*I(t)
%
%   Factores (primero parámetros, después condiciones iniciales):
%       pars(1) -> beta  (1/(persona*día))
%       pars(2) -> gamma (1/día)
%       pars(3) -> S0    (personas)
%       pars(4) -> I0    (personas)
%   R0 = 0 al inicio.
%
%   Salida: sol.x (1 x Nt) = domain(1):domain(2)
%           sol.y (3 x Nt) = [S; I; R]

beta  = pars(1);
gamma = pars(2);

t = domain(1):domain(2);
Y = zeros(3, numel(t));          % 3 salidas -> 3 filas
Y(:, 1) = [pars(3); pars(4); 0];

for k = 1:numel(t)-1
    S = Y(1, k);  I = Y(2, k);  R = Y(3, k);
    % TODO: escribe aquí las tres ecuaciones de actualización.
    %       Mientras tanto la plantilla deja S, I y R constantes.
    Y(:, k+1) = [S; I; R];
end

sol.x = t;
sol.y = Y;
end
