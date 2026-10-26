function sol = sir_discreto_solucion(pars, domain, ~)
% SIR_DISCRETO_SOLUCION Solución del ejercicio del módulo 0 (SIR discreto).
%
%   Ábrela después de intentar SIR_DISCRETO_PLANTILLA por tu cuenta.

beta  = pars(1);
gamma = pars(2);

t = domain(1):domain(2);
Y = zeros(3, numel(t));
Y(:, 1) = [pars(3); pars(4); 0];

for k = 1:numel(t)-1
    S = Y(1, k);  I = Y(2, k);  R = Y(3, k);
    nuevos     = beta*S*I;       % contagios del día
    recuperados = gamma*I;       % recuperaciones del día
    Y(:, k+1) = [S - nuevos; I + nuevos - recuperados; R + recuperados];
end

sol.x = t;
sol.y = Y;
end
