function solution = sir_sencillo(factores, domain, input)
%SIR_SENCILLO Modelo SIR discreto (SIRS) con probabilidad de infeccion 1-exp(-beta*I/N).
%
%   solution = sir_sencillo(factores, domain, input)
%
%   Modelo de ejemplo usado por el curso de Modelacion Experimental (EAFIT) y por
%   las pruebas de identificabilidad practica. Poblacion cerrada, N = S0+I0+R0
%   constante. alpha es la perdida de inmunidad (R -> S).
%
%   factores = [S0 I0 R0 beta gamma alpha]
%     S0,I0,R0 <-- condiciones iniciales (susceptibles, infectados, recuperados)
%     beta     <-- transmisibilidad
%     gamma    <-- probabilidad de recuperacion por paso (I -> R)
%     alpha    <-- probabilidad de perdida de inmunidad por paso (R -> S)
%   domain  = [t0 tf]  (se simula en pasos enteros t0, t0+1, ..., tf)
%   input   <-- ignorado (firma Kind 5 de gsua_userdefined)
%
%   solution.y <-- 3 x numel(t0:tf) con las filas [S; I; R]
%   solution.x <-- t0:tf
%
%   Ecuaciones de actualizacion (P = 1 - exp(-beta*I/N)):
%     S(k+1) = S(k)*(1-P) + alpha*R(k)
%     I(k+1) = S(k)*P     + I(k)*(1-gamma)
%     R(k+1) = I(k)*gamma + R(k)*(1-alpha)

simulation_steps = domain(1):(domain(2)-1);
numSteps = numel(simulation_steps) + 1;
estados = zeros(3, numSteps);

%Asignamos las condiciones iniciales
estados(1:3,1) = factores(1:3);

%Asignamos parámetros
beta = factores(4);
gamma = factores(5);
alpha = factores(6);

%simulamos

for k = simulation_steps
    N = sum(estados(:,k)); %población total
    % Probabilidad de que un susceptible se infecte: 1 - exp(-beta*I/N).
    % (El original usaba beta*S*I/N, que cuenta S dos veces: con S del orden de
    % miles casi toda la población se infecta en un solo paso.)
    prob_inf = 1 - exp(-beta*estados(2,k)/N);

    estados(1,k+1) = estados(1,k)*(1-prob_inf) + alpha*estados(3,k);
    estados(2,k+1) = estados(1,k)*prob_inf + estados(2,k)*(1-gamma);
    estados(3,k+1) = estados(2,k)*gamma + estados(3,k)*(1-alpha);
end

solution.y = estados;
solution.x = domain(1):domain(2);

end
