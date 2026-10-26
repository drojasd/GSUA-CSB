%% Módulo 00: Cómo se escribe un modelo para GSUA-CSB
% Curso guiado de GSUA-CSB
%
% Todo el curso (muestreo, incertidumbre, sensibilidad, estimación) trabaja
% sobre un *modelo de usuario*: una función de MATLAB que tú escribes. La
% toolbox no ve tus ecuaciones; solo llama a esa función con un vector de
% factores y lee lo que devuelve. Este módulo enseña ese contrato paso a paso:
%
% # la firma |sol = modelo(pars, domain, opt)| y la forma de |sol.x| y |sol.y|;
% # el orden de los factores: primero los parámetros, después las condiciones iniciales;
% # un modelo discreto (un bucle) y uno continuo (|ode45| + |deval| en una rejilla);
% # cómo declararlo con |gsua_dataprep| (rangos, |names|, |out_names|), qué es la tabla |T| y cómo se evalúa con |gsua_eval|;
% # los errores típicos (firma, orientación, rejilla, unidades);
% # un ejercicio: programar un SIR discreto desde cero.
%
% El código completo de cada modelo aparece en la página tal como está en su archivo.

%% 1. Configuración del entorno
clc; clear; close all;

% --- Localización portable de la toolbox GSUA-CSB ---------------------------
% Si la toolbox ya está instalada (.mltbx), sus funciones están en el path.
% Si no, se busca dentro del repositorio clonado y en ubicaciones habituales.
if isempty(which('gsua_dataprep'))
    candidatos = {
        fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', 'Functions')  % dentro del repo clonado
        fullfile(getenv('USERPROFILE'), 'MATLAB Drive', 'Toolbox', 'Functions')
        fullfile(getenv('USERPROFILE'), 'Documents', 'MATLAB', 'GSUA-CSB', 'Functions')
    };
    for k = 1:numel(candidatos)
        if isfolder(candidatos{k}), addpath(candidatos{k}); break; end
    end
    if isempty(which('gsua_dataprep'))
        error(['No se encontró la toolbox GSUA-CSB en el path de MATLAB. ' ...
               'Instala GSUA-CSB.mltbx o añade su carpeta Functions con addpath(...).']);
    end
end

% Los modelos de este módulo viven junto a este script (y en errores/).
aqui = fileparts(mfilename('fullpath'));
addpath(aqui, fullfile(aqui, 'errores'));

%% 2. El contrato de un modelo de usuario
% Un modelo de usuario es un archivo |.m| propio (por ejemplo |mi_modelo.m|)
% que contiene una función con esta firma:
%
%   sol = mi_modelo(pars, domain, opt)
%
% *Entradas*
%
% * |pars|: vector columna con los factores. La posición manda: |pars(1)| es
%   el factor de la fila 1 de los rangos, |pars(2)| el de la fila 2, etc.
%   Convención del curso: *primero los parámetros, después las condiciones
%   iniciales*.
% * |domain|: |[t0, tf]|, el intervalo de simulación. Es el mismo que
%   declaras con |'domain'| en |gsua_dataprep|; la toolbox no lo cambia.
% * |opt|: cualquier cosa que pases con |'opt'| en |gsua_dataprep|
%   (por defecto, una estructura que puedes ignorar con |~|).
%
% *Salida*: una estructura con dos campos
%
% * |sol.x|: vector de tiempos, 1 x Nt, de |t0| a |tf|;
% * |sol.y|: matriz *salidas x tiempo* (nSalidas x Nt): una fila por salida,
%   en el mismo orden de |out_names|.
%
% Después, |gsua_eval| interpola linealmente |sol.y| en los tiempos que le
% pidas (|xdata|). Por eso |sol.x| debe cubrir todo el dominio y ser lo
% bastante fino.

%% 3. Paso 1 — un modelo discreto: la logística
% $N(t+1) = N(t) + r\,N(t)\,(1 - N(t)/K)$, con tres factores: |r|, |K| (parámetros)
% y |N0| (condición inicial). Un bucle avanza un paso de tiempo por iteración.
% Este es el archivo |logistica_discreta.m| completo:
%
% <include>logistica_discreta.m</include>
%
% Antes de dárselo a la toolbox, llámalo tú a mano y mira las formas:

sol = logistica_discreta([0.3; 500; 10], [0 50], []);
fprintf('sol.x: %d x %d   (tiempos t0:tf)\n', size(sol.x));
fprintf('sol.y: %d x %d   (1 salida x %d tiempos)\n', size(sol.y), numel(sol.x));
fprintf('N(0) = %.1f,  N(50) = %.1f\n\n', sol.y(1), sol.y(end));

%% 4. Paso 2 — declararlo con gsua_dataprep: rangos, names y out_names
% |gsua_dataprep| recibe el *nombre* del archivo (texto, no un handle) y una
% matriz de rangos |Np x 2| con |[mínimo, máximo]| por fila. La fila i
% corresponde a |pars(i)| y a |names{i}|: las tres listas deben ir en el
% mismo orden. |out_names| nombra las filas de |sol.y|, también en orden.

rangos_log = [
    0.05,  0.50;   % r   (parámetro)
    100,   1000;   % K   (parámetro)
    5,     50      % N0  (condición inicial, al final)
];

[T_log, ~] = gsua_dataprep('logistica_discreta', rangos_log, ...
    'domain',    [0 50], ...
    'names',     {'r', 'K', 'N0'}, ...
    'out_names', {'Poblacion'});

%%
% El resultado |T| es una *tabla*: una fila por factor libre, con su rango
% y su valor nominal (por defecto, el punto medio del rango).
disp(T_log)

%%
% Lo demás que necesita la toolbox va guardado en las propiedades de la tabla.
% Todas las funciones del curso (|gsua_dmatrix|, |gsua_ua|, |gsua_sa|,
% |gsua_pe|...) leen el modelo desde aquí:
P = T_log.Properties.CustomProperties;
fprintf('Solver : %s\n', func2str(P.Solver));
fprintf('Kind   : %d   (5 = función de usuario con dominio)\n', P.Kind);
fprintf('Domain : [%g %g]\n', P.Domain);
fprintf('Vars   : %s   (out_names)\n', strjoin(P.Vars, ', '));
fprintf('output : %s   (salidas activas)\n\n', mat2str(P.output));

%% 5. Paso 3 — evaluar con gsua_eval
% |y = gsua_eval(par, T, xdata, ydata, parallel, show, verbose)|.
% |par| es *Np x N*: una columna por simulación. Los cuatro últimos argumentos
% son opcionales; aquí los ponemos para no graficar ni imprimir progreso.
% La salida tiene forma *N x numel(xdata) x nSalidasActivas*.

xdata = 0:50;
y_nom = gsua_eval(T_log.Nominal, T_log, xdata, [], false, false, false);
fprintf('Una simulación:  size(y) = %s\n', mat2str(size(y_nom)));

par3 = [T_log.Nominal, [0.45; 800; 10], [0.10; 300; 40]];   % 3 columnas = 3 simulaciones
y3 = gsua_eval(par3, T_log, xdata, [], false, false, false);
fprintf('Tres simulaciones: size(y) = %s\n\n', mat2str(size(y3)));

figure('Name', 'Logística discreta evaluada con gsua_eval');
plot(xdata, y3', 'LineWidth', 1.5);
xlabel('Paso de tiempo t'); ylabel('Población N(t)'); grid on;
legend('nominal', 'r=0.45, K=800, N0=10', 'r=0.10, K=300, N0=40', 'Location', 'southeast');
title('Cada columna de par es una simulación');

%% 6. Paso 4 — un modelo continuo: el SIR con ode45
% Un modelo continuo resuelve sus EDO *dentro* de la función. El patrón es:
%
% # integrar con |ode45| sobre |domain|;
% # evaluar la solución con |deval| en una rejilla fina y regular;
% # devolver esa rejilla en |sol.x| y los estados en |sol.y| (estados x tiempo).
%
% Factores: |beta|, |gamma| (parámetros), |S0|, |I0| (condiciones iniciales).
% Salidas: S, I y R, en ese orden de filas. Archivo |sir_continuo.m| completo:
%
% <include>sir_continuo.m</include>

sol = sir_continuo([4.5e-4; 0.125; 950; 25], [0 80], struct('dt', 0.5));
fprintf('sol.x: %d x %d,  sol.y: %d x %d  (S, I, R x tiempo)\n\n', size(sol.x), size(sol.y));

%%
% Al declararlo, |out_names| tiene tres nombres porque |sol.y| tiene tres filas.
% El argumento |'opt'| llega tal cual como tercer argumento del modelo.

rangos_sir = [
    1e-4,  8e-4;   % beta  (1/(persona*día))
    0.05,  0.20;   % gamma (1/día)
    900,   999;    % S0    (personas)
    1,     50      % I0    (personas)
];

[T_sir, ~] = gsua_dataprep('sir_continuo', rangos_sir, ...
    'domain',    [0 80], ...
    'names',     {'beta', 'gamma', 'S0', 'I0'}, ...
    'out_names', {'S', 'I', 'R'}, ...
    'opt',       struct('dt', 0.1));
disp(T_sir)

xs = linspace(0, 80, 161);
y_sir = gsua_eval(T_sir.Nominal, T_sir, xs, [], false, false, false);
fprintf('Tres salidas activas: size(y) = %s  (1 simulación x tiempos x 3 salidas)\n', mat2str(size(y_sir)));

%%
% Para quedarte con una sola salida (por ejemplo, la que tienes medida),
% cambia |output| en la tabla. Afecta a todas las funciones que reciban |T|.
T_I = T_sir;
T_I.Properties.CustomProperties.output = 2;     % fila 2 de sol.y = I
y_I = gsua_eval(T_I.Nominal, T_I, xs, [], false, false, false);
fprintf('Solo I: size(y) = %s,  pico de I = %.1f personas\n\n', mat2str(size(y_I)), max(y_I));

figure('Name', 'SIR continuo evaluado con gsua_eval');
plot(xs, squeeze(y_sir), 'LineWidth', 1.8);
xlabel('Tiempo (días)'); ylabel('Personas'); grid on;
legend('S', 'I', 'R', 'Location', 'east');
title('Las tres salidas: y(1,:,j) es la salida j');

%% 7. Fijar un factor: rango degenerado
% Si un factor no te interesa, dale un rango con mínimo = máximo. Desaparece
% de |T| (ya no se muestrea ni se estima), pero el modelo lo sigue recibiendo
% en su posición: |pars| siempre llega completo. Lo que cambia es |par| en
% |gsua_eval|, que solo lleva los factores *libres*.

rangos_fijo = rangos_sir;
rangos_fijo(4, :) = [10 10];                    % I0 fijo en 10
[T_fijo, ~] = gsua_dataprep('sir_continuo', rangos_fijo, 'domain', [0 80], ...
    'names', {'beta', 'gamma', 'S0', 'I0'}, 'out_names', {'S', 'I', 'R'});
disp(T_fijo)
y_fijo = gsua_eval(T_fijo.Nominal, T_fijo, xs, [], false, false, false);   % par: 3 x 1
fprintf('par con %d factores libres; el modelo recibió I0 = 10 en pars(4). I(0) = %.1f\n\n', ...
    height(T_fijo), y_fijo(1, 1, 2));

%% 8. Errores típicos
% *(a) Firma con dos entradas.* La guía de usuario muestra modelos
% |modelo(pars, domain)| y |gsua_dataprep| los acepta, pero |gsua_eval| llama
% siempre con tres argumentos. Declara el tercero aunque no lo uses (|~|).
[T_mal, ~] = gsua_dataprep('logistica_dos_entradas', rangos_log, 'domain', [0 50]);
try
    gsua_eval(T_mal.Nominal, T_mal, xdata, [], false, false, false);
catch err
    fprintf('logistica_dos_entradas -> error: %s\n\n', err.message);
end

%%
% *(b) out_names no coincide con las filas.* |out_names| no reordena nada:
% solo pone nombre a las filas de |sol.y|, en orden. Si declaras una sola
% salida llamada 'I' para un modelo cuya fila 1 es S, la toolbox analiza S
% con la etiqueta 'I' y no da ningún error.
[T_etiq, ~] = gsua_dataprep('sir_continuo', rangos_sir, 'domain', [0 80], ...
    'out_names', {'I'});
y_etiq = gsua_eval(T_etiq.Nominal, T_etiq, xs, [], false, false, false);
fprintf('Salida etiquetada ''I'': valor inicial %.1f (es S0, no I0)\n', y_etiq(1));
fprintf('Arreglo: out_names = {''S'',''I'',''R''} y output = 2.\n\n');

%%
% *(c) Orientación.* La convención es |sol.y| = *salidas x tiempo*. Desde la
% versión 1.13, si devuelves |tiempo x salidas| la toolbox lo transpone,
% porque reconoce el eje del tiempo por |numel(sol.x)|. Si las dos
% dimensiones coinciden (tantas salidas como tiempos) no puede saberlo y la
% deja como está: no dependas de la corrección, devuelve |sol.y| bien orientada.

%%
% *(d) Rejilla: fuera del dominio no hay datos.* |gsua_eval| interpola en
% |sol.x|; un |xdata| que se sale de |[t0, tf]| devuelve |NaN|.
y_fuera = gsua_eval(T_I.Nominal, T_I, 0:100, [], false, false, false);
fprintf('xdata = 0:100 con domain = [0 80]: %d valores NaN\n', sum(isnan(y_fuera)));

%%
% *(e) Rejilla: devolver ode45 tal cual.* La estructura de |ode45| tiene
% campos |x| e |y| con la orientación correcta, así que "funciona"; pero
% |sol.x| son los pasos adaptativos del integrador (pocos y desiguales), y
% |gsua_eval| interpola *linealmente* entre ellos. El modelo
% |errores/sir_continuo_crudo.m| hace exactamente eso:
%
% <include>errores/sir_continuo_crudo.m</include>

[T_crudo, ~] = gsua_dataprep('sir_continuo_crudo', rangos_sir, 'domain', [0 80], ...
    'out_names', {'S', 'I', 'R'}, 'output', 2);
xf = 0:0.05:80;
I_bien  = squeeze(gsua_eval(T_I.Nominal,     T_I,     xf, [], false, false, false));
I_crudo = squeeze(gsua_eval(T_crudo.Nominal, T_crudo, xf, [], false, false, false));
ode = sir_continuo_crudo(T_crudo.Nominal, [0 80]);
fprintf('ode45 crudo: %d pasos en 80 días. Error máximo en I: %.2f personas (pico %.1f)\n\n', ...
    numel(ode.x), max(abs(I_bien - I_crudo)), max(I_bien));

figure('Name', 'ode45 crudo frente a deval en rejilla');
plot(xf, I_bien, 'b-', 'LineWidth', 2); hold on;
plot(xf, I_crudo, 'r--', 'LineWidth', 1.2);
plot(ode.x, ode.y(2, :), 'ro', 'MarkerSize', 4);
xlim([0 40]); grid on;
xlabel('Tiempo (días)'); ylabel('Infectados I(t)');
legend('deval en rejilla de 0.1 días', 'ode45 crudo + interpolación lineal', ...
    'pasos de ode45', 'Location', 'northwest');
title('Devuelve la solución evaluada en una rejilla, no los pasos del integrador');

%%
% *(f) Unidades.* Rangos, condiciones iniciales y dominio deben estar en las
% mismas unidades. En |sir_continuo|, |beta| va en 1/(persona*día) porque
% el término es |beta*S*I| con S e I en personas. Si escribes S0 e I0 como
% fracciones de la población (0.95, 0.025) sin cambiar |beta|, el número
% reproductivo |R0 = beta*S0/gamma| se hunde y la epidemia desaparece;
% si el dominio está en semanas y |gamma| en 1/día, el tiempo de
% recuperación queda multiplicado por 7.
b = 4.5e-4; g = 0.125;
fprintf('R0 con S0 en personas  (950)  : %.2f\n', b*950/g);
fprintf('R0 con S0 en fracción  (0.95) : %.4f   <- mismas ecuaciones, otra escala\n\n', b*0.95/g);

%% 9. Comprobar tu modelo antes de usarlo
% |verificar_modelo.m| (ayuda del curso, no de la toolbox) llama a tu modelo
% una vez y revisa la firma, los campos, la orientación de |sol.y|, que
% |sol.x| cubra el dominio y que no haya |NaN|. Úsalo cada vez que escribas
% un modelo nuevo.
verificar_modelo('logistica_discreta', [0.3; 500; 10], [0 50], {'Poblacion'});
verificar_modelo('sir_continuo', T_sir.Nominal, [0 80], {'S', 'I', 'R'});
verificar_modelo('logistica_dos_entradas', [0.3; 500; 10], [0 50], {'Poblacion'});

%% 10. Ejercicio: programa un SIR discreto desde cero
% Con un paso de un día, el SIR discreto es
%
% $$S_{t+1} = S_t - \beta S_t I_t,\quad I_{t+1} = I_t + \beta S_t I_t - \gamma I_t,\quad R_{t+1} = R_t + \gamma I_t.$$
%
% # Copia |sir_discreto_plantilla.m| como |mi_sir_discreto.m| y cambia el
%   nombre de la función en su primera línea.
% # Completa el |TODO| del bucle con las tres ecuaciones.
% # Pásalo por |verificar_modelo| con |out_names = {'S','I','R'}|.
% # Decláralo con |gsua_dataprep| usando los mismos rangos, |names| y
%   |out_names| que |sir_continuo|, y |'domain', [0 80]|.
% # Evalúa el nominal con |gsua_eval| y grafica I(t) junto a la de |sir_continuo|.
%
% *Preguntas.* ¿Por qué no coinciden las dos curvas de I aunque las
% ecuaciones "son las mismas"? ¿Qué pasa con S si |beta*S0 > 1|? ¿Se conserva
% S + I + R en tu modelo?
%
% La plantilla (ya cumple el contrato; solo le falta la dinámica):
%
% <include>sir_discreto_plantilla.m</include>
%
% La solución está en |solucion/sir_discreto_solucion.m|. Ábrela cuando
% hayas terminado.

verificar_modelo('sir_discreto_plantilla', T_sir.Nominal, [0 80], {'S', 'I', 'R'});
[T_ej, ~] = gsua_dataprep('sir_discreto_plantilla', rangos_sir, 'domain', [0 80], ...
    'names', {'beta', 'gamma', 'S0', 'I0'}, 'out_names', {'S', 'I', 'R'}, 'output', 2);
I_ej = gsua_eval(T_ej.Nominal, T_ej, 0:80, [], false, false, false);
if all(abs(I_ej - I_ej(1)) < 1e-12)
    fprintf('La plantilla cumple el contrato, pero I(t) sigue constante: completa el TODO.\n');
else
    fprintf('I(t) ya cambia: pico de %.1f personas en el día %d.\n', max(I_ej), find(I_ej == max(I_ej), 1) - 1);
end

%% 11. Resumen: lista de chequeo
% * Archivo |.m| propio, en el path, con *tres* entradas: |sol = f(pars, domain, opt)|.
% * |pars(i)| = fila i de los rangos = |names{i}|: parámetros primero, condiciones iniciales después.
% * |sol.x| de 1 x Nt que cubre |[t0, tf]|; |sol.y| de nSalidas x Nt, filas en el orden de |out_names|.
% * Discreto: un bucle de |t0| a |tf|. Continuo: |ode45| sobre |domain| y |deval| en una rejilla fina.
% * |gsua_dataprep('f', rangos, 'domain', ..., 'names', ..., 'out_names', ...)| crea |T|.
% * |gsua_eval(par, T, xdata)| con |par| de Np x N devuelve N x numel(xdata) x nSalidas.
% * Una sola salida: |T.Properties.CustomProperties.output = k|. Un factor fijo: rango [v v].
%
% En el módulo 1 estos mismos dos modelos (|discrete_logistic_userdefined| y
% |sir_continuous_userdefined|) se exploran con ensambles de trayectorias.
