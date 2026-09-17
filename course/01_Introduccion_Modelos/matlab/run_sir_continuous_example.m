% RUN_SIR_CONTINUOUS_EXAMPLE Ejemplo completo de Modelo Epidémico SIR Continuo
% Curso guiado de GSUA-CSB
%
% Este script muestra el flujo completo en MATLAB:
% 1. Definición del modelo continuo en EDO ('sir_continuous_userdefined')
% 2. Construcción de la TABLA DE CONTROL (T) mediante gsua_dataprep
% 3. Simulación de múltiples salidas (S, I, R) y filtrado mediante gsua_eval

clc; clear; close all;

% Asegurar la inclusión de funciones del toolbox
toolbox_path = fullfile('C:', 'Users', 'Lanre', 'MATLAB Drive', 'Toolbox', 'Functions');
if exist(toolbox_path, 'dir'), addpath(toolbox_path); end

fprintf('=== EJEMPLO 2: MODELO EPIDÉMICO SIR CONTINUO ===\n\n');

%% PASO 1: Especificar Rangos, Dominio, Nombres de Factores (Parámetros + C.I.) y Salidas
ranges = [
    0.0001, 0.0008;  % beta  (tasa de transmisión)
    0.05,   0.20;    % gamma (tasa de recuperación)
    900,    999;     % S0    (Susceptibles iniciales)
    1,      50       % I0    (Infectados iniciales)
];
domain = [0, 80]; % 80 días
factor_names = {'beta', 'gamma', 'S0', 'I0'};
output_names = {'Susceptibles', 'Infectados', 'Recuperados'};

%% PASO 2: Crear la Tabla de Control del Sistema (T) con gsua_dataprep
[T, solver] = gsua_dataprep('sir_continuous_userdefined', ranges, ...
    'domain', domain, ...
    'names', factor_names, ...
    'out_names', output_names);

fprintf('1. Tabla resumen T creada por gsua_dataprep:\n');
disp(T);

%% PASO 3: Simular y Evaluar Todas las Salidas con gsua_eval
xdata = linspace(0, 80, 200);

% T tiene las 3 salidas activas [1 2 3] por defecto
y_all = gsua_eval(T.Nominal, T, xdata);

if size(y_all, 1) == 3
    S = y_all(1,:); I = y_all(2,:); R = y_all(3,:);
else
    S = y_all(:,1)'; I = y_all(:,2)'; R = y_all(:,3)';
end

%% PASO 4: Simular Únicamente la Salida de Infectados (Salida 2)
% Modificar la propiedad CustomProperties.output en T para filtrar salidas
T.Properties.CustomProperties.output = 2;
y_infectados_only = gsua_eval(T.Nominal, T, xdata);

%% PASO 5: Graficar los Resultados de gsua_eval
figure('Name', 'Modelo SIR Continuo - gsua_eval', 'Position', [100 100 900 400]);

subplot(1, 2, 1);
plot(xdata, S, 'b-', 'LineWidth', 2, 'DisplayName', 'S(t) Susceptibles');
hold on;
plot(xdata, I, 'r-', 'LineWidth', 2, 'DisplayName', 'I(t) Infectados');
plot(xdata, R, 'g-', 'LineWidth', 2, 'DisplayName', 'R(t) Recuperados');
xlabel('Tiempo (Días)');
ylabel('Población');
title('Todas las Salidas Evaluadas con gsua\_eval');
grid on; legend('Location', 'east');

subplot(1, 2, 2);
plot(xdata, y_infectados_only, 'r-', 'LineWidth', 2, 'DisplayName', 'Infectados I(t)');
xlabel('Tiempo (Días)');
ylabel('Infectados I(t)');
title('Salida Única de Infectados (output = 2)');
grid on; legend('Location', 'northeast');

fprintf('¡Ejemplo de Modelo SIR Continuo completado exitosamente!\n');
