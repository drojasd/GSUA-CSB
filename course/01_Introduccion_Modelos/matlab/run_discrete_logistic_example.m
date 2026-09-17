% RUN_DISCRETE_LOGISTIC_EXAMPLE Ejemplo completo de Modelo Logístico Discreto
% Curso guiado de GSUA-CSB
%
% Este script muestra el flujo completo en MATLAB:
% 1. Definición del modelo en tiempo discreto ('discrete_logistic_userdefined')
% 2. Construcción de la TABLA DE CONTROL (T) mediante gsua_dataprep
% 3. Simulación y evaluación usando gsua_eval

clc; clear; close all;

% Asegurar la inclusión de funciones del toolbox
toolbox_path = fullfile('C:', 'Users', 'Lanre', 'MATLAB Drive', 'Toolbox', 'Functions');
if exist(toolbox_path, 'dir'), addpath(toolbox_path); end

fprintf('=== EJEMPLO 1: MODELO LOGÍSTICO DISCRETO ===\n\n');

%% PASO 1: Especificar Rangos, Dominio, Nombres de Factores y Salidas
ranges = [
    0.05, 0.50;  % r  (tasa de crecimiento)
    100,  1000;  % K  (capacidad de carga)
    5,    50     % N0 (población inicial)
];
domain = [0, 50]; % t0 = 0, tf = 50
factor_names = {'r', 'K', 'N0'};
output_names = {'Poblacion'};

%% PASO 2: Crear la Tabla de Control del Sistema (T) con gsua_dataprep
[T, solver] = gsua_dataprep('discrete_logistic_userdefined', ranges, ...
    'domain', domain, ...
    'names', factor_names, ...
    'out_names', output_names);

fprintf('1. Tabla resumen T creada por gsua_dataprep:\n');
disp(T);

fprintf('2. Propiedades CustomProperties en T que controlan a gsua_csb:\n');
fprintf('   - Solver: %s\n', func2str(T.Properties.CustomProperties.Solver));
fprintf('   - Kind: %d\n', T.Properties.CustomProperties.Kind);
fprintf('   - Domain: [%d, %d]\n', T.Properties.CustomProperties.Domain(1), T.Properties.CustomProperties.Domain(2));
fprintf('   - Output Activo: %s\n\n', num2str(T.Properties.CustomProperties.output));

%% PASO 3: Simular y Evaluar el Modelo con gsua_eval
xdata = 0:50; % Puntos de tiempo para la evaluación

% Evaluación en el punto nominal [r_mid, K_mid, N0_mid]
y_nominal = gsua_eval(T.Nominal, T, xdata);

% Evaluación en un punto de prueba alternativo
test_params = [0.40; 800; 10]; % r=0.40, K=800, N0=10
y_test = gsua_eval(test_params, T, xdata);

%% PASO 4: Graficar los Resultados de gsua_eval
figure('Name', 'Modelo Logístico Discreto - gsua_eval');
plot(xdata, y_nominal, 'b-o', 'LineWidth', 1.5, 'DisplayName', 'Nominal (r=0.275, K=550)');
hold on;
plot(xdata, y_test, 'r-s', 'LineWidth', 1.5, 'DisplayName', 'Prueba (r=0.40, K=800)');
xlabel('Paso de Tiempo (t)');
ylabel('Población N(t)');
title('Simulación con gsua\_eval para Modelo Logístico Discreto');
grid on;
legend('Location', 'southeast');

fprintf('¡Ejemplo de Modelo Logístico Discreto completado exitosamente!\n');
