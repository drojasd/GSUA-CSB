%% RECETA: UA + SA sobre un SIR continuo (CI-primero)
% Flujo mínimo declarar -> UA (ensamble) -> SA (con chequeo de confiabilidad).
% Usa el modelo sir_continuo del curso. Ver ../references/{flujos,lineamientos}.md
clc; clear; close all;

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % .../Toolbox
addpath(fullfile(root,'Functions'));
addpath(fullfile(root,'course','00_Escribir_Modelos','matlab'));

% --- Declarar (CI-primero: fila i de rangos = pars(i) = names{i}) -------------
ranges = [900 999;  1 50;  1e-4 8e-4;  0.05 0.20];   % S0, I0, beta, gamma
[T, ~] = gsua_dataprep('sir_continuo', ranges, 'domain', [0 80], ...
    'names', {'S0','I0','beta','gamma'}, 'out_names', {'S','I','R'}, 'opt', struct('dt',0.1));
T.Properties.CustomProperties.output = 2;            % analizar Infectados
xdata = linspace(0,80,161);

% --- UA: ensamble de trayectorias -------------------------------------------
M = gsua_dmatrix(T, 200);                            % 200 × Np
Yua = gsua_eval(M', T, xdata, [], false, false, false);
fprintf('UA: ensamble size = %s\n', mat2str(size(Yua)));

% --- SA: índices de sensibilidad + confiabilidad ----------------------------
Msa = gsua_dmatrix(T, 500);
[T, ~, ~] = gsua_sa(Msa, T, 'parallel', false, 'verbose', false);
Si = T.Si;                                           % primer orden por factor
reliability = sum(Si) / sum(abs(Si));                % -> 1 cuando N es suficiente
fprintf('SA: sum(Si)=%.3f  sum(|Si|)=%.3f  confiabilidad=%.3f (sube N si <~0.95)\n', ...
    sum(Si), sum(abs(Si)), reliability);
% NOTA: no existe umbral publicado (p.ej. sum(Si)>0.65). sum(Si)->1 = aditividad;
%       sum(Si)/sum(|Si|)->1 = N suficiente. Ver lineamientos.md §4.
