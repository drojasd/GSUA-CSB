%% RECETA: estimación (PE) + identificabilidad práctica (IA)
% Datos sintéticos (verdad + ruido) -> gsua_pe -> bootstrap -> gsua_ia -> clustersummary.
% Fija I0 (rango degenerado) y estima S0, beta, gamma. Ver ../references/lineamientos.md §5-6.
clc; clear; close all;

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % .../Toolbox
addpath(fullfile(root,'Functions'));
addpath(fullfile(root,'course','01_Introduccion_Modelos','matlab'));

% --- Declarar (CI-primero; I0 fijo en 5) ------------------------------------
ranges = [900 999;  5 5;  1e-4 8e-4;  0.05 0.20];   % S0(libre), I0(FIJO), beta, gamma
[T, ~] = gsua_dataprep('sir_continuous_userdefined', ranges, 'domain', [0 80], ...
    'names', {'S0','I0','beta','gamma'}, 'out_names', {'S','I','R'});
T.Properties.CustomProperties.output = 2;            % solo Infectados
names = T.Properties.RowNames;                        % {'S0','beta','gamma'}
xdata = linspace(0,80,200);
opt   = optimoptions('fmincon','Display','off');

% --- Datos "observados" sintéticos (verdad conocida + 5% ruido) -------------
% ydata va como SALIDAS x TIEMPO (fila, para una sola salida); así la devuelve
% gsua_deval. Con el solver por defecto (fmincon) un ydata en COLUMNA dispara un
% error críptico 'badsubscript' (solo lsqc se autocorrige). Ver ../references/errores.md.
true_free = [950; 5e-4; 0.12];                        % S0, beta, gamma (orden de los libres)
y_clean   = squeeze(gsua_eval(true_free, T, xdata, [], false, false, false));   % 1 x Nt (fila)
rng(3); ydata = y_clean .* (1 + 0.05*randn(size(y_clean)));

% --- Estimación multicomienzo (solver por defecto: fmincon; margin crudo) ----
[T, res] = gsua_pe(T, xdata, ydata, 'N', 10, 'solver', 'fmincon', 'margin', 0.1, ...
    'opt', opt, 'save', false, 'timer', false);

% --- Identificabilidad por ensamble -----------------------------------------
Est = gsua_bootstrap(T, xdata, ydata, T.Estfmincon, res, 'B', 20, 'solver', 'fmincon', ...
    'N', 3, 'noise', 'gaussian', 'opt', opt, 'seed', 7, 'timer', false);
[T, ci] = gsua_ia(T, Est, false, false, true, false);
S = gsua_clustersummary(T, ci);
disp(S(:, {'Cluster','Parameter','PctWidthRel','Verdict'}));   % lee la banda, no el IC de la mediana

% La combinación identificable: R0 = beta*S0/gamma  (libres {S0,beta,gamma})
R0 = (Est(2,:).*Est(1,:))./Est(3,:);
fprintf('R0 = %.2f ± %.2f  (verdad %.2f)\n', mean(R0), std(R0), (true_free(2)*true_free(1))/true_free(3));
