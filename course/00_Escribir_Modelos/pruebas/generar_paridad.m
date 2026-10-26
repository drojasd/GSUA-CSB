%GENERAR_PARIDAD  Archivos de paridad MATLAB -> Python de los modelos del módulo 0.
%   Para cada modelo de la lista escribe en esta carpeta:
%     pars_<modelo>.csv     20 vectores de factores (hipercubo latino, semilla fija),
%                           una fila por vector, columnas en el orden de los factores;
%     x_<modelo>.csv        la rejilla sol.x que devuelve MATLAB (una fila);
%     salidas_<modelo>.csv  una fila por vector: sol.y aplanado por filas
%                           (toda la fila 1 de sol.y, luego la fila 2, ...).
%   paridad.py corre los mismos vectores en Python sobre la misma rejilla y compara.

aqui = fileparts(mfilename('fullpath'));
if isempty(aqui), aqui = pwd; end
mod0 = fileparts(aqui);
addpath(fullfile(mod0, 'matlab'), fullfile(mod0, 'matlab', 'solucion'));

rangos_log = [0.05 0.50; 100 1000; 5 50];
rangos_sir = [1e-4 8e-4; 0.05 0.20; 900 999; 1 50];
modelos = {                                          % nombre, rangos, domain
    'logistica_discreta',     rangos_log, [0 50]
    'sir_continuo',           rangos_sir, [0 80]
    'sir_discreto_plantilla', rangos_sir, [0 80]
    'sir_discreto_solucion',  rangos_sir, [0 80] };
nVec = 20;

for k = 1:size(modelos, 1)
    [nombre, R, domain] = modelos{k, :};
    rng(20261026 + k, 'twister');
    Np = size(R, 1);
    U = zeros(nVec, Np);
    for j = 1:Np                                     % hipercubo latino simple
        U(:, j) = (randperm(nVec).' - rand(nVec, 1))/nVec;
    end
    P = R(:, 1).' + U.*(R(:, 2) - R(:, 1)).';
    salidas = [];
    for i = 1:nVec
        sol = feval(nombre, P(i, :).', domain, struct());
        fila = reshape(sol.y.', 1, []);
        if isempty(salidas), salidas = zeros(nVec, numel(fila)); end
        salidas(i, :) = fila;
    end
    escribir(fullfile(aqui, ['pars_' nombre '.csv']), P);
    escribir(fullfile(aqui, ['x_' nombre '.csv']), sol.x);
    escribir(fullfile(aqui, ['salidas_' nombre '.csv']), salidas);
    fprintf('  %-24s %d vectores x %d factores -> %d columnas de salida\n', ...
        nombre, nVec, Np, size(salidas, 2));
end

function escribir(ruta, M)
% CSV sin encabezado, 17 cifras significativas (ida y vuelta exacta de double).
fid = fopen(ruta, 'w');
formato = [repmat('%.17g,', 1, size(M, 2) - 1) '%.17g\n'];
fprintf(fid, formato, M.');
fclose(fid);
end
