function ok = run_checks()
%RUN_CHECKS Verificación ejecutable del entorno y de las convenciones GSUA-CSB.
%   ok = run_checks()  devuelve true si todo pasa. Imprime un resumen.
%   Comprueba: (1) la toolbox y sus funciones clave están en el path;
%   (2) el invariante de factores fijos (el modelo recibe el vector completo);
%   (3) el desempaque CI-primero de un modelo del curso vía gsua_eval.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % .../Toolbox
addpath(fullfile(root,'Functions'));
addpath(fullfile(root,'course','00_Escribir_Modelos','matlab'));
ok = true;
fprintf('==== run_checks GSUA-CSB ====\n');

% (1) funciones clave en el path
need = {'gsua_dataprep','gsua_eval','gsua_ua','gsua_sa','gsua_pe','gsua_bootstrap', ...
        'gsua_ia','gsua_clustersummary','gsua_fullpars','gsua_dmatrix'};
for i = 1:numel(need)
    present = ~isempty(which(need{i}));
    ok = ok & present;
    fprintf('  [%s] %s\n', tick(present), need{i});
end

% (2) invariante: un factor fijo se reinserta en su posición (el modelo recibe Np completo)
tmp = fullfile(tempdir,'pb_probe_model.m');
fid = fopen(tmp,'w');
fprintf(fid,'function sol = pb_probe_model(params, domain, ~)\n');
fprintf(fid,'  assignin(''base'',''PB_LEN'', numel(params));\n');
fprintf(fid,'  t = domain(1):domain(2); sol.x = t; sol.y = params(1)*ones(1,numel(t));\n');
fprintf(fid,'end\n'); fclose(fid);
addpath(tempdir); rehash;
[Tp,~] = gsua_dataprep('pb_probe_model',[1 1;2 4;5 5],'domain',[0 3], ...
    'names',{'a','b','c'},'out_names',{'out'});
evalc('gsua_eval(Tp.Nominal, Tp, 0:3, [], false, false, false)');
full_vec = (evalin('base','PB_LEN') == 3);
ok = ok & full_vec;
fprintf('  [%s] factor fijo reinsertado (modelo recibe Np=3 con 1 libre)\n', tick(full_vec));

% (3) desempaque CI-primero: sir_continuo con [S0,I0,beta,gamma]
[Ts,~] = gsua_dataprep('sir_continuo',[900 999;1 50;1e-4 8e-4;0.05 0.20],'domain',[0 80], ...
    'names',{'S0','I0','beta','gamma'},'out_names',{'S','I','R'},'opt',struct('dt',0.5));
y = gsua_eval(Ts.Nominal, Ts, 0:80, [], false, false, false);   % 1 x Nt x 3
S0n = Ts.Nominal(1); I0n = Ts.Nominal(2);
ci_ok = abs(y(1,1,1)-S0n) < 1e-6 && abs(y(1,1,2)-I0n) < 1e-6;
ok = ok & ci_ok;
fprintf('  [%s] CI-primero: y(t0)=[%.1f %.1f] == [S0 I0]=[%.1f %.1f]\n', ...
    tick(ci_ok), y(1,1,1), y(1,1,2), S0n, I0n);

fprintf('==== RESULTADO: %s ====\n', ternary(ok,'TODO OK','HAY FALLAS'));
end

function s = tick(b),   if b, s='OK'; else, s='X '; end, end
function s = ternary(b,a,c), if b, s=a; else, s=c; end, end
