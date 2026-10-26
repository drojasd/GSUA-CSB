function ok = verificar_modelo(nombre, pars, domain, out_names, opt)
% VERIFICAR_MODELO Revisa que un modelo de usuario cumpla el contrato de GSUA-CSB.
%
%   ok = verificar_modelo(nombre, pars, domain, out_names)
%   ok = verificar_modelo(nombre, pars, domain, out_names, opt)
%
%   Llama al modelo UNA vez, a mano, antes de entregárselo a gsua_dataprep,
%   y comprueba lo que la toolbox da por hecho:
%     1. el archivo nombre.m está en el path y declara 3 entradas;
%     2. devuelve una estructura con campos x e y;
%     3. sol.y tiene una fila por salida (numel(out_names)) y una columna
%        por punto de sol.x;
%     4. sol.x cubre el dominio completo [t0, tf];
%     5. no hay NaN ni Inf.
%
%   Es una ayuda del curso, no una función de la toolbox.

if nargin < 5, opt = struct(); end
ok = true;
fprintf('Revisando %s ...\n', nombre);

if exist(nombre, 'file') ~= 2
    fprintf('  FALLA  %s.m no está en el path.\n', nombre); ok = false; return
end
n_in = nargin(nombre);
if n_in ~= 3
    fprintf('  FALLA  declara %d entradas; gsua_eval lo llama con 3: (pars, domain, opt).\n', abs(n_in));
    ok = false; return
end
fprintf('  ok     firma sol = %s(pars, domain, opt)\n', nombre);

sol = feval(nombre, pars, domain, opt);
if ~isstruct(sol) || ~isfield(sol, 'x') || ~isfield(sol, 'y')
    fprintf('  FALLA  la salida debe ser una estructura con campos x e y.\n'); ok = false; return
end

nx = numel(sol.x);
[nf, nc] = size(sol.y);
ns = numel(out_names);
fprintf('  info   sol.x: %d puntos;  sol.y: %d x %d;  out_names: %d\n', nx, nf, nc, ns);
if nc == nx && nf == ns
    fprintf('  ok     sol.y es (salidas x tiempo)\n');
elseif nf == nx && nc == ns
    fprintf('  aviso  sol.y viene (tiempo x salidas); la toolbox la transpone, pero mejor devuélvela (salidas x tiempo)\n');
else
    fprintf('  FALLA  sol.y debería ser %d x %d (una fila por salida, una columna por tiempo)\n', ns, nx);
    ok = false;
end

if abs(sol.x(1) - domain(1)) > 1e-9 || abs(sol.x(end) - domain(end)) > 1e-9
    fprintf('  FALLA  sol.x va de %g a %g, pero el dominio es [%g, %g]\n', ...
        sol.x(1), sol.x(end), domain(1), domain(end));
    ok = false;
else
    fprintf('  ok     sol.x cubre el dominio [%g, %g]\n', domain(1), domain(end));
end

if any(~isfinite(sol.y(:)))
    fprintf('  FALLA  sol.y tiene NaN o Inf\n'); ok = false;
else
    fprintf('  ok     sin NaN ni Inf\n');
end

if ok, fprintf('  => cumple el contrato\n\n'); else, fprintf('  => revisa los puntos marcados\n\n'); end
end
