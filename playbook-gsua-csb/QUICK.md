# GSUA-CSB — guía rápida

La receta principal (lista para copiar) y las reglas que más se incumplen. El detalle
(flujos, errores, API) está en `references/` — empieza por `lineamientos.md`.

## Receta: estimación con identificabilidad práctica

Corre tal cual (requiere `gsua_*` y tu modelo en el path); el ejemplo usa el SIR del curso con
datos de conteo sintéticos. Reemplaza el modelo, los rangos y `ydata` por los tuyos.
**No copies la plantilla vieja del kit.**

```matlab
% Modelo Kind 5: sol = f(pars,domain,opt) con 3 entradas SIEMPRE;
% sol.x (1xNt), sol.y (salidas x tiempo). ode45 + deval en malla fina (no pasos crudos).
modelo = 'sir_continuous_userdefined';
ranges = [900 999;  5 5;  1e-4 8e-4;  0.05 0.20];   % CI primero: S0, I0(fijo), beta, gamma
[T,~]  = gsua_dataprep(modelo, ranges, 'domain',[0 80], ...
    'names',{'S0','I0','beta','gamma'}, 'out_names',{'S','I','R'});
T.Properties.CustomProperties.output = 2;            % serie observada (Infectados)
xdata  = linspace(0,80,200);

% TUS datos como FILA (salidas x tiempo). Aquí, conteos sintéticos de demostración:
rng(3); mu = squeeze(gsua_eval([950;5e-4;0.12],T,xdata,[],false,false,false));
ydata = poissrnd(max(mu,0));

% 1) Estimación multicomienzo (default fmincon; margin crudo > 0)
[T,res] = gsua_pe(T, xdata, double(ydata), 'N',10, 'solver','fmincon', 'margin',0.1);

% 2) Umbral de aceptación calibrado por ruido  (NO uses costf < k*min)
nf = gsua_noisefloor(T, xdata, double(ydata), T.Estfmincon, res, 'noise','quasipoisson', 'Show',true);
% nf.threshold, nf.accepted (= res < nf.threshold). Agrega 'cumulative',[..true..] si la serie es acumulada.

% 3) Ensamble paramétrico  (NO pases el pool del multistart a gsua_ia)
Est = gsua_bootstrap(T, xdata, double(ydata), T.Estfmincon, res, 'B',50, 'solver','fmincon', ...
    'noise','quasipoisson');   % 'gaussian' si los datos son continuos; agrega 'cumulative' si es acumulada

% 4) Identificabilidad + intervalos  (lee la BANDA, no el IC de la mediana)
[T,ci] = gsua_ia(T, Est, false, false, true, false);
S = gsua_clustersummary(T, ci);     % PctWidthRel <0.25 identificable / >0.5 no / dudoso
disp(S(:, {'Cluster','Parameter','PctWidthRel','Verdict'}));
```

## Reglas (haz / no hagas) — con fuente (`references/lineamientos.md`)

- **Haz** CI primero: fila i de rangos = `pars(i)` = `names{i}` (regla dura = consistencia posicional). §2
- **Haz** el modelo continuo con `ode45` + `deval` en malla fina; 3 entradas; `sol.y` salidas×tiempo. §1
- **Haz** `margin > 0` crudo (`0.1` = 10 %) y pasa `ydata` como FILA (salidas×tiempo); `gsua_pe` ya la orienta sola, pero es lo claro. §3, §6, errores §10b
- **Haz** el ruido coherente con los datos: `'gaussian'` continuo, `'poisson'`/`'quasipoisson'` conteos, `'cumulative'` acumuladas (en `gsua_bootstrap` y `gsua_noisefloor`). §5
- **No** selecciones con `costf < k·min`: usa `gsua_noisefloor` para el umbral. §7
- **No** pases el pool del multistart a `gsua_ia`: usa el ensamble de `gsua_bootstrap` (su dispersión es la incertidumbre real). §7
- **No** concluyas identificabilidad desde `T.Range` / IC de la mediana: usa la banda `PctWidthRel` de `gsua_clustersummary`. §5
- **No** leas el índice como umbral (es un ranking); `|corr|` alto = combinación (p. ej. `R0 = β·S0/γ`), no identificabilidad individual. §5
- **No** tomes un STi > 1 (o Si < 0) como bug: es ruido del estimador con N pequeño — sube N hasta que `sum(Si)/sum(|Si|)` ≈ 1. Para factores acoplados usa `'SensMethod','Shapley'` (`'Given'` sobre un ensamble). §4, `metodos_sa.md`

Lo demás está en `references/`: `lineamientos.md` (decisiones + citas), `flujos.md`, `errores.md`, `api.md`/`usage.md`; recetas completas en `recipes/`, verificación en `checks/run_checks.m`.
