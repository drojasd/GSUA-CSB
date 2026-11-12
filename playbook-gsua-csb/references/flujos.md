# Flujos de trabajo GSUA-CSB

Flujos end-to-end. Snippets mínimos; versiones ejecutables en `../recipes/`. Convenciones y citas en
`lineamientos.md`. `gsua_eval` es **posicional**:
`gsua_eval(par, T, xdata, ydata, parallel, show, verbose)` (los cuatro últimos opcionales).

## 0. Declarar un modelo

```matlab
% CI-primero: fila i de rangos = pars(i) = names{i}
ranges = [900 999;  1 50;  1e-4 8e-4;  0.05 0.20];   % S0, I0, beta, gamma
[T, solver] = gsua_dataprep('sir_continuo', ranges, 'domain', [0 80], ...
    'names', {'S0','I0','beta','gamma'}, 'out_names', {'S','I','R'}, 'opt', struct('dt',0.1));
```
Factor fijo: rango `[v v]` (sale de `T`, pero el modelo recibe el vector completo con `v` reinsertado).
Una sola salida: `T.Properties.CustomProperties.output = k`.

Python:
```python
import gsua_csb as gc
model = gc.UserFunctionModel(sir_continuo, ranges=ranges,
                             names=['S0','I0','beta','gamma'], out_names=['S','I','R'])
```

## 1. Análisis de incertidumbre (UA)

```matlab
M  = gsua_dmatrix(T, 500);                 % 500 × Np (filas = muestras)
[T, Y] = gsua_ua(T, M', 'xdata', xdata);   % gsua_ua / gsua_eval quieren Np × N
gsua_plot('UncertaintyAnalysis', T, xdata, Y);
```
Multi-salida (Kind 6, v1.18): `gsua_ua` lleva dimensión nativa de series; `xdata` por defecto = largo
de la salida nominal.

## 2. Análisis de sensibilidad (SA)

```matlab
M = gsua_dmatrix(T, N);
[T, Y, J] = gsua_sa(M, T, 'parallel', false);   % firma real: [T, Y, J]
```
- **Confiabilidad:** sube N hasta que `sum(Si)/sum(|Si|) → 1` (ver `lineamientos.md §4`). `sum(Si) → 1`
  es aditividad del modelo, cosa distinta. No hay umbral publicado.
- Multi-serie (v1.18): `Si_vec`/`STi_vec` de forma `(Np × Nd × nSeries)`; la escalarización
  multi-objetivo antigua ahora es opt-in con `'multiObjective', true`.
- Python: `res = gc.sensitivity_analysis(model, M)` (con `output_index=None` devuelve una por serie).

## 3. UA + filtrado Monte Carlo / métricas de cobertura

`gsua_covmetric` da `cost_data` (exactitud: mediana vs datos) y `cost_band` (precisión: ancho P5–P95).
`margin` aquí es banda dual-mode (`<1` crudo, `>=1` ya `1+fracción`). Con una previa sin ajustar es
normal que `cost_band >> 1` y que `cost_data` sea grande (los datos aún no son alcanzables); el módulo
05 lo reduce estimando.

## 4. Estimación de parámetros (PE)

```matlab
[T, res] = gsua_pe(T, xdata, ydata, 'N', 10, 'solver', 'fmincon', 'margin', 0.1);
[bestCost, ix] = min(res);
best = T.Estfmincon(:, ix);     % columnas = estimaciones ordenadas por costo
```
- Default: `fmincon`, `margin > 0` (crudo, `0.1 = 10 %`). `TypicalX` automático para magnitudes dispares.
- `margin` solo se registra con solvers de función de costo (`lsqc`/`lsqn` lo ignoran).
- Traza: `T.IPfmincon` (punto inicial por estimación) y `Retriesfmincon` (reinicios).
- Python: `gc.parameter_estimation(model, xdata, ydata, solver='least_squares', margin=0.1)`.

## 5. Identificabilidad práctica (IA)

```matlab
thr = gsua_noisefloor(T, xdata, ydata, best, 'Show', true);   % umbral (reemplaza costf<1.5·costf(1))
Est = gsua_bootstrap(T, xdata, ydata, T.Estfmincon, res, 'B', 50, 'solver', 'fmincon', ...
                     'noise', 'gaussian');                     % ensamble paramétrico
[T, ci] = gsua_ia(T, Est, false, false, true, false);         % cluster=false si una sola cuenca
S = gsua_clustersummary(T, ci);                                % veredicto por parámetro
disp(S(:, {'Cluster','Parameter','PctWidthRel','Verdict'}));
```
- **Lee la identificabilidad por la banda `PctWidthRel`** (<0.25 identificable / >0.5 no / dudoso),
  realización por ensamble del CSB. El **índice es un ranking**, no un umbral.
- `|corr|` alto entre factores = **combinación** identificable (p. ej. `R0 = beta*S0/gamma`), no
  identificabilidad individual.
- **Ruido del bootstrap coherente con los datos:** `'gaussian'` para series continuas; `'poisson'`/
  `'quasipoisson'`/`'nbglobal'` para conteos (con `'cumulative'` por fila para series acumuladas).
  Ruido por fila: cellstr p. ej. `{'quasipoisson','none'}` para un ajuste conjunto conteo+fracción.
- `gsua_bootstrap` corre en silencio por defecto (Display off) desde v1.17.1.

## 6. Modelo multi-salida de usuario (Kind 6, v1.18)

Una función de un argumento `y = f(pars)` que devuelve `[nSeries × nSteps]` fluye por
`gsua_ua`/`gsua_sa`/`gsua_pe`/`gsua_dia` sin adaptador escalar. La forma se descubre de la corrida
nominal o se fija con `'output_size', [nSeries nSteps]`. `gsua_fullpars(T, libres)` reinserta los
factores fijos para reconstruir el vector completo. Ejemplo: `Examples/twoScaleUserModel.m`.

## Diferencias MATLAB ↔ Python (resumen)
| | MATLAB | Python |
|---|---|---|
| `domain` al modelo | `[t0, tf]` | la rejilla de puntos pedidos |
| interpola | `gsua_eval` (lineal sobre `sol.x`) | nadie: el modelo devuelve los puntos |
| varias simulaciones | `par` Np × N (columnas) | `M` N × Np (filas) |
| forma de salida | N × tiempo × salidas | N × salidas × tiempo |
| salida transpuesta | se corrige sola (salvo caso cuadrado) | no se corrige |
