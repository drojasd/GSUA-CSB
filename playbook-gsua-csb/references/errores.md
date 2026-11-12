# Errores típicos en GSUA-CSB y cómo diagnosticarlos

Ordenados por frecuencia. Muchos vienen del módulo 00 del curso (`errores/`), donde hay modelos mal
escritos a propósito.

## Al escribir el modelo (Kind 5)

1. **Firma de dos entradas.** `gsua_dataprep` acepta `f(pars, domain)`, pero `gsua_eval` llama al
   modelo con **tres** argumentos y falla con *"Too many input arguments"*. Declara el tercero aunque
   no lo uses: `function sol = f(pars, domain, ~)`. (Modelo de muestra: `logistica_dos_entradas.m`.)

2. **Salida transpuesta.** La convención es `sol.y` = **salidas × tiempo**. Si devuelves
   tiempo × salidas, MATLAB lo corrige (reconoce el eje del tiempo por `numel(sol.x)`) **salvo** cuando
   hay tantas salidas como tiempos (matriz cuadrada). Python **no** corrige. No dependas de la
   corrección: orienta bien la salida. (Muestra: `sir_transpuesto.m`.)

3. **Devolver los pasos crudos del integrador.** `sol = ode45(...)` "funciona" (tiene `x`/`y`), pero
   `sol.x` son pasos adaptativos pocos y desiguales y `gsua_eval` interpola **linealmente** entre ellos
   → error numérico. Patrón correcto: `sol.x = t0:dt:tf; sol.y = deval(ode, sol.x);`. En Python,
   `solve_ivp` sin `t_eval` tiene el mismo problema. (Muestra: `sir_continuo_crudo.m`.)

4. **Rejilla de evaluación y `t0`.** La condición inicial vale en `t0 = domain(1)`, no en el primer
   punto que pidas. En MATLAB `t0` queda fijo en `T`; `xdata` fuera de `[t0, tf]` sale `NaN`. En Python
   el riesgo es el contrario: integrar desde `domain[0]` hace "arrancar" la dinámica en el primer punto
   pedido (muestra Python `sir_arranca_en_domain0.py`).

5. **Nombres de salida que no coinciden con las filas.** `out_names` solo **nombra** las filas de
   `sol.y`, en orden; no reordena nada. Declarar `{'I'}` para un modelo cuya fila 1 es `S` analiza `S`
   con la etiqueta `I`, **sin error**.

6. **Unidades.** Rangos, condiciones iniciales y dominio deben estar en las mismas unidades. En un SIR,
   `beta` en 1/(persona·día) porque el término es `beta*S*I` con `S`,`I` en personas; escribir `S0`,`I0`
   como fracciones sin reescalar `beta` hunde `R0 = beta*S0/gamma`.

## De convención / orden

7. **Orden de factores inconsistente.** El único requisito duro: fila *i* de rangos = `pars(i)` =
   `names{i}` = lo que el modelo desempaqueta en la posición *i*. Si reordenas el modelo, reordena
   también `ranges`, `names`, los vectores de parámetros literales y **las lecturas posicionales del
   nominal** (`T.Nominal(k)`) y de las estimaciones. Un desalineo **no da error**: da números malos
   (p. ej. un `R0` absurdo). Convención del curso: CI-primero (ver `lineamientos.md §2`).

8. **Factor fijo "corre" el desempaque — no.** Al fijar un factor con rango `[v v]` sale de `T`, pero
   el modelo **recibe el vector completo** con `v` reinsertado en su posición. No añadas defaults por
   `length(pars)` pensando que el fijo desaparece del vector. Usa `gsua_fullpars` si necesitas
   reconstruir el vector ordenado a partir de los libres.

9. **`margin` con la convención vieja.** Desde v1.18 `gsua_likelihood` recibe la fracción **cruda**
   (`0.08` = 8 %), no el valor desplazado (`1.08`). Un `margin >= 1` dispara warning `MarginConvention`.
   En `gsua_covmetric`/`gsua_rcostf`/`gsua_costfMulti` la banda sigue siendo dual-mode. Ver
   `lineamientos.md §3`.

10. **`gsua_eval` con nombre–valor.** `gsua_eval` es **posicional**
    (`gsua_eval(par, T, xdata, ydata, parallel, show, verbose)`). Pasarle `'xdata', ...` como par
    nombre–valor no funciona.

10b. **`gsua_pe` con el solver por defecto (`fmincon`) y `ydata` en columna.** `gsua_pe` espera
    `ydata` como **salidas × tiempo** (fila, para una sola salida), igual que devuelve `gsua_deval`.
    `lsqc`/`lsqn` se autocorrigen (capturan el *size mismatch* y transponen `ydata`), pero
    **`fmincon`** (el default ratificado) y los metaheurísticos **no**: un `ydata` en **columna**
    dispara `MATLAB:badsubscript` (*"logical indices contain a true value outside of the array
    bounds"*), se reintenta 100 veces desde puntos aleatorios y termina en *"Every estimation
    failed"* — un mensaje que no apunta a la causa real (la orientación). **Pasa `ydata` como fila**
    (salidas × tiempo). Ojo: el patrón del curso `squeeze(gsua_eval(...)).'` produce una **columna**;
    funciona con `lsqc` (como en `demo_05`), pero no con el default `fmincon`. **Pendiente de Daniel**
    (mejora de robustez propuesta): orientar `ydata` de forma genérica al entrar a `gsua_pe`, como ya
    hace `gsua_noisefloor` desde v1.17.1.

## De análisis

11. **Ruido del bootstrap incoherente con los datos.** `'poisson'`/`'quasipoisson'` asumen conteos
    enteros; aplicarlos a una serie continua (o a una fila de fracciones 0–1) introduce ruido y ceros
    espurios. Usa `'gaussian'` para continuo; para un ajuste conjunto conteo+fracción, ruido por fila
    (`{'quasipoisson','none'}`). Serie acumulada: `'cumulative'`.

12. **Concluir identificabilidad desde el IC de la mediana** (`gsua_ia.Range`). Usa la banda
    `PctWidthRel` de `gsua_clustersummary` (el **índice es un ranking**, no un umbral). `|corr|` alto =
    combinación, no no-identificabilidad individual. Ver `lineamientos.md §5`.

13. **Fijar N con una regla inventada.** No existe `sum(Si) > 0.65` ni umbral publicado. Sube N hasta
    que `sum(Si)/sum(|Si|)` se estabilice cerca de 1. Ver `lineamientos.md §4`.

## Históricos (resueltos en v1.18)

14. **Kind 6 (modelo matriz) reventaba `gsua_pardeval`** dentro de `parfor` con *"indices on the left
    side are not compatible"*. Resuelto: la forma se descubre del nominal (o `'output_size'`) y
    `gsua_ua`/`gsua_sa` llevan dimensión de series.

15. **Costo con regulador cero.** Una fila de datos todo-cero/constante daba regulador cero y hacía
    fallar a `fmincon` "objective undefined at initial point" sin nombrar la fila. Resuelto: se excluye
    la fila con un warning `ZeroRegulator`.

## Entorno

16. **MATLAB MCP "failed to attach".** Relanzar `matlab.exe` con `matlab.engine.shareEngine`
    (ver memoria `reference_matlab_mcp_reattach`).
17. **`twine` rompe `streamlit`** al instalarse en el entorno base (conflicto `packaging`). Instalar
    `twine` por `pipx` o en un venv aislado, no en base (ver `python/PUBLISHING.md`).
