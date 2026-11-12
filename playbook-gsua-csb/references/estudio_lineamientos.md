# Estudio de lineamientos GSUA-CSB (cómo escribir modelos y cómo usar la toolbox)

Convención de citas:
- **[TB]** = repo de la toolbox, `D:\MATLAB Drive.source-of-truth-20260718T095207\Toolbox` (HEAD `81276f0`, "v1.18.0 / 0.4.0", 2026-11-05).
- **[SM]** = repo del semillero, `D:\Other-projects\Semillero Matemáticas` (HEAD `0174225`). `EPI/` abrevia `cursos/Modelado epidemiológico/`.
- `archivo:N-M` = líneas; `§` = sección/encabezado. Las citas a `Functions/*.m` son al código actual (fuente de verdad cuando la documentación discrepa).
- No ejecuté MATLAB. Lo marcado "(no ejecutado)" es inferencia por lectura de código.

---

## 1. Cómo escribir un modelo para GSUA

### 1.1 Tipos de modelo y cómo entran a la toolbox
| Tipo | Llamada | Fuente |
|---|---|---|
| Simulink | `[T,tout] = gsua_dataprep(model_name, Ranges, Parameter_names)` | [TB] `docs/gsua_userguide.m:10,239` |
| EDO simbólicas | `[T,solver] = gsua_dataprep(odes, vars, domain, modelName)` | [TB] `docs/gsua_userguide.m:11,254` |
| Función de usuario | `[T,sol] = gsua_dataprep(func, Ranges, 'domain', domain)` (`func` = nombre del archivo, **char**, no handle) | [TB] `docs/gsua_userguide.m:12,273`; `course/00_Escribir_Modelos/matlab/demo_00_escribir_modelos.m:92` |
- `gsua_dataprep` despacha por `exist(p1,'file')`: 4 → Simulink, 2 → función de usuario; si no es texto → simbólico ([TB] `Functions/gsua_dataprep.m:39-65`). La función debe ser un archivo `.m` propio **en el path** (`demo_00...m:46-47`, `:372`).
- El 2.º output para función de usuario es el *handle* del solver, no una solución (`Functions/gsua_userdefined.m:1,44`).

### 1.2 Firma obligatoria (función de usuario con dominio, Kind 5)
- **MATLAB:** `sol = f(pars, domain, opt)` — **tres entradas siempre**, aunque `opt` no se use (`~`) ([TB] `course/00_Escribir_Modelos/README.md:9-11`; `demo_00...m:45-67,220-224,372`; `logistica_discreta.m:14-15`).
  - Motivo en código: con 3 entradas declaradas se construye `solver=@(pars,domain,opt)` (`Functions/gsua_userdefined.m:166-167,180-181`), y la evaluación Kind 5 llama siempre `fun(par,domain,opt)` (`Functions/gsua_deval.m:57-60`; `Functions/gsua_eval.m:144-160`). Con 2 entradas el handle es `@(pars,domain)` y la llamada con 3 argumentos falla.
  - El número de entradas declaradas decide el "Kind": 1 entrada → Kind 6 (o 2 si `'vectorized',true`), 2 → handle de 2 args, 3 → Kind 5 (`gsua_userdefined.m:15-18,151-183`).
- **Python:** `func(pars, domain, **opt)` que devuelve un arreglo `(n_salidas, len(domain))` evaluado en los puntos de `domain` ([TB] `course/00_Escribir_Modelos/README.md:10-11`).
- Modelo "libre de tiempo": `out = f(factors)` que devuelve escalar o vector; si devuelve matriz `nSeries x nSteps` (v1.18, Kind 6) funciona en `gsua_ua/sa/pe/dia` sin adaptador; forma detectada con el nominal o `'output_size'` ([TB] `docs/gsua_userguide.m:272,281`; `CHANGELOG.md:19-23`).

### 1.3 Entradas
| Argumento | Regla | Fuente |
|---|---|---|
| `pars` | La **posición manda**: `pars(i)` = fila *i* de `Ranges` = `names{i}`. Convención del curso: **primero parámetros, después condiciones iniciales**. | [TB] `course/00.../README.md:12-13`; `demo_00...m:53-56,93-95,373` |
| `pars` (forma) | Llega como vector **columna**: el handle hace `user_func(pars',...)` sobre una fila (`gsua_userdefined.m:181`; `gsua_eval.m:59,159`). Los modelos del hands-on hacen `pars(:).'` y aceptan fila o columna ([SM] `EPI/hands_on_epi/matlab/modelos/sir_continuo.m:29`). | código |
| `pars` con factores fijos | El modelo **siempre recibe el vector completo** (el fijo se reinserta en su posición) aunque `T` y `par` solo lleven los libres. | `demo_00...m:200-204`; `gsua_userdefined.m:153-168,271-274` |
| `domain` | `[t0, tf]`, el que se declara con `'domain'` en `gsua_dataprep`; la toolbox no lo cambia → la condición inicial vale en `t0 = domain(1)`. Default de `'domain'` = `1`. | `demo_00...m:57-59`; `gsua_userdefined.m:25-26,53` |
| `opt` | Lo que se pase con `'opt'` llega tal cual como 3.er argumento. Default `struct('action',1)`. | `demo_00...m:60-61,161-163`; `gsua_userdefined.m:62`; `docs/gsua_userguide.m:283` |
| `opt` (hands-on) | Solo **elecciones estructurales**, nunca parámetros a estimar (p. ej. `opt.prob`, `opt.N`). | [SM] `EPI/hands_on_epi/CONTRATO.md:46-48` |

### 1.4 Salida
- Estructura con `sol.x` (1 × Nt, de `t0` a `tf`) y `sol.y` (**nSalidas × Nt**: una fila por salida, en el orden de `out_names`) ([TB] `course/00.../README.md:9-11`; `demo_00...m:63-67,374`).
- `gsua_eval` **interpola linealmente** `sol.y` en `xdata`; por eso `sol.x` debe cubrir todo el dominio y ser lo bastante fino (`demo_00...m:69-71`; `Functions/gsua_intrp.m:113-119`, `extra=2` usado por Kind 5 en `gsua_deval.m:60`).
- Puntos de `xdata` fuera de `[t0,tf]` salen `NaN` (interp1 sin extrapolación) (`demo_00...m:261-264`).
- Orientación: si `sol.y` viene `tiempo × salidas`, desde v1.13 la toolbox la transpone usando `numel(sol.x)`; si es cuadrada (nSalidas = Nt) **no** puede y la deja; también acepta `[t,y]=ode45` (`gsua_intrp.m:35-68`; `demo_00...m:246-250`). Regla: **no depender de la corrección** (`demo_00...m:250`).
- Una forma incompatible (ninguna dimensión = `numel(sol.x)`) lanza `GSUA:gsua_intrp:BadShape` (`gsua_intrp.m:60-64`).
- Hands-on: las filas de `sol.y` son primero los estados y al final las **salidas observables**; `sol.solver = ''` ([SM] `EPI/hands_on_epi/CONTRATO.md:31-34,40-42`).

### 1.5 Manejo del tiempo
| Caso | Patrón | Fuente |
|---|---|---|
| Discreto | Bucle de `t0` a `tf`, `t = domain(1):domain(2)` (un punto por paso). | [TB] `course/00.../matlab/logistica_discreta.m:24-29`; `demo_00...m:375` |
| Continuo | `ode45` **dentro** de la función sobre `domain`, luego `deval` en una rejilla fina y regular; devolver esa rejilla. Curso: `RelTol 1e-6, AbsTol 1e-8`, `dt` por `opt.dt` (0.1 por omisión). | [TB] `course/00.../matlab/sir_continuo.m:27-42`; `demo_00...m:146-150,375` |
| Continuo (hands-on) | `ode45` con `RelTol = AbsTol = 1e-6` y `deval` en `t0:tf` paso 1 (misma rejilla que los discretos); Python `solve_ivp(method="RK45", t_eval=..., rtol=1e-6, atol=1e-6)`. | [SM] `EPI/hands_on_epi/CONTRATO.md:31-39` |
| Incidencia (hands-on) | Convención hacia atrás `casos(t) = C(t) − C(t−1)`; `casos(t0)` = infectados iniciales; en EDO se agrega el estado acumulado `C' = incidencia`. | [SM] `CONTRATO.md:43-45`; `matlab/modelos/sir_continuo.m:51` |
| Determinismo (hands-on) | Sin estado global ni aleatoriedad; paridad MATLAB/Python 1e-8 (discretos), 1e-4 (continuos). | [SM] `CONTRATO.md:49-51` |

### 1.6 Nombres y salidas
- `'names'`: nombres de factores (cellstr, string o char; default `1..Np`; TeX/LaTeX se renderiza) ([TB] `docs/gsua_userguide.m:279`; `gsua_userdefined.m:23-24,67,140-144`).
- `'out_names'`: nombran las filas de `sol.y` **en orden; no reordenan nada** (default `{'out'}`) (`demo_00...m:95,236-239`; `gsua_userguide.m:280`).
- `'output'`: índices de salidas activas (default todas); también `T.Properties.CustomProperties.output = k`, que afecta a todas las funciones que reciben `T` (`demo_00...m:186-191,378`; `gsua_userguide.m:32,278`).
- Factor fijo: rango degenerado `[v v]`; **desaparece de `T`** y del muestreo/estimación (`demo_00...m:200-204,378`; `gsua_userguide.m:13`).
- Kind 6 multi-serie: si `out_names` no coincide con nSeries, warning `OutNamesMismatch` y auto-numeración `out1..` (`gsua_userdefined.m:202-210`).

### 1.7 Qué revisa `verificar_modelo` (ayuda del curso, no de la toolbox)
MATLAB, `verificar_modelo(nombre, pars, domain, out_names[, opt])` ([TB] `course/00_Escribir_Modelos/matlab/verificar_modelo.m:7-16`):
1. `nombre.m` está en el path (`exist == 2`) y declara **3 entradas** (`:22-30`).
2. Devuelve estructura con campos `x` e `y` (`:32-35`).
3. `sol.y` es `numel(out_names) × numel(sol.x)`; si viene transpuesta da **aviso** (no falla) (`:37-48`).
4. `sol.x` va de `domain(1)` a `domain(end)` (tolerancia 1e-9) (`:50-56`).
5. Sin `NaN` ni `Inf` (`:58-62`).
- **No** detecta falta de precisión: el modelo crudo de ode45 pasa (`demo_00...m:326-327`).
- Python (`python/verificar_modelo.py:7-10,46-54`) además exige que el valor en un tiempo no dependa de la rejilla pedida (recorta a la segunda mitad y compara, `rtol 1e-4, atol 1e-6`) y trata la transpuesta como **FALLA** (`:33-35`).

### 1.8 Errores deliberados del módulo 00
| # | Archivo | Qué está mal | Cómo se detecta | Fuente |
|---|---|---|---|---|
| a | `matlab/errores/logistica_dos_entradas.m` | Firma con 2 entradas `(pars, domain)`; `gsua_dataprep` la acepta pero `gsua_eval` falla ("Too many input arguments"). | `verificar_modelo` FALLA en la firma; `gsua_eval` lanza error. | `logistica_dos_entradas.m:1-6`; `demo_00...m:220-233` |
| a-Py | `python/errores/logistica_sin_opt.py` | Sin argumento para `opt`: con `opt={"t0":0}` Python falla "unexpected keyword argument". | `verificar_modelo.py` FALLA al llamar. | `logistica_sin_opt.py:4-6` |
| b | (en el demo, sin archivo) | `out_names={'I'}` para un modelo cuya fila 1 es S: se analiza S con etiqueta I **sin error**. | Solo inspección (valor inicial = S0). | `demo_00...m:236-243` |
| c | `matlab/errores/sir_transpuesto.m` | Devuelve `sol.y` Nt × 3. MATLAB la corrige sola (salvo caso cuadrado); Python no (`set_output("I")` toma un instante). | `verificar_modelo`: aviso (MATLAB) / FALLA (Python). | `sir_transpuesto.m:4-9`; `python/errores/sir_transpuesto.py:4-11`; `demo_00...m:246-258` |
| d-Py | `python/errores/sir_arranca_en_domain0.py` | Integra desde `domain[0]`: con `xdata` desde el día 20 la epidemia "arranca" el día 20. Sin equivalente MATLAB (allí `t0` queda fijo en `T`; fuera de `[t0,tf]` sale `NaN`). | `verificar_modelo.py` (prueba de rejilla recortada). | `sir_arranca_en_domain0.py:4-6,16-17`; `demo_00...m:261-272`; `README.md:59-61` |
| e | `matlab/errores/sir_continuo_crudo.m` (+ `.py`) | Devuelve `ode45` tal cual (tolerancias por defecto): `sol.x` = pasos adaptativos pocos y desiguales; `gsua_eval` interpola linealmente entre ellos → error de precisión. | **No** lo detecta `verificar_modelo` (cumple el contrato); se ve comparando contra la versión con `deval`. | `sir_continuo_crudo.m:1-15`; `demo_00...m:284-299,326-327` |
| f | (en el demo) | Unidades inconsistentes: S0/I0 como fracción sin reescalar `beta` hunde `R0 = beta*S0/gamma`; dominio en semanas con `gamma` en 1/día multiplica el tiempo de recuperación por 7. | Solo razonamiento/inspección. | `demo_00...m:312-321` |

Plantilla y solución del ejercicio: `sir_discreto_plantilla.m` ya cumple el contrato (3 entradas, `sol.x = domain(1):domain(2)`, `sol.y` 3 × Nt) y solo le falta la dinámica (TODO deja S, I, R constantes) (`sir_discreto_plantilla.m:22-37`); la solución actualiza `[S - βSI; I + βSI − γI; R + γI]` (`solucion/sir_discreto_solucion.m:13-18`). El demo comprueba la plantilla mirando si I(t) sigue constante (`demo_00...m:360-369`). Preguntas del ejercicio: por qué no coincide con el continuo, qué pasa si `beta*S0 > 1`, si se conserva S+I+R (`demo_00...m:349-351`).

---

## 2. Flujo de trabajo

### 2.1 Secuencia canónica
Curso: Escribir modelo → Modelo → Muestreo → Incertidumbre → Sensibilidad → Estimación + Identificabilidad ([TB] `course/README.md:18-34`). README: `gsua_dataprep → gsua_dmatrix → gsua_sa/gsua_ua/gsua_pe → gsua_ia` ([TB] `README.md:139-157`). Guía: dataprep, simular, muestrear, UA, SA, estimar, identificabilidad, intervalos ([TB] `docs/gsua_userguide.m:4`). Ciclo semiautomático: alcanzabilidad → multistart → identificabilidad → delgadez del IC (`coverage_metric`) → acción correctiva → repetir ([TB] `python/USERGUIDE.md:447-461`).

| Pregunta | Función | Entradas y formas | Fuente |
|---|---|---|---|
| ¿Simula `T`? | `[y,xdata] = gsua_eval(values, T, xdata, ydata, parallel, show, verbose)` | `values` **Np × N** (columnas = simulaciones); argumentos finales **posicionales**; salida **N × numel(xdata) × nSalidasActivas** | `gsua_userguide.m:19-25`; `demo_00...m:122-126,377`; `Functions/gsua_eval.m:30-39` |
| Diseño de experimentos | `[M,T2] = gsua_dmatrix(T, N, ...)` | `M` **N × Np** (filas = muestras) → para `gsua_eval` pasar `M'` | `gsua_userguide.m:35-36`; `demo_00...m:182-183` |
| Incertidumbre Monte Carlo (+ MCF automático) | `Y = gsua_ua(M, T, ...)` | `M` N × Np; `Y` N × len (× nout) | `gsua_userguide.m:56-57`; `course/03.../demo_03_incertidumbre.m:79-80` |
| Cobertura y delgadez de la banda | `[cost_data,cost_band,P5,P50,P95] = gsua_covmetric(Y, ydata, ...)` | `Y` reps × len (× nout); `ydata` nout × len | `gsua_userguide.m:176`; `Functions/gsua_covmetric.m:24-28` |
| MCF por separado | `[Jsup,Jinf] = gsua_MCF(T, M, Y, y_exp, tIn, tOut)` | | `gsua_userguide.m:228` |
| Sensibilidad | `[T,Y,J] = gsua_sa(M, T, ...)` | índices en `T.Si`, `T.STi`, `T.Si_vec`, `T.STi_vec` (`Np × Nd [× nSeries]`) | `gsua_userguide.m:69,78`; `course/04.../demo_04_sensibilidad.m:60-62` |
| Estimación | `[T,res] = gsua_pe(T, xdata, ydata, ...)` (3.er output `IP`) | `ydata` salidas × tiempo; resultados `T.Est<solver>` **Np × N** ordenados de mejor a peor, `T.IP<solver>`, propiedad `Retries<solver>` | `gsua_userguide.m:91-92`; `Functions/gsua_pe.m:1-16` |
| Filtrar corridas malas | `keep = gsua_costcutoff(cost, method, rtol, atol, gapRatio, minKeep)` (6 posicionales obligatorios) | | `gsua_userguide.m:189` |
| Umbral de aceptación calibrado por ruido | `out = gsua_noisefloor(T, xdata, ydata, Est, res, ...)` | `Est` Np × N; `out.threshold`, `out.accepted`, `out.byModel` | `gsua_userguide.m:198,206` |
| Bootstrap paramétrico (incertidumbre honesta) | `[Est,costs,out] = gsua_bootstrap(T, xdata, ydata, Est0, res, ...)` | `Est` Np × B → a `gsua_ia` | `gsua_userguide.m:136`; `Functions/gsua_bootstrap.m:62-66` |
| Identificabilidad práctica | `[T,clusterInfo] = gsua_ia(T, T_est, correction, outlier, show, cluster, maxK, silThreshold, ...)`; `gsua_dia` = sin figuras | `T_est` **Np × N** (columna = corrida) | `gsua_userguide.m:109-110,224`; `python/USERGUIDE.md:477` |
| Resumen por cluster/parámetro | `S = gsua_clustersummary(T, clusterInfo)` | veredicto <25 % identificable, >50 % no, en medio dudoso | `gsua_userguide.m:124,131`; `CHANGELOG.md:99-102` |
| Muestreo conjunto tras IA | `gsua_dmatrix(Tia, N, 'Method','Joint', ...)` | pool `Tia.Est` (Np × nPool) o `'Pool'` | `gsua_userguide.m:41-52` |
| IC por perfil de verosimilitud | `Tfinal = gsua_likelihood(T, xdata, ydata, alpha, step, margin, tol1, tol2, limit, reps, show, parallel, saver, pars[, verbose])` (todo posicional) | | `gsua_userguide.m:216-219` |
| Rango OAT / CSB | `T = gsua_oatr(T, ...)`; `[T,New_range,J_test,Y_test,sup] = gsua_csb(T, N, ...)` | | `gsua_userguide.m:148-163` |
| Reinsertar fijos | `[full,names] = gsua_fullpars(T, freeValues)` | | `gsua_userguide.m:230` |
| Persistencia | `gsua_save` / `gsua_load` (reconstruye el handle) | | [TB] `README.md:204-208` |

Flujo de identificabilidad recomendado en la versión actual (1.17+): `gsua_pe` (multistart, con `'margin'` y solver no-lsq para el costo `gsua_costf`) → `gsua_bootstrap` → `gsua_ia(T, Eb, ...)` → `gsua_clustersummary` (`gsua_userguide.m:142-146`; `Functions/gsua_bootstrap.m:71-74`; `course/05.../demo_05_estimacion.m:97-104,137`). Alternativa: `gsua_pe` → `gsua_noisefloor` → `gsua_ia(T, Est(:,out.accepted))` (`gsua_userguide.m:208-212`).

### 2.2 Opciones nombre-valor y defaults (guía, verificados contra código/`functionSignatures.json`)
**`gsua_dataprep` / `gsua_userdefined` (función de usuario)** — `'domain'` 1; `'rMethod'` `'range'` (`'range'|'std'|'percent'|'normal'`); `'nominal'` punto medio (con `'normal'`, la 1.ª columna); `'output'` todas; `'names'` `1..Np`; `'out_names'` `{'out'}`; `'output_size'` `[]`; `'vectorized'` `false` (solo 1 entrada); `'opt'` `struct('action',1)` (`gsua_userguide.m:274-283`; `Functions/gsua_userdefined.m:52-75,107-139`; `functionSignatures.json` `gsua_userdefined`). Simbólico: `'output'` `1:model_order`, `'solver'` `'ode45'|'ode4'` (ode4: 3.ª columna de `domain` = paso; solo `NonNegative` en `opt`) (`gsua_userguide.m:256-260`). Simulink: `'modelkind'` `'Dynamic'|'Static'`, `'Step'`, una sola salida (`gsua_userguide.m:239-244`).

**`gsua_dmatrix`** — `'Method'` `'LatinHypercube'` (`'Uniform'|'Sobol'|'Joint'`); `'Show'` `'off'`; `'Seed'` `[]`; solo Joint: `'JointType'` `'Bootstrap'` (`'SmoothBootstrap'|'Gaussian'`), `'Pool'` `[]`, `'MinPoolN'` 5 (≥2), `'Bandwidth'` `[]` (Silverman), `'Clip'` `[]` (`gsua_userguide.m:38-46`; `Functions/gsua_dmatrix.m:79-97`).

**`gsua_ua`** — `'xdata'` `[]` (grilla del dominio), `'ynom'` `[]` (salida en `T.Nominal`; mismo largo que `xdata`), `'parallel'` **true**, `'verbose'` true, `'plot'` true (`gsua_userguide.m:59-63`; `Functions/gsua_ua.m:32-36`).

**`gsua_sa`** — `'SensMethod'` `'Xiao'` (`'Saltelli'|'Jansen'` mismo costo `N*(Np/2+1)`; `'Sobol'` y `'OAT'` `N*(Np+1)`; `'brute-force'` `N+2*Np*N^2`); `'ynom'` `[]` (no cambia los índices); `'xdata'` `[]`; `'parallel'` **true**; `'verbose'` true; `'xAlpha'` 1 en (0,2]; `'bandwidth'` 1; `'multiObjective'` false (`gsua_userguide.m:71-78`; `Functions/gsua_sa.m:38-63`).

**`gsua_pe`** — `'N'` 1; `'solver'` `'lsqc'` (`lsqn|fmincon|ga|particle|psearch|surrogate|annealing`); `'opt'` `[]` (optimoptions del optimizador correcto; si no trae `TypicalX`, en `fmincon/lsqc/lsqn` se llena con `max(abs(T.Nominal),eps)`); `'Multistart'` 1 (solo lsqc); `'ipoint'` `[]` (N filas × libres); `'margin'` 0; `'alpha'` 2; `'A','B','Aeq','Beq','nonlcon'` `[]`; `'Show'` `'off'` (`'on'`, `'on2'` llama a `gsua_ia`); `'save'` **false** (true antes de v1.15); `'timer'` true (`gsua_userguide.m:94-103`; `Functions/gsua_pe.m:52-88`).

**`gsua_ia`** (posicionales) — `correction` false, `outlier` false, `show` false, `cluster` false, `maxK` 5, `silThreshold` 0.5; nombre-valor: `'cost'` `[]`, `'CostMethod'` `'rtol'` (`'gap'`), `'CostRtol'` 0.1, `'CostAtol'` 1e-8, `'CostGapRatio'` 3, `'MinKeep'` 3, `'IndexMethod'` `'spread'` (`'ci'` = legado), `'truth'`, `'xdata'`, `'ydata'` `[]` (`gsua_userguide.m:110-122`; `Functions/gsua_ia.m:178-193`). Efecto: `T.Range` ← IC de la **mediana** (`gsua_medianCI`), `T.Nominal` ← mediana, `T.index` (0 bien, 1 mal), `T.Est` ← estimaciones usadas (`gsua_userguide.m:109`).

**`gsua_bootstrap`** — `'B'` 200; `'solver'` **`'fmincon'`**; `'N'` 20; `'margin'`/`'alpha'` recuperados de `T` (si no hay, `margin` 0 = MSE); `'noise'` `'quasipoisson'` o cellstr por fila con `'gaussian'`/`'none'`; `'cumulative'` todo false; `'seed'` `[]`; `'parallel'` false; `'opt'` Display-off por defecto; `'timer'` true (`gsua_userguide.m:137-140`; `Functions/gsua_bootstrap.m:84-94,151-160`).

**`gsua_noisefloor`** — `'B'` 2000; `'noise'` `'quasipoisson'` (`'poisson'|'nbglobal'`); `'quantile'` 0.95; `'cumulative'` todo false (uno por fila); `'seed'`; `'margin'`/`'alpha'` recuperados de `T` (si faltan → error `NoMarginAlpha`); `'Show'` false (`gsua_userguide.m:199-205`; `Functions/gsua_noisefloor.m:118-125,154-176`).

**`gsua_covmetric`** — `'margin'` 0.1, `'alpha'` 2, `'plevels'` `[5 95]` (`gsua_userguide.m:177-179`; `Functions/gsua_covmetric.m:56-58`).

**`gsua_oatr`** — `'lim'` 0.3, `'parallel'` false, `'correct'` false, `'show'` false, `'titerlimit'` inf. **`gsua_csb`** — `'reps'` 100, `'stop'` 0.95, `'lim'` 0.3, `'select'` 0.5, `'switch'` 0.75, `'ynom'` `[]`, `'correct'` true, `'parallel'` false, `'show'` false; además en código `'protect'` false, `'breaking'` 1, `'stretch'` 3 (no documentados en la guía) (`gsua_userguide.m:151-163`; `Functions/gsua_csb.m:30-54`).

**`gsua_medianCI`** `alpha` 0.05; **`gsua_eval`** `xdata []`, `ydata []`, `parallel false`, `show true`, `verbose true` (`gsua_userguide.m:21-25,226`; `gsua_eval.m:35-39`).

---

## 3. Reglas y convenciones

### Modelo y tabla
1. Archivo `.m` propio, en el path, con 3 entradas; `pars(i)` = fila *i* de rangos = `names{i}`; parámetros primero, CI después ([TB] `demo_00...m:372-373`).
2. Rangos `Np × 2` `[mín, máx]`; nominal por defecto = punto medio (`demo_00...m:93,109-110`; `gsua_userguide.m:13`).
3. Unidades coherentes entre rangos, CI y dominio (`demo_00...m:312-318`).
4. Continuos: `ode45` + `deval` en rejilla fina; no devolver los pasos del integrador (`demo_00...m:146-150,284-288`; `sir_continuo.m:38-42`).
5. Comprobar cada modelo nuevo con `verificar_modelo` antes de usarlo (`demo_00...m:323-328`).
6. Antes de entregarlo a la toolbox, llamarlo a mano y mirar las formas (`demo_00...m:84-89`).

### Muestreo y rangos
7. GSA con **rangos amplios** (todo el rango plausible); UA de "incertidumbre plausible" con rangos más estrechos (~±20 % de la verdad en el curso) ([TB] `course/04.../demo_04_sensibilidad.m:40-41`; `course/03.../demo_03_incertidumbre.m:42-44`).
8. Parámetros que abarcan varios órdenes de magnitud: muestrear en log (en MATLAB solo se ilustra a mano; `log_scale` es solo Python) (`course/02.../demo_02_muestreo.m:100-112`; `python/USERGUIDE.md:178-185,494`).
9. Muestreo **conjunto** (`'Joint'`) cuando hay factores correlacionados (heatmap de `gsua_ia` o `cost_band ≫ cost_data`); no recortar (`'Clip'`) al `T.Range` de `gsua_ia` porque es un IC de la mediana; pool < `MinPoolN` (5) se rechaza y cae a LHS (`gsua_userguide.m:42-48`).
10. Tamaños usados en el curso: 250 (muestreo), 300 (UA), 500 (SA); la guía usa 500 y 2000 (`demo_02...m:62`; `demo_03...m:71`; `demo_04...m:56`; `gsua_userguide.m:47,182`). No hay regla numérica explícita de N; sí la de convergencia: **`Si` negativos (Saltelli/distancia) indican N demasiado pequeño → subirlo** (`gsua_userguide.m:88`); y probar que los índices se estabilizan al subir N (`demo_04...m:109`).

### Sensibilidad
11. `STi` es la medida clave; `sum(STi) ≥ 1`, `sum(Si) ≤ 1`; si `sum(Si) > 0.65` interacciones despreciables y AS local/analítico es apropiado (`gsua_userguide.m:88`).
12. Estimar los de `STi` alto; considerar **fijar** los de `STi ≈ 0`; `STi ≫ Si` = interacciones (`demo_04...m:72,82`).
13. Comparar estimadores (Xiao vs Saltelli) para robustez, partiendo de una tabla fresca (reutilizar `T` con propiedades de Xiao da un warning inofensivo) (`demo_04...m:90-99`).
14. En MCF: mayor separación entre curvas bajo/alto = factor más importante; con `'ynom'` específico la lectura es más sutil y la AS es más limpia (`gsua_userguide.m:65`).

### Estimación y costo
15. `N > 1` siempre: una sola estimación no dice si el optimizador halló *la* respuesta (`python/USERGUIDE.md:263-265`).
16. `'margin'` es la banda relativa **cruda** (`0.3` = 30 %), usada internamente como `m+1`; solo reescala el costo, no mueve el óptimo; `margin=0` → MSE; `margin<0` → NLL gaussiana; con costo `gsua_costf`, `cost < 1` ≈ residuo dentro de esa fracción del RMS de los datos (`gsua_userguide.m:99`; `Functions/gsua_pe.m:35-46`).
17. El costo `gsua_costf` (`margin` ≠ 0) **solo** lo usan `ga, particle, psearch, surrogate, annealing, fmincon`; `lsqc/lsqn` lo ignoran y su `res` es la suma de cuadrados del optimizador (`Functions/gsua_pe.m:35-36,115,229-276`). Es el que conviene "cuando un modelo debe ajustar varias señales a la vez" (`gsua_userguide.m:99`).
18. Restricciones `A,B,Aeq,Beq,nonlcon` solo llegan a `ga`, `psearch`, `fmincon` (`gsua_pe.m:277-290`).
19. Silenciar el optimizador con `optimoptions(..., 'Display','off')`; `'save'` false para no sobrescribir `Estimations.mat` (`demo_05...m:59,69-70`; `gsua_userguide.m:102`; `CHANGELOG.md:79-84`).
20. Filas de datos constantes o todo-cero se excluyen del costo con warning `ZeroRegulator` (v1.18) (`CHANGELOG.md:50-53`; `gsua_pe.m:250-266`).

### Aceptación de corridas, piso de ruido, bootstrap
21. Filtrar corridas mal convergidas antes de cualquier estadística (`'cost', res` en `gsua_ia`; `rtol` 10 % o `gap` ≥ 3× la brecha mediana; `MinKeep` 3) (`gsua_userguide.m:118-120,188-189`).
22. La regla `res < 1.5*res(1)` es dependiente de escala y "al revés"; reemplazarla por `gsua_noisefloor` (`gsua_userguide.m:197`; `Functions/gsua_noisefloor.m:8-12`).
23. Modelo de ruido: **no** usar Poisson por defecto (casi siempre demasiado estrecho); `quasipoisson` recomendado; `nbglobal` reportar solo como cota; comparar los tres (`out.byModel`) (`gsua_userguide.m:200`; `gsua_noisefloor.m:74-90`).
24. Si el mejor ajuste excede el piso, se emite warning: el modelo se rechaza como descripción de los datos (`gsua_userguide.m:206`).
25. **Series acumuladas:** marcar `'cumulative'` por fila; el ruido se estima y genera sobre la serie incidente (diferenciada) y se re-acumula; ruido independiente sobre una suma acumulada "es incorrecto" (`gsua_userguide.m:202`; `gsua_noisefloor.m:57-61`).
26. **Multistart vs bootstrap:** multistart responde "¿se halló el óptimo global?" y colapsa a un óptimo; bootstrap responde "¿cuánto cambiaría la estimación con otro ruido?". Pasar a `gsua_ia` el ensamble **bootstrap**, no el pool multistart (`gsua_userguide.m:135`; `Functions/gsua_bootstrap.m:9-13,25-28`).
27. Ajuste conjunto conteo + fracción: `'noise'` por fila, p. ej. `{'quasipoisson','none'}` (`gsua_userguide.m:139`; `CHANGELOG.md:30-32`).
28. El piso es una calibración de banda, no una región de confianza formal; condicional a que θ* sea aproximadamente correcto (optimista bajo mala especificación); el 95 % es convención (`gsua_noisefloor.m:24-30`).

### Identificabilidad e intervalos
29. `T.Range` tras `gsua_ia` es IC de la **mediana**, no del parámetro; el hands-on y los retos reportan percentiles del ensamble (`gsua_userguide.m:109`; [SM] `EPI/hands_on_epi/README.md:50-53`; [SM] `retos/README.md:71-74`).
30. `IndexMethod='spread'` (default desde 1.17): ancho 2.5–97.5 % sobre el rango **original**; `'ci'` encoge con el tamaño de muestra (`gsua_userguide.m:121`; `CHANGELOG.md:121-124`).
31. Veredicto por ancho relativo al rango: <25 % identificable, >50 % no, en medio dudoso (`CHANGELOG.md:99-102`; [SM] `hands_on_epi/README.md:50-51`).
32. `|corr| > 0.9` → los datos identifican la combinación, no cada parámetro; arreglo estándar: fijar el de intervalo más ancho y re-estimar (`demo_05...m:118-121`; `python/USERGUIDE.md:300-306`).
33. `cluster=true` si se sospechan mínimos múltiples; con split real todo `T` refleja el cluster dominante y `clusterInfo` describe cada cuenca; requiere Statistics Toolbox; con <3 corridas no se evalúa (`ClusteringEvaluated`) (`gsua_userguide.m:114,124`; `CHANGELOG.md:36-38`; `demo_05...m:102-104`).
34. `outlier=true` es agresivo con datos limpios y agrupados: inspeccionar `T.Est` (`gsua_userguide.m:112`).
35. Perfil de verosimilitud: más lento, pero mejor cuando hay correlación porque re-ajusta los demás; `margin` crudo (`0.08` = 8 %) desde v1.18 (`gsua_userguide.m:215-216`).
36. `cost_band ≫ cost_data` = firma de confusión de factores → muestreo conjunto; `cost_data < 1` = datos alcanzables; `cost_data<1 && cost_band<1` = convergido (`gsua_userguide.m:175`; `Functions/gsua_covmetric.m:44-48`).
37. "Un buen ajuste no garantiza parámetros confiables" (lección central) ([TB] `course/05.../demo_05_estimacion.m:4-5`; [TB] `README.md:71-77`).

### Ejecución
38. `gsua_csb` **no** dentro de un Live Script; correr en la ventana de comandos (`gsua_userguide.m:156`).
39. Desactivar paralelismo (`'parallel', false`) si no hay Parallel Toolbox; `'plot', false` en lotes (`demo_03...m:75-79`; `gsua_userguide.m:63`).

---

## 4. Errores comunes y trampas
1. **Modelo de 2 entradas** → `gsua_eval` falla; dentro de `gsua_pe`/`gsua_likelihood` el error se traga: `gsua_deval` emite un `warning` y devuelve `Inf` (`Functions/gsua_deval.m:62-65`), y `gsua_pe` reintenta hasta 100 errores (`gsua_pe.m:131,296-300`). Síntoma: costos `Inf`/`ErrorLimit` sin error claro (no ejecutado).
2. **Archivo fuera del path**: `gsua_dataprep` hace `disp(p1,' Is not recognized...')` (`gsua_dataprep.m:58`), que con dos argumentos probablemente falla con un mensaje opaco en vez del aviso previsto (no ejecutado).
3. `out_names` mal ordenados: análisis silenciosamente equivocado (`demo_00...m:236-243`).
4. Salida transpuesta cuadrada (nSalidas = Nt) no se corrige (`gsua_intrp.m:54-55`).
5. `xdata` fuera de `[t0,tf]` → `NaN`; en Python el riesgo inverso (arrancar en `domain[0]`) (`demo_00...m:261-267`). `gsua_covmetric` excluye corridas con no-finitos con warning (`gsua_covmetric.m:76-81`).
6. Devolver `ode45` crudo → error de interpolación (curso: ~2 personas en un pico de ~356 con RelTol 1e-6) (`course/01.../demo_01_introduccion.m:49-53`).
7. Orientaciones: `gsua_eval` quiere **Np × N**, `gsua_dmatrix` da **N × Np** (`demo_00...m:182-183`); `gsua_ia` `T_est` Np × N vs Python N × Np (cuadrado falla en silencio) (`python/USERGUIDE.md:477`).
8. Con factor fijo, `T`, `M` y `T.Est*` solo llevan los libres (índices se corren); usar `gsua_fullpars` para reconstruir (`gsua_userguide.m:230`; `python/USERGUIDE.md:482-485`).
9. Reutilizar una `T` ya pasada por `gsua_sa` para otro método da warning; usar tabla fresca (`demo_04...m:92-93`).
10. `gsua_eval` grafica por defecto (`show=true`) e imprime progreso; `gsua_ua` abre figuras salvo `'plot',false`; `gsua_sa`/`gsua_ua` paralelos por defecto (`gsua_eval.m:38-39`; `python/USERGUIDE.md:491-494,503`).
11. `'margin'` con `lsqc`/`lsqn` no hace nada y no registra `Margin` en `T` → `gsua_noisefloor` sin override falla con `NoMarginAlpha` (`gsua_pe.m:115`; `gsua_noisefloor.m:154-166`).
12. Escala de `res`: con `lsqc` es suma de cuadrados absoluta (dependiente de unidades); con `gsua_costf` es normalizada (<1 = dentro del margen) (`gsua_pe.m:239-276`). Mezclarlos invalida umbrales.
13. Nombres de opción: `gsua_pe 'B'` es la matriz de restricción lineal; `gsua_bootstrap 'B'` es el número de réplicas (`gsua_pe.m:81`; `gsua_bootstrap.m:84`).
14. Usar ruido de conteo (`poisson/quasipoisson/nbglobal`) en filas que no son conteos (fracciones 0–1) dibuja ceros y anula la escala del costo → usar `'gaussian'`/`'none'` por fila (`gsua_bootstrap.m:47-54`).

### Bugs corregidos (usuarios con versiones viejas pueden toparlos)
| Versión | Problema | Fuente |
|---|---|---|
| <1.18 | `gsua_likelihood` tomaba `margin` desplazado en 1 (`1.08` = 8 %); ahora crudo, `≥1` da warning `MarginConvention`. **Rompe** scripts viejos. | `CHANGELOG.md:11-15` |
| <1.18 | Modelo de usuario de 1 entrada que devuelve matriz rompía `gsua_pardeval` con error de `parfor`. | `CHANGELOG.md:19-23` |
| <1.18 | Fila de datos todo-cero/constante: `fmincon` fallaba "objective undefined" 100 veces sin nombrar la fila. | `CHANGELOG.md:50-53` |
| <1.18 | `gsua_ua` con modelos matriz tomaba `xdata` de un punto. | `CHANGELOG.md:46` |
| <1.17.1 | Gráficas de UA transpuestas cuando N == numel(xdata). | `CHANGELOG.md:64-71` |
| <1.17.1 | `.mltbx` limitado a R2025b ("not supported for this release"). | `CHANGELOG.md:72-78` |
| <1.17.1 | `gsua_bootstrap` lento (Display on en cada refit). | `CHANGELOG.md:79-84` |
| <1.17 | `gsua_pe 'psearch'` pasaba `Np` como punto inicial y `'annealing'` omitía `x0`. | `CHANGELOG.md:131-133` |
| <1.17 | `gsua_noisefloor/bootstrap/noisemodel` confundían `len` con número de salidas si `ydata` venía en columna. | `CHANGELOG.md:134-137` |
| <1.17 | Índice de identificabilidad basado en IC de la mediana (factor esparcido por todo el rango salía "identificado"). | `CHANGELOG.md:121-124` |
| <1.15 | `gsua_pe 'save'` era true y sobrescribía `Estimations.mat`. | `gsua_userguide.m:102`; `CHANGELOG.md:55-56` |
| <1.13 | Salida `tiempo × salidas` no se transponía. | `demo_00...m:246-248`; `errores/sir_transpuesto.m:4-6` |
(El CHANGELOG local solo cubre 1.17.0–1.18.0; versiones anteriores no están documentadas allí.)

### Trampas del hands-on / CONTRATO
- Contrato: `pars` en el orden de `catalogo.json` (fuente única de nombres, rangos, unidades y salidas); los modelos validan `numel(pars)` y fallan con mensaje si no coincide ([SM] `EPI/hands_on_epi/CONTRATO.md:15-17,27-29`; `matlab/modelos/sir_discreto.m:40-43`).
- `sir_discreto`/`sir_cuarentena_cl`: el número de parámetros depende de `opt.prob` (exp: 4; G1/G2: 6; G3: 5) (`sir_discreto.m:11-15,30-43`; `sir_cuarentena_cl.m:26-28`).
- Escenarios con dos picos (`dengue_cali_2015_2016`, `covid_colombia_tres_olas`) a propósito: un modelo de un brote no los reproduce ([SM] `hands_on_epi/README.md:82-84`).
- Los rangos de Medellín/Cali/Itagüí se escalaron desde Bello por población: "ELECCIÓN DEL CURSO, no un valor publicado" (`escenarios.json`, campo `notas` de cada escenario dengue).
- Advertencia de transmisión: `handson.m` avisa si el modelo es de otra transmisión que los datos (`matlab/handson.m:36-38`).

---

## 5. Contrato del hands-on y de los retos

### 5.1 `hands_on_epi` ([SM] `EPI/hands_on_epi/`)
- **Estructura:** `datos/` (CSV reales + `FUENTES.md`), `modelos/catalogo_*.json` (fuente única), `matlab/modelos/<id>.m` y `python/modelos/<id>.py` (mismo id), `matlab/handson.m`, `python/handson.py`, `pruebas/` (`CONTRATO.md:11-23`).
- **Interfaz MATLAB:** `sol = <id>(pars, domain, opt)`, Kind 5 vía `gsua_dataprep`; `domain=[t0 tf]`; `sol.x` = `t0:tf` paso 1; `sol.y` Nsalidas × Nt; continuos con `ode45`+`deval` en `t0:tf`, `sol.solver=''` (`CONTRATO.md:6,27-34`).
- **Interfaz Python:** `def <id>(pars, domain, **opt) -> np.ndarray (Nsalidas, Nt)`; `domain` puede ser `[t0,tf]` o vector de tiempos (`CONTRATO.md:35-39`).
- **Filas:** estados primero, luego salidas observables; el catálogo da el índice (base 1) (`CONTRATO.md:40-42`). Esquema del catálogo: `id, nombre, tipo, unidad_tiempo, transmision, parametros[{id,nombre,descripcion,rango,unidad}], condiciones_iniciales, estados, salidas[{id,fila}], opt_defecto, fuente` (`CONTRATO.md:76-94`).
- **Pruebas obligatorias:** (1) nominal sin NaN ni negativos, conserva población; (2) paridad modelos con N=20 vectores; (3) cada modelo se construye y `gsua_pe`/`parameter_estimation` con N=2 termina; (4) paridad de escenarios (`SOLO_PREPARAR=true`, regenerar el JSON al cambiar `escenarios.json`); (5) paridad de preparación de datos (`CONTRATO.md:99-116`).
- **Escenarios** (`escenarios.json`, 18 entradas): campos `id, nombre, transmision, datos, columna, filtro{columna,desde,hasta}, acumular, salidas, opt, rangos{<modelo>:{<factor>:[lo,hi]}}, modelos_sugeridos, notas`. Solo `dengue_bello_acumulado` usa `acumular: true` (`escenarios.json:200`).
- **Recorrido** (`handson.m`): 1 datos (`x = 0:numel(y)-1`) → 2 modelo (`gsua_dataprep` con `'domain',[0 x(end)]`, `'output', fila`, `'opt'`) → 3 nominal → 4 UA (`N_UA=200`, LHS, % de datos dentro de P5–P95) → 5 SA Xiao (`N_SA=200`) → 6 `gsua_pe` lsqc (`N_MULTI` = 10 si ≤4 libres, 60 si más) → 7 **bootstrap de residuos** (`N_BOOT=20`, re-muestreo de residuos con reemplazo, recorte a ≥0, refits con `max(3, floor(N_MULTI/6))` arranques) → 8 `gsua_ia` + IC por percentiles 2.5–97.5 del bootstrap + veredicto por ancho relativo (<25 % sí, >50 % no) + par más correlacionado (`matlab/handson.m:13-23,46-187`). Semilla 2026; `RAPIDO=true` reduce tamaños (`:19-23`).
- **Salida:** JSON en `resultados/<escenario>__<modelo>__matlab.json` con `costo, datos_dentro_banda_ua, sensibilidad_STi, parametros[{parametro, estimacion, ic_inferior, ic_superior, ancho_relativo, indice_gsua, identificable_sugerido}], tamanos, semilla, segundos` + 3 figuras (`handson.m:41,195-201`).
- **Requisito declarado:** GSUA-CSB **v1.13** (`hands_on_epi/README.md:16`; `handson.m:30`).

### 5.2 HandsOn_SIR ([SM] `EPI/HandsOn_SIR/`)
- Modelo `SIRSemilleroSV(pars, domain, in)` con `pars = [S0 I0 R0 beta gamma alpha]` (CI **primero**), `sol.x = t0:tf`, `sol.y` 3 × Nt (`SIRSemilleroSV.m:1-23,33-53`).
- Datos `states_export(curva, tiempo, salida)`, `xdata=0:100`; curvas 1–3 solo I (`num_outs=2`), 4–6 S+I (`[1 2]`, `ydata` transpuesta a salidas × tiempo) (`README.md:38-47`; `plantilla_estudiante.m:106-118`).
- Flujo de la plantilla: `gsua_dataprep` → `gsua_dmatrix(…,100)`+`gsua_eval` (factibilidad) → `gsua_pe` N=30 fmincon → seleccionar `costf < 1.2·costf(1)` → `gsua_ia(…, correction=true, outlier=false, show=true)` → `gsua_dmatrix`+`gsua_eval` sobre los nuevos rangos (`plantilla_estudiante.m:99-157`; `README.md:49-57`).

### 5.3 Retos ([SM] `retos/`)
- **Pregunta por caso y parámetro:** ¿identificable?, estimación, IC (`README.md:3-11`).
- **Archivos:** `reto_K.json` (modelo, rangos, casos, nivel, sello SHA-256 de la verdad), `reto_K_datos.csv` (`caso,t,salida,valor`), `.mat` (`casos.<id>.t/.y/.observadas`, `rangos`, `parametros`), `entrega_plantilla.json`, `modelos/`, `llenar_entrega`, `guardar_entrega`, `calificar` (público) (`README.md:13-23`).
- **Formato de entrega:** un `entrega.json` por equipo; por caso `c1..` y parámetro `{estimacion, ic_inferior, ic_superior, identificable}` + `justificacion` por caso; ids de `reto_K.json`; campos `equipo`, `integrantes` (`README.md:25-39`; `reto_2/entrega_plantilla.json:1-45`).
- **`llenar_entrega(E, caso, R, ...)`** acepta: tabla de `gsua_ia` (usa `Nominal`, `Range`, `index`; veredicto `index < umbral`, `'umbral'` 0.5), matriz **Np × N** de estimaciones (mediana + cuantiles con interpolación lineal tipo numpy, `'nivel'` 0.95; sin veredicto automático), o struct propio; las filas deben ir en el orden de la plantilla y `numel == Np` (`llenar_entrega.m:9-29,50-81,99-111`).
- **Calificación:** puntaje de intervalo de Gneiting–Raftery `IS = (u−l) + (2/α)(l−x)[x<l] + (2/α)(x−u)[x>u]`, dividido por el ancho del rango `W`; intervalo ausente = todo el rango. Veredicto: acierto 0; "no" siendo sí (o sin respuesta) 1; "sí" siendo no con IC que contiene la verdad 2; "sí" siendo no sin contenerla **5**. Total = Σ IS/W + veredicto; gana el menor; desempate: error relativo medio |x̂−x|/W en los identificables (`README.md:76-105`; `calificar.m:13-27,76-97`).
- **Sello:** SHA-256 de la verdad publicado al lanzar; el calificador se niega si no coincide (`README.md:107-118`).
- **Reto 2:** `SIRSemilleroSV`, versión `exponencial-2026-10-25` (entregas previas inválidas), dominio `[0,100]`, nivel 0.95, 6 casos (c1–c3 solo I con ruido 2/5/10 %; c4–c6 S e I con 2/5/8 %), 36 parámetros; envío por WhatsApp directo a la coordinación (`reto_2/README.md:133-167`; `reto_2/reto_2.json`; `reto_2/LEEME_kit.md:5-7,33-36`).

---

## 6. Diferencias MATLAB / Python (breve)
| Tema | MATLAB | Python | Fuente |
|---|---|---|---|
| `domain` del modelo | `[t0 tf]`; `gsua_eval` interpola | la rejilla pedida; nadie interpola | [TB] `course/00.../README.md:26-33` |
| Varias simulaciones | `par` Np × N; salida N × tiempo × salidas | `M` N × Np; salida N × salidas × tiempo | idem |
| Factor fijo | desaparece de `T` | se queda (`model.fixed`); `full_params` | idem; `python/USERGUIDE.md:482-485` |
| Transpuesta | se corrige (salvo cuadrada) | no | `course/00.../README.md:33` |
| Orden de argumentos | `gsua_sa(M,T)` | `sensitivity_analysis(model, M)` | `python/USERGUIDE.md:475` |
| `T_est`/pool | Np × N | N × Np | `python/USERGUIDE.md:477-478` |
| `noise_floor` margin/alpha | recuperados de `T` | hay que pasarlos (`pe.margin`, `pe.alpha`) | `python/USERGUIDE.md:398-400` |
| Solvers | lsqc/lsqn, fmincon, ga, annealing, particle, psearch, surrogate | least_squares, minimize, differential_evolution, dual_annealing | `python/USERGUIDE.md:518-526` |
| Sin equivalente Py | restricciones `A..nonlcon`, `Multistart`, `Show`, `timer`, `gsua_costcutoff` standalone, `brute-force`, `bandwidth` | — | `python/USERGUIDE.md:507-516` |
| Solo Python | `log_scale`, PEtab | — | `python/USERGUIDE.md:178-185,113-152,494` |
| Defaults | `parallel` true en sa/ua; gráficos implícitos | `n_jobs=1`; nunca grafica implícitamente | `python/USERGUIDE.md:491-503` |
| Percentiles | Hazen | lineal NumPy | `python/USERGUIDE.md:535-536` |
| NaN en `rcostf` | no tolerante | tolerante | `python/USERGUIDE.md:537-538` |
| Simulink | sí | no (usar MATLAB Engine) | [TB] `README.md:32-34` |
| Escalado PE | `TypicalX` = nominal | `x_scale='jac'` | `CHANGELOG.md:42-45` |

---

## 7. Contradicciones o huecos
1. **Firma de 2 entradas.** Guía Apéndice C: "`sol = user_function(factors, domain)` (or `(factors, domain, opt)`)" ([TB] `docs/gsua_userguide.m:272`; también `Functions/gsua_userdefined.m:15-18` dice que 1, 2 o 3 entradas valen) vs curso: con 2 entradas `gsua_eval` falla (`course/00.../matlab/demo_00...m:220-224`; `errores/logistica_dos_entradas.m:4-6`) y código Kind 5 siempre llama con 3 (`gsua_deval.m:57-60`; `gsua_eval.m:144-160`). El propio curso nombra la contradicción (`demo_00...m:221-222`). **Preguntar:** ¿se corrige la guía o el código?
2. **Salida cruda de `ode*`.** Guía: "The native output of MATLAB's `ode*`/`dde*` solvers already meets this" (`gsua_userguide.m:272`) vs curso: no devolver los pasos del integrador, usar `deval` (`demo_00...m:284-288`); además el módulo 01 y su modelo devuelven `ode45` crudo (`course/01.../matlab/sir_continuous_userdefined.m:37-41`), reconocido como subóptimo en `demo_01_introduccion.m:49-53`.
3. **`'margin'` "unificado" vs dual-mode.** CHANGELOG 1.18: margin crudo "everywhere" (`CHANGELOG.md:11-13`) vs `gsua_covmetric`/`gsua_rcostf`/`gsua_costfMulti`: "values <1 are treated as a fraction and converted to 1+margin" (≥1 se usa tal cual) (`Functions/gsua_covmetric.m:29-31`; `gsua_userguide.m:177`). En `gsua_pe` un `1.1` significa 110 %; en `gsua_covmetric`, 10 %.
4. **Ejemplo §7.3 de la guía roto.** `gsua_pe(...,'solver','lsqc',...,'margin',0.1)` seguido de `gsua_noisefloor(T,...)` sin override (`gsua_userguide.m:209-211`): con `lsqc` no se registra `Margin` (`gsua_pe.m:115`) y `gsua_noisefloor` lanza `NoMarginAlpha` (`gsua_noisefloor.m:159-166`). Además `res` de lsqc no está en la escala de `gsua_costf`. La guía §5 no dice que `margin` solo aplica a los solvers no-lsq (el código sí: `gsua_pe.m:35-36`).
5. **`gsua_bootstrap` "must carry Margin/Alpha"** (`gsua_bootstrap.m:32-33`) pero el código cae a `margin=0` si falta (`:153-158`). El curso módulo 5 la usa con `lsqc` sin margin (`demo_05...m:97-98`).
6. **Ruido de conteo sobre datos no-conteo.** El módulo 5 genera ruido gaussiano multiplicativo 5 % sobre datos continuos (`demo_05...m:63-66`) y luego usa `'noise','quasipoisson'` (`:97-98`), mientras la ayuda dice que esos modelos "assume integer-count data" (`gsua_bootstrap.m:47-50`). **Preguntar:** ¿cuándo es aceptable quasipoisson con datos no enteros, y cuándo usar `'gaussian'`?
7. **Firma de `gsua_sa` en el README del módulo 04:** `SA = gsua_sa(M, T, Y)` ([TB] `course/04_Analisis_Sensibilidad/README.md:12`) vs real `[T,Y,J] = gsua_sa(M, T, ...)` (`gsua_userguide.m:69`; `functionSignatures.json`).
8. **Nombres de método del módulo 02:** README MATLAB usa `'latin_hypercube'` (`course/02_Muestreo/README.md:6-9,12`) pero MATLAB solo acepta `'LatinHypercube'|'Uniform'|'Sobol'|'Joint'` vía `validatestring` (`Functions/gsua_dmatrix.m:82-84`); el demo usa los nombres correctos (`demo_02_muestreo.m:65-71`). El README lista `joint` pero el demo no lo cubre.
9. **Módulo 05 README** promete profile likelihood y refinamiento de rangos (`course/05.../README.md:8-10,13`); el demo no llama `gsua_likelihood` ni `gsua_oatr`.
10. **Rutas absolutas.** `course/README.md:80-83` afirma "No hay rutas absolutas codificadas" pero `course/01.../matlab/run_discrete_logistic_example.m:12` y `run_sir_continuous_example.m:12` codifican `C:\Users\Lanre\...`; [SM] `HandsOn_SIR/plantilla_estudiante.m:78-79` y `verificacion.m:9-10` codifican `/home/drojasd/...`.
11. **Versión exigida desactualizada.** [SM] `hands_on_epi/README.md:16`, `handson.m:30` y `retos/reto_2/LEEME_kit.md:11` piden **v1.13**; el kit remite al módulo 5 (`LEEME_kit.md:30-32`) que usa `gsua_bootstrap` y `gsua_clustersummary` (nuevos en 1.17, `CHANGELOG.md:93,99`). Con 1.13 tampoco existe `IndexMethod='spread'`.
12. **Bootstrap del hands-on ≠ bootstrap de la toolbox.** `handson.m:152-161` hace bootstrap **de residuos** (no paramétrico, iid, recorte a ≥0) en lugar de `gsua_bootstrap`; aplicado a `dengue_bello_acumulado` (`escenarios.json:200`) re-muestrea residuos sobre una serie **acumulada**, lo que la guía declara incorrecto (`gsua_userguide.m:202`; `gsua_noisefloor.m:57-61`).
13. **HandsOn_SIR contra las reglas actuales:** (a) orden CI-primero `[S0 I0 R0 beta gamma alpha]` (`SIRSemilleroSV.m:11`) vs convención "parámetros primero" (`course/00.../README.md:12`; [SM] `CONTRATO.md:28-29`); (b) selección `costf < 1.2*costf(1)` (`plantilla_estudiante.m:141`; `verificacion.m:30`) = el tipo de regla de escala que §7.3 desaconseja; (c) pasa el **pool multistart** a `gsua_ia` (`plantilla_estudiante.m:151`) contra "feed `gsua_ia` this bootstrap ensemble, not the raw multistart pool" (`gsua_userguide.m:135`); (d) `RESPUESTAS.md:183-191` concluye identificabilidad a partir de `gsua_ia.Range` (IC de la mediana), que el propio semillero advierte que no es IC del parámetro (`retos/README.md:71-74`); (e) `gsua_eval(M',T,'ydata',ydata,'xdata',xdata,'parallel',false)` (`plantilla_estudiante.m:131,157`) usa estilo nombre-valor con una función de argumentos **posicionales** (`gsua_userguide.m:20`; `gsua_eval.m:33-39`) → probable error o asignación equivocada (no ejecutado).
14. **`pars` "vector fila"** en [SM] `CONTRATO.md:28` vs la toolbox entrega columna (`gsua_userdefined.m:181`); inocuo porque los modelos hacen `pars(:).'`, pero el contrato describe mal lo que llega.
15. **Tolerancias ODE:** curso `AbsTol 1e-8` (`course/00.../matlab/sir_continuo.m:36`) vs hands-on `AbsTol 1e-6` (`CONTRATO.md:39`). No es contradicción formal, pero no hay regla única. Rejilla: curso pide rejilla fina (`dt=0.1`) vs contrato `t0:tf` paso 1 (válido solo si `xdata` es entero).
16. **`functionSignatures.json` incompleto:** `gsua_dataprep` solo declara `p1..p4` genéricos (sin los nombre-valor de usuario); `gsua_noisefloor`/`gsua_noisedraw` no listan `'gaussian'/'none'` (correcto para noisefloor, que no los acepta: `gsua_noisefloor.m:121`); `gsua_csb` expone `'protect','breaking','stretch'` que la guía no documenta (`gsua_userguide.m:157-163`).
17. **Guía: `[T,sol] = gsua_dataprep(...)`** (`gsua_userguide.m:12,273`) pero el 2.º output es el handle del solver (`gsua_userdefined.m:44`). Apéndice C dice `factors` "Np x N" y `sol.x` "1 x N" (N con dos significados) (`gsua_userguide.m:272`).
18. **USERGUIDE Python** dice que `gsua_ia` abre cuatro figuras y que `gsua_pe/gsua_ia` imprimen incondicionalmente (`python/USERGUIDE.md:491-493,501`); en código `gsua_ia show` es false por defecto (`gsua_ia.m:180`) y `gsua_pe` tiene `'timer'` (`gsua_pe.m:88`).
19. **Huecos para preguntar a Daniel:**
   - Regla numérica de N para UA/SA/PE (solo hay ejemplos 250–2000 y la heurística de `Si` negativos).
   - Criterio oficial "identificable": umbral de `index` (`llenar_entrega` usa 0.5, `retos/llenar_entrega.m:23`) vs ancho relativo 25/50 % (`CHANGELOG.md:99-102`) vs `|corr|>0.9`; y cómo se define "identificable" en la verdad del reto (no público).
   - ¿El orden "parámetros primero, luego CI" es obligatorio o solo convención del curso? (la toolbox no lo exige; HandsOn_SIR y Apéndice B simbólico usan CI primero: `python/USERGUIDE.md:110`).
   - ¿Qué solver y `margin` recomienda por defecto para estimación con identificabilidad (curso usa `lsqc` sin margin; guía §5.2 y `gsua_bootstrap` usan `fmincon` con `margin 0.3`)?
   - ¿Se actualizará el hands-on/retos de v1.13 a ≥1.17 (bootstrap paramétrico, `IndexMethod='spread'`)?
   - `docs/gsua_userguide.pdf` está sin versionar en el repo (git status `??`); ¿es la versión publicada?

---

## 8. Fuentes leídas
**[TB] Toolbox** (`D:\MATLAB Drive.source-of-truth-20260718T095207\Toolbox`, HEAD 81276f0):
| Archivo | Líneas |
|---|---|
| README.md | 262 |
| CHANGELOG.md | 137 |
| docs/gsua_userguide.m | 295 |
| Functions/resources/functionSignatures.json | 650 |
| python/README.md | 108 |
| python/USERGUIDE.md | 571 |
| course/README.md | 130 |
| course/00_Escribir_Modelos/README.md | 63 |
| course/00_Escribir_Modelos/matlab/demo_00_escribir_modelos.m | 381 |
| course/00_Escribir_Modelos/matlab/logistica_discreta.m | 33 |
| course/00_Escribir_Modelos/matlab/sir_continuo.m | 43 |
| course/00_Escribir_Modelos/matlab/sir_discreto_plantilla.m | 38 |
| course/00_Escribir_Modelos/matlab/verificar_modelo.m | 65 |
| course/00_Escribir_Modelos/matlab/errores/{logistica_dos_entradas, sir_continuo_crudo, sir_transpuesto}.m | 17 / 16 / 10 |
| course/00_Escribir_Modelos/matlab/solucion/sir_discreto_solucion.m | 22 |
| course/00_Escribir_Modelos/python/errores/{logistica_sin_opt, sir_arranca_en_domain0, sir_continuo_crudo, sir_transpuesto}.py | 19 / 17 / 23 / 11 |
| course/00_Escribir_Modelos/python/verificar_modelo.py | 57 |
| course/01_Introduccion_Modelos/README.md + matlab/{demo_01_introduccion, discrete_logistic_userdefined, sir_continuous_userdefined, run_discrete_logistic_example, run_sir_continuous_example}.m | 26 + 215/45/53/63/74 |
| course/02_Muestreo/README.md + matlab/demo_02_muestreo.m | 13 + 132 |
| course/03_Incertidumbre_MonteCarlo/README.md + matlab/demo_03_incertidumbre.m | 12 + 130 |
| course/04_Analisis_Sensibilidad/README.md + matlab/demo_04_sensibilidad.m | 13 + 111 |
| course/05_Estimacion_Identificabilidad/README.md + matlab/demo_05_estimacion.m | 14 + 150 |
| Functions/gsua_userdefined.m (completo), gsua_deval.m (completo), gsua_intrp.m (completo), gsua_dataprep.m (completo) | 275 / 65 / 119 / 66 |
| Functions/gsua_eval.m (1-60, 95-170), gsua_pe.m (1-135, 225-300), gsua_bootstrap.m (1-94 + grep margin), gsua_noisefloor.m (1-176), gsua_covmetric.m (1-80) | 205 / 408 / 283 / 279 / 122 (parcial) |
| Functions/gsua_dmatrix.m, gsua_sa.m, gsua_ua.m, gsua_ia.m, gsua_csb.m (solo bloques de inputParser vía grep) | 268 / 547 / 108 / 516 / 368 (parcial) |
| Examples/user_dependent.m, user_free.m (solo firma) | — |

**No leídos del Toolbox:** demos Python del curso (`demo_0N_*.py`, `.ipynb`), modelos Python del curso salvo `errores/` y `verificar_modelo.py`, `course/00.../pruebas/*`, `doc/GettingStarted.mlx`, `gsua_userguide.mlx`/`.pdf`, el resto de `Functions/*.m` y `Examples/*`.

**[SM] Semillero** (`D:\Other-projects\Semillero Matemáticas`, HEAD 0174225; copia local completa, no hizo falta `gh api`):
| Archivo | Líneas |
|---|---|
| EPI/hands_on_epi/README.md | 157 |
| EPI/hands_on_epi/CONTRATO.md | 116 |
| EPI/hands_on_epi/escenarios.json (estructura completa y `notas` de las 18 entradas, vía script) | 1061 |
| EPI/hands_on_epi/matlab/handson.m | 215 |
| EPI/hands_on_epi/matlab/modelos/{sir_continuo, sir_discreto, sir_cuarentena_cl, dengue_lizarralde, dengue_cl_discreto}.m | 65 / 113 / 117 / 89 / 102 |
| EPI/HandsOn_SIR/{README, README_estudiante, RESPUESTAS}.md | 90 / 56 / 62 |
| EPI/HandsOn_SIR/{SIRSemilleroSV, plantilla_estudiante, verificacion, generar_datos, fig_datasets}.m | 54 / 118 / 49 / 95 / 22 |
| retos/README.md | 127 |
| retos/reto_2/README.md | 46 |
| retos/reto_2/LEEME_kit.md | 39 |
| retos/reto_2/reto_2.json (campos públicos, vía script) | 88 |
| retos/reto_2/entrega_plantilla.json (1-60) | — |
| retos/llenar_entrega.m | 111 |
| retos/calificar.m (cabecera y fórmula vía grep) | 193 (parcial) |

**No leídos del Semillero:** `hands_on_epi/python/*`, `hands_on_epi/modelos/catalogo_*.json`, `hands_on_epi/datos/FUENTES.md`, `hands_on_epi/pruebas/*`, `HandsOn_SIR/*.mlx`, `retos/privado/*` (contiene la verdad del reto; omitido a propósito), `retos/modelos/*`.
