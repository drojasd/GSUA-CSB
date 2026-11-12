# Lineamientos GSUA-CSB (ratificados)

Decisiones autoritativas que resuelven las ambigüedades del *estudio de lineamientos*
(`estudio_lineamientos.md`), **ratificadas por Daniel Rojas Díaz** (autor de GSUA-CSB) el
2026-11-12. Aplican a la toolbox MATLAB (`GSUA-CSB`, v1.18.0) y al port Python (`gsua_csb`, v0.4.0).

Cada lineamiento da: la **regla**, el **porqué**, **dónde** vive en el código y la **cita**.
Las cuatro preguntas abiertas del estudio (Q1–Q4) quedan cerradas aquí; las contradicciones del
estudio se resuelven en la §8.

---

## 1. Contrato de los modelos (Kinds)

GSUA-CSB clasifica el modelo en un *Kind* a partir de su firma y su salida:

| Kind | Firma | Qué es |
|---|---|---|
| 1 | — | Simulink |
| 2 | `y = f(par)` (vectorizada, 1 arg) | devuelve una métrica escalar por muestra |
| 3 / 4 | `f(par(1:ni), par(ni+1:end))` | EDO simbólica (llamada de **2 entradas**) |
| 5 | `sol = f(pars, domain, opt)` | **función de usuario con dominio; SIEMPRE 3 entradas** |
| 6 | `y = f(pars)` (1 arg, matriz) | devuelve `[nSeries × nSteps]` (v1.18) |

**Reglas duras del Kind 5 (el del curso):**
- La firma es **siempre de tres argumentos** `sol = f(pars, domain, opt)`, aunque no uses `opt`
  (decláralo con `~`). `gsua_eval` llama al modelo con tres argumentos; una firma de dos falla al
  evaluar aunque `gsua_dataprep` la acepte.
- Salida: `sol.x` (1×Nt, de `t0=domain(1)` a `tf=domain(2)`) y `sol.y` (**nSalidas × tiempo**).
  Un modelo continuo integra con `ode45` y **evalúa con `deval` en una rejilla fina y regular**;
  no devuelvas los pasos crudos del integrador (pocos y desiguales) porque `gsua_eval` interpola
  linealmente entre ellos.
- **Factores fijos:** si un factor tiene rango degenerado `[v v]`, sale de la tabla `T`, pero el
  modelo **sigue recibiendo el vector completo** con el valor fijo reinsertado en su posición
  (verificado empíricamente: con `ranges=[1 1; 2 4; 5 5]` el modelo recibe `[1 3 5]`). Por eso el
  desempaque posicional del modelo nunca se "corre" al fijar un factor, y `gsua_fullpars(T, libres)`
  reconstruye el vector ordenado completo (v1.18).

En Python el contrato es `func(pars, domain, **opt)`, `domain` es la rejilla de puntos pedidos, la
función devuelve `(n_salidas, len(domain))` y **nadie interpola** por ti; `pars` lleva siempre los Np
factores.

---

## 2. Orden de los factores — Q3 (ratificado: **condiciones iniciales primero**)

- **Convención recomendada:** primero las **condiciones iniciales**, después los **parámetros**.
- **Regla DURA (lo único obligatorio):** consistencia posicional — fila *i* de los rangos =
  `pars(i)` = `names{i}` = lo que el modelo desempaqueta en la posición *i*. Cualquier orden funciona
  mientras esos tres concuerden; CI-primero es la **convención de estilo** del curso, no un requisito
  del motor.
- Dónde: `course/00_Escribir_Modelos/README.md` y `demo_00`; los modelos del curso (`logistica_discreta`,
  `sir_continuo`, `discrete_logistic_userdefined`, `sir_continuous_userdefined`, plantillas, soluciones
  y demos 00–05) están reordenados a CI-primero en MATLAB y Python.
- Consecuencia para el semillero: `SIRSemilleroSV.m` con `[S0 I0 R0 beta gamma alpha]` **ya está en el
  orden recomendado** (no reordenar). Un factor `R(t0)=0`/`R0_pop` que es constante puede quedar al
  final como argumento opcional (análogo a `R(t0)=0` en `sir_continuo`).

---

## 3. Convención de `margin` (ítem 3 — ratificado: unificada, solo documentación)

- **Camino de estimación/verosimilitud — RAW:** en `gsua_pe` / `gsua_costf` / `gsua_likelihood` y en
  el port Python, `margin` es la **fracción cruda** (`0.1 = 10 %`; internamente `m+1`). `gsua_likelihood`
  **cambió** en v1.18: antes recibía el valor ya desplazado (`1.08` = 8 %); ahora pasa la fracción cruda
  (`0.08`). Un valor `>= 1` dispara un warning `MarginConvention`.
- **Helpers de banda — dual-mode:** `gsua_covmetric` / `gsua_rcostf` / `gsua_costfMulti` mantienen una
  banda de doble modo — un valor `< 1` es la misma fracción cruda (así `0.1` = 10 % y concuerda con
  `gsua_pe`), mientras que un valor `>= 1` se lee como ya `1+fracción`.
- **Default:** `margin > 0` (típicamente `0.1`). Solo las rutas con función de costo registran
  `Margin` en `T`; `lsqc`/`lsqn` lo ignoran.

---

## 4. Confiabilidad del SA y tamaño de muestra N — Q1 (ratificado)

- **No hay umbral numérico publicado** para N. El estilo del estudio (`sum(Si) > 0.65`) **era
  incorrecto**: mezclaba dos cosas distintas y fijaba un corte sin fuente.
- **Indicador de confiabilidad del estimador:** `sum(Si) / sum(|Si|)`. Cuando se acerca a **1**, N es
  suficiente: los `Si` espurios negativos (típicos de Saltelli/Xiao con N bajo) se cancelan. Procedimiento
  práctico: subir N hasta que el ratio se estabilice cerca de 1.
- **Aditividad del modelo (cosa aparte):** `sum(Si) → 1` significa que los efectos de primer orden
  explican casi toda la varianza (modelo casi aditivo). No confundir con lo anterior.
- Cita: Lizarralde-Bejarano, D.P., Rojas-Díaz, D., *et al.* (2020), **PLOS ONE**,
  [10.1371/journal.pone.0229668](https://doi.org/10.1371/journal.pone.0229668).

---

## 5. Criterio de identificabilidad práctica — Q2 (ratificado)

- **Criterio primario:** la identificabilidad práctica de un factor es que la función de
  costo/verosimilitud tenga un **mínimo único** (bien determinado) en esa dirección. Base teórica:
  profile likelihood (Raue *et al.* 2009).
- **Operacionalización en la toolbox:** `gsua_clustersummary` → **banda de ancho relativo**
  `PctWidthRel` (ancho del intervalo de percentiles 2.5–97.5 del ensamble, relativo al rango original):
  **< 0.25 identificable**, **> 0.5 no identificable**, intermedio **dudoso**. Es la realización por
  **ensamble** del **CSB** (Confidence Subcontour Box).
- **El índice de identificabilidad es un RANKING, no un umbral.** `index = (2*boxin + corrin + extrin)/4`
  ordena los factores de mejor a peor determinados; no existe un corte tipo `index < 0.5` ratificado.
- **`|corr|` alto = combinación / acoplamiento estructural** entre factores (p. ej. `beta` y `S0` vía
  `R0 = beta*S0/gamma`), **no** identificabilidad por sí solo. Los datos identifican la *combinación*,
  no cada factor por separado.
- **Identificabilidad estructural** (¿es teóricamente posible?) es un análisis **separado y previo** a
  la práctica.
- **No** concluir identificabilidad desde el IC de la mediana de `gsua_ia.Range`.
- Citas: CSB — Rojas-Díaz, D., *et al.* (2024), **Computational Statistics**,
  [10.1007/s00180-023-01362-4](https://doi.org/10.1007/s00180-023-01362-4); profile likelihood — Raue, A.,
  *et al.* (2009), **Bioinformatics**, [10.1093/bioinformatics/btp358](https://doi.org/10.1093/bioinformatics/btp358).

---

## 6. Solver y opciones por defecto — Q4 (ratificado)

- **Solver por defecto:** `fmincon`.
- **`margin > 0`** por defecto (p. ej. `0.1`).
- **Escalado automático:** `gsua_pe` rellena `TypicalX` con las magnitudes nominales para los solvers
  de gradiente (`fmincon`/`lsqc`/`lsqn`) cuando el usuario no lo fija (v1.18); Python usa `x_scale='jac'`.
- **Trazabilidad:** `gsua_pe` registra `T.IP<solver>` (punto inicial usado por estimación) y la
  propiedad `Retries<solver>` (reinicios silenciosos).

---

## 7. Flujo recomendado (detalle en `flujos.md`)

`gsua_dataprep` → `gsua_ua` → `gsua_sa` → `gsua_pe` (multistart) → `gsua_noisefloor` (umbral de
aceptación, reemplaza `costf < 1.5·costf(1)`) → `gsua_bootstrap` (ensamble paramétrico; `'cumulative'`
por fila para series acumuladas) → `gsua_ia` → `gsua_clustersummary`. `gsua_eval` usa argumentos
**posicionales**: `gsua_eval(par, T, xdata, ydata, parallel, show, verbose)`.

---

## 8. Contradicciones del estudio en la toolbox — resueltas (nivel documentación)

Todas verificadas contra el código antes de editar; v1.18.0 (HEAD `81276f0`).

| # | Dónde | Resolución |
|---|---|---|
| 1, 2 | `docs/gsua_userguide.m` (Apéndice C) | Modelo tiempo-dependiente requiere **3 entradas** siempre; patrón `ode45`+`deval` en rejilla fina (no pasos crudos). |
| 3 | `CHANGELOG.md` | Convención `margin` unificada (ver §3). **Solo doc**: `covmetric`/`rcostf`/`costfMulti` conservan banda dual-mode. |
| 4 | guía §7.3 | Solo los solvers con función de costo registran `Margin`; ejemplo cambiado a `fmincon`. |
| 5 | `Functions/gsua_bootstrap.m` (ayuda) | Documentado el fallback `margin=0`. Sin cambio de código. |
| 6 | `course/05 .../demo_05_estimacion.m` | Ruido de bootstrap `'quasipoisson'` → `'gaussian'` (datos continuos). |
| 7 | `course/04 .../README.md` | Firma real `[T, Y, J] = gsua_sa(M, T, ...)` (y Python). |
| 8 | `course/02 .../README.md` | Métodos MATLAB capitalizados (`'LatinHypercube'`…); Python en minúsculas. |
| 9 | `course/05 .../README.md` | Contenidos alineados al demo real (`gsua_bootstrap`, `gsua_clustersummary`). |
| 10 | `course/01 .../run_*.m` | Rutas absolutas → relativas (`fileparts(mfilename('fullpath'))`). |
| 16 | guía §6 | Línea de opciones avanzadas de `gsua_csb`. |
| 17 | guía | `[T, sol]` → `[T, solver]` (nombre real de la 2.ª salida de `gsua_dataprep`). |
| 18 | `python/USERGUIDE.md` | `gsua_ia` abre figuras solo con `show=true`; `gsua_pe` con `timer`. |
| Q1 | guía (línea del SA) | `sum(Si)>0.65` reemplazado por el ratio de confiabilidad + aditividad (ver §4). |
| Q3 | `course/00…05` | Curso reordenado a **CI-primero** (ver §2). |
| — | `functionSignatures.json` | `gsua_csb` ya tenía `protect/breaking/switch/stretch`: sin cambio. |

Las 6 que caen en el repo del semillero se enviaron por el buzón de Octopus al
`semillero-dengue-agent` (mensaje `20261112-52b38faf` + seguimiento `20261112-d6a08f40` que retira el
ítem 13a y ratifica Q1–Q4).

---

## 9. Citas verificadas

- **PLOS ONE 2020** — Lizarralde-Bejarano, Rojas-Díaz, *et al.*, *A study of the parameters' sensitivity
  and identifiability…*, [10.1371/journal.pone.0229668](https://doi.org/10.1371/journal.pone.0229668).
- **Computational Statistics 2024** — Rojas-Díaz, *et al.*, Confidence Subcontour Box (CSB),
  [10.1007/s00180-023-01362-4](https://doi.org/10.1007/s00180-023-01362-4).
- **Bioinformatics 2009** — Raue, *et al.*, *Structural and practical identifiability analysis of partially
  observed dynamical models by exploiting the profile likelihood*,
  [10.1093/bioinformatics/btp358](https://doi.org/10.1093/bioinformatics/btp358).
- **IJERPH 2023** — [10.3390/ijerph20186740](https://doi.org/10.3390/ijerph20186740) (identificabilidad
  estructural vs. práctica en contexto epidemiológico).
- **Zenodo** — GSUA-CSB software, [10.5281/zenodo.3383316](https://doi.org/10.5281/zenodo.3383316).
