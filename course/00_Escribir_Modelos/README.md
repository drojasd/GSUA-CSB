# Módulo 00: Cómo se escribe un modelo para GSUA-CSB

El resto del curso da por hecho que el modelo ya está programado. Este módulo enseña a
programarlo: qué firma debe tener la función, en qué orden van los factores, qué forma tiene
la salida y cómo se declara ante la toolbox. Va antes del módulo 1 porque todo lo demás
(muestreo, incertidumbre, sensibilidad, estimación) llama a esa función.

## Contenidos
* **El contrato.** MATLAB: `sol = f(pars, domain, opt)` en su propio archivo `.m`, con `sol.x`
  (1 × Nt, de `t0` a `tf`) y `sol.y` (nSalidas × Nt). Python: `func(pars, domain, **opt)` que
  devuelve un arreglo (n_salidas, len(domain)) evaluado en los puntos de `domain`.
* **Orden de los factores:** primero los parámetros, después las condiciones iniciales; fila *i*
  de los rangos = `pars(i)` = `names{i}`.
* **Modelo discreto** (logística, con un bucle) y **modelo continuo** (SIR con `ode45` + `deval`
  en una rejilla / `solve_ivp` con `t_eval`).
* **Declarar y evaluar:** `gsua_dataprep` (rangos, `names`, `out_names`, `opt`), la tabla `T` y
  sus propiedades, `gsua_eval`; en Python `UserFunctionModel`, `evaluate`, `evaluate_batch`,
  `set_output`. Factores fijos con rango degenerado.
* **Errores típicos:** firma con dos entradas (MATLAB), `opt` sin argumento que lo reciba
  (Python), `out_names` que no coincide con las filas, orientación de la salida, rejilla (fuera
  del dominio, `ode45` sin `deval`, `t0` atado al primer punto pedido, `solve_ivp` sin `t_eval`)
  y unidades.
* **Ejercicio:** programar un SIR discreto desde una plantilla y comprobarlo con `verificar_modelo`.

## Diferencias MATLAB ↔ Python que conviene tener presentes

| | MATLAB | Python |
|---|---|---|
| `domain` que recibe el modelo | `[t0, tf]` | la rejilla de puntos pedidos |
| Quién interpola | `gsua_eval` (lineal sobre `sol.x`) | nadie: el modelo devuelve los puntos pedidos |
| Varias simulaciones | `par` de Np × N (columnas) | `M` de N × Np (filas) |
| Forma de la salida | N × tiempo × salidas | N × salidas × tiempo |
| Factor fijo | desaparece de `T`; `par` lleva solo los libres | `pars` siempre lleva los Np factores |
| Salida transpuesta | se corrige sola (salvo caso cuadrado) | no se corrige |

## Archivos
* `matlab/demo_00_escribir_modelos.m` — script guía (se publica como página).
* `matlab/logistica_discreta.m`, `matlab/sir_continuo.m` — los dos modelos del módulo.
* `matlab/errores/` — modelos escritos mal a propósito para la sección de errores.
* `matlab/verificar_modelo.m` — revisa que un modelo cumpla el contrato (ayuda del curso).
* `matlab/sir_discreto_plantilla.m` — plantilla del ejercicio; `matlab/solucion/` tiene la solución.
* `python/demo_00_escribir_modelos.py` — script guía (los modelos están definidos dentro).
* `python/verificar_modelo.py`, `python/sir_discreto_plantilla.py`, `python/solucion/`.
