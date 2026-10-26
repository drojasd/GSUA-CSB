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
* **Errores típicos:** (a) firma, (b) nombres de salida que no coinciden con las filas,
  (c) orientación, (d) la rejilla de evaluación y `t0`, (e) devolver los pasos del integrador,
  (f) unidades.
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

## Paridad MATLAB / Python

Las dos versiones tienen las mismas 11 secciones, en el mismo orden, con los mismos modelos,
ejercicio y cuatro figuras (logística, SIR, rejilla de evaluación, pasos del integrador). Donde
la toolbox se comporta distinto en cada lenguaje (firma, orientación, rejilla), las dos páginas
explican los dos comportamientos. Los modelos dan el mismo resultado numérico:

```bash
matlab -batch "cd course/00_Escribir_Modelos/pruebas; generar_paridad"   # escribe los CSV de MATLAB
python course/00_Escribir_Modelos/pruebas/paridad.py                      # compara con Python
```

Tolerancia relativa: 1e-8 en los modelos discretos y 1e-4 en el continuo (`ode45` frente a
`solve_ivp` con las mismas tolerancias). Resultado actual: 0 en los discretos, 2.6e-6 en el SIR.

## Archivos

Cada archivo existe en los dos lenguajes, con el mismo nombre (`.m` / `.py`):

* `demo_00_escribir_modelos` — guion del módulo (se publica como página).
* `logistica_discreta`, `sir_continuo` — los dos modelos del módulo.
* `verificar_modelo` — revisa que un modelo cumpla el contrato (ayuda del curso).
* `sir_discreto_plantilla` — plantilla del ejercicio; `solucion/sir_discreto_solucion` es la solución.
* `errores/` — modelos escritos mal a propósito: `sir_transpuesto` y `sir_continuo_crudo` en los
  dos lenguajes; el error de firma de cada lenguaje (`logistica_dos_entradas.m`,
  `logistica_sin_opt.py`) y el de `t0` de Python (`sir_arranca_en_domain0.py`), que no tiene
  equivalente en MATLAB.
* `python/mostrar_codigo.py` — muestra un archivo en la página (lo que hace `<include>` en MATLAB).
* `pruebas/` — prueba de paridad numérica (`generar_paridad.m`, `paridad.py` y sus CSV).
