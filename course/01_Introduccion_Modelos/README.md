# Módulo 01: Introducción a los Modelos Dinámicos

Punto de partida del curso: cómo describir un modelo para la toolbox GSUA-CSB, declarando
explícitamente sus factores (parámetros + condiciones iniciales) y las salidas, y cómo
simularlo y explorar su variabilidad.

> **Antes de este módulo:** el [módulo 00](../00_Escribir_Modelos/) enseña a escribir un modelo
> de usuario paso a paso (firma, orden de los factores, forma de la salida, errores típicos).
> Aquí los dos modelos ya están escritos; la página muestra su código completo al inicio.

## Contenidos
* Modelo en **tiempo discreto**: crecimiento logístico $N(t+1)=N(t)+rN(t)(1-N(t)/K)$.
* Modelo en **tiempo continuo**: sistema epidémico SIR resuelto con `ode45` / `solve_ivp`.
* Simulación nominal, cálculo del número reproductivo básico $R_0$ y plano de fase.
* Ensamble de trayectorias como primera mirada a la incertidumbre.

## Funciones Clave
* **MATLAB**: `[T,~] = gsua_dataprep(func, ranges, 'names', ..., 'out_names', ...)`,
  `y = gsua_eval(par, T, xdata)`, `M = gsua_dmatrix(T, N)`
* **Python**: `model = create_*_model()` (envuelve `gc.UserFunctionModel`),
  `model.evaluate(par, domain)`, `gc.design_matrix(model, n=N)`

## Archivos
* `matlab/demo_01_introduccion.m` — script guía del módulo.
* `matlab/discrete_logistic_userdefined.m`, `matlab/sir_continuous_userdefined.m` — modelos.
* `python/demo_01_introduccion.py` / `.ipynb` y los modelos `*_userdefined.py`.
