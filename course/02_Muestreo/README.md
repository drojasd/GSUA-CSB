# Módulo 02: Muestreo de Espacios de Parámetros

En este módulo abordamos el diseño de experimentos computacionales mediante el muestreo del hipercubo de parámetros declarados en el modelo.

## Contenidos
* Muestreo aleatorio uniforme (`uniform`).
* Muestreo estratificado por Hipercubo Latino (`latin_hypercube` / `LHS`).
* Muestreo por secuencias de baja discrepancia de Sobol (`sobol`).
* Muestreo conjunto preservando correlación (`joint`).

## Funciones Clave
* **MATLAB**: `M = gsua_dmatrix(T, N, 'Method', 'LatinHypercube')`
* **Python**: `M = gc.design_matrix(model, n=N, method="latin_hypercube")`

> Los nombres de método difieren por lenguaje: en MATLAB se escriben con mayúsculas — `'LatinHypercube'`, `'Uniform'`, `'Sobol'`, `'Joint'` (`gsua_dmatrix`); en Python en minúsculas — `latin_hypercube`, `uniform`, `sobol`, `joint`.
