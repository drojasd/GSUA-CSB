# Módulo 04: Análisis de Sensibilidad Global (GSUA)

Estudio de la importancia relativa de los parámetros sobre la varianza de la salida del modelo.

## Contenidos
* Índices de Sobol de primer orden ($S_i$) y de orden total ($ST_i$).
* Métodos de Jansen, Saltelli y Xiao.
* Análisis One-At-a-Time (OAT).
* Selección de parámetros críticos para estimación y fijación de parámetros no sensibles.

## Funciones Clave
* **MATLAB**: `[T, Y, J] = gsua_sa(M, T, ...)` — los índices quedan en `T.Si`, `T.STi`, `T.Si_vec`, `T.STi_vec`; `Y` es la matriz de salida y `J` su costo escalar.
* **Python**: `res_sa = gc.sensitivity_analysis(model, M)`
