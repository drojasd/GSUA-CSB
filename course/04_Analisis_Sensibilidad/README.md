# Módulo 04: Análisis de Sensibilidad Global (GSUA)

Estudio de la importancia relativa de los parámetros sobre la varianza de la salida del modelo.

## Contenidos
* Índices de Sobol de primer orden ($S_i$) y de orden total ($ST_i$).
* Métodos de Jansen, Saltelli y Xiao.
* Análisis One-At-a-Time (OAT).
* Selección de parámetros críticos para estimación y fijación de parámetros no sensibles.

## Funciones Clave
* **MATLAB**: `SA = gsua_sa(M, T, Y)`
* **Python**: `res_sa = gc.sensitivity_analysis(model, M, Y)`
