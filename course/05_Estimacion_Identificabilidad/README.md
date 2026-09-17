# Módulo 05: Estimación de Parámetros e Identificabilidad Práctica

En este módulo abordamos el ciclo de calibración de modelos y diagnóstico de identificabilidad práctica.

## Contenidos
* Estimación multicomienzo (Multistart Parameter Estimation - PE) para detectar mínimos globales y locales.
* Matriz de correlación entre estimaciones y compensaciones de factores.
* Profile Likelihood para intervalos de confianza prácticos.
* Detección de múltiples mínimos globales y agrupamiento (Clustering IA).
* Refinamiento de rangos y fijación de parámetros no identificables.

## Funciones Clave
* **MATLAB**: `[Tpe, res] = gsua_pe(T, xdata, ydata)`, `gsua_likelihood`, `gsua_ia`
* **Python**: `pe = gc.parameter_estimation(model, xdata, ydata)`, `gc.profile_likelihood`, `gc.identifiability_analysis`
