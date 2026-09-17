# Módulo 03: Análisis de Incertidumbre y Filtrado Monte Carlo

En este módulo estudiamos la propagación de incertidumbre desde los parámetros hacia las salidas del modelo, evaluando la alcanzabilidad y el ajuste de bandas de percentiles.

## Contenidos
* Evaluación de ensamble de Monte Carlo sobre la matriz de diseño.
* Cálculo de percentiles $P_5$ y $P_{95}$ y métrica de contención (`coverage_metric`).
* Filtrado de Monte Carlo (MCF) para clasificar combinaciones de parámetros aceptables vs. no aceptables.

## Funciones Clave
* **MATLAB**: `Y = gsua_ua(M, T)`, `gsua_covmetric(Y, ydata)`, `gsua_MCF`
* **Python**: `Y = gc.uncertainty_analysis(model, M)`, `gc.coverage_metric(Y, ydata)`, `gc.monte_carlo_filter`
