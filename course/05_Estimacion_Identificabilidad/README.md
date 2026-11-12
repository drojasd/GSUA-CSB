# Módulo 05: Estimación de Parámetros e Identificabilidad Práctica

En este módulo abordamos el ciclo de calibración de modelos y diagnóstico de identificabilidad práctica.

## Contenidos
* Estimación multicomienzo (Multistart Parameter Estimation - PE) para detectar mínimos globales y locales.
* Matriz de correlación entre estimaciones y compensaciones de factores.
* Bootstrap paramétrico (`gsua_bootstrap`) para una incertidumbre honesta del ensamble.
* Identificabilidad práctica e índice por parámetro; detección de múltiples mínimos globales (agrupamiento) y resumen por cluster (`gsua_clustersummary`).

## Funciones Clave
* **MATLAB**: `[Tpe, res] = gsua_pe(T, xdata, ydata)`, `gsua_bootstrap`, `gsua_ia`, `gsua_clustersummary`
* **Python**: `pe = gc.parameter_estimation(model, xdata, ydata)`, `gc.parametric_bootstrap`, `gc.identifiability_analysis`, `gc.cluster_summary`

> Herramientas relacionadas, no cubiertas en el demo de este módulo: intervalos por perfil de verosimilitud (`gsua_likelihood` / `gc.profile_likelihood`) y refinamiento de rangos (`gsua_oatr`). Ver la guía de usuario (§7.4 y §6).
