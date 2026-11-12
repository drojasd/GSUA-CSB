# Recetas GSUA-CSB

Plantillas ejecutables mínimas por tarea. Las recetas `.m` de aquí asumen que la toolbox
(`Functions/`) y los modelos del curso están en el path (lo resuelven solas). Para flujos completos y
comentados, usa los **demos del curso**, que son las recetas canónicas en los dos lenguajes.

| Tarea | Receta mínima (aquí) | Receta completa (curso, MATLAB + Python) |
|---|---|---|
| Escribir/declarar un modelo | — | `course/00_Escribir_Modelos/{matlab,python}/demo_00_*` |
| Introducción + ensamble | — | `course/01_Introduccion_Modelos/{matlab,python}/demo_01_*` |
| Muestreo | — | `course/02_Muestreo/...` |
| Incertidumbre (UA) + MCF | `recipe_ua_sa.m` | `course/03_Incertidumbre_MonteCarlo/...` |
| Sensibilidad (SA) | `recipe_ua_sa.m` | `course/04_Analisis_Sensibilidad/...` |
| Estimación + identificabilidad | `recipe_pe_ia.m` | `course/05_Estimacion_Identificabilidad/...` |

- **API Python:** nombres y firmas en `../references/usage.md` y `../references/api.md`; recetas Python
  ejecutables = los `demo_0X_*.py` del curso.
- **Convenciones** (orden CI-primero, `margin`, criterio de identificabilidad): `../references/lineamientos.md`.
