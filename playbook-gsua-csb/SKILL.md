---
name: gsua-csb
description: >-
  Experticia para trabajar en GSUA-CSB (toolbox MATLAB) y su port Python gsua_csb:
  análisis global de sensibilidad e incertidumbre, identificabilidad por CSB
  (Confidence Subcontour Box), estimación de parámetros y escritura de modelos de
  usuario. Use when working with GSUA-CSB / gsua_csb — gsua_dataprep, gsua_ua, gsua_sa,
  gsua_pe, gsua_bootstrap, gsua_ia, gsua_clustersummary, gsua_csb, gsua_eval; user-defined
  models (Kinds 1–6), factor order, margin convention, practical identifiability,
  sensitivity reliability, or the guided course (modules 00–05).
---

# GSUA-CSB — paquete de experticia

**GSUA-CSB** (Global Sensitivity and Uncertainty Analysis – Confidence Subcontour Box) es una
toolbox de MATLAB con port Python (`gsua_csb`) para análisis de sensibilidad/incertidumbre,
identificabilidad y calibración de modelos dinámicos. Versión actual: **1.19.1 (MATLAB) / 0.5.0
(Python)**. Autor: Daniel Rojas Díaz (EAFIT). Repo: `drojasd/GSUA-CSB`.

## Cuándo usar este paquete
Cualquier trabajo sobre GSUA-CSB / `gsua_csb`: escribir o depurar modelos de usuario, correr el flujo
UA→SA→PE→IA, interpretar identificabilidad, tocar el curso guiado (módulos 00–05), o responder dudas
de convención (orden de factores, `margin`, criterio de identificabilidad, N del muestreo).

## Mapa del paquete
- **`references/lineamientos.md`** — **empieza aquí.** Decisiones ratificadas por Daniel (contrato de
  Kinds, orden CI-primero, convención `margin`, confiabilidad del SA, criterio de identificabilidad,
  solver por defecto) con citas verificadas. Es la fuente de verdad para zanjar ambigüedades.
- `references/api.md` — tarjeta de API (funciones, firmas, argumentos) verbatim.
- `references/usage.md` — tarjeta de uso (flujos de ejemplo) verbatim.
- `references/flujos.md` — flujos de trabajo end-to-end (UA, SA, UA+MCF, PE+IA, multi-salida Kind 6).
- `references/errores.md` — errores típicos y cómo diagnosticarlos.
- `references/metodos_sa.md` — qué estima cada método de `gsua_sa` (Xiao α, Saltelli, Shapley
  independiente/dependiente), por qué un STi estimado puede pasar de 1, y cuál usar.
- `references/estudio_lineamientos.md` — el estudio independiente que originó estas decisiones (contexto).
- `recipes/` — plantillas ejecutables mínimas por tarea (derivadas del curso 00–05).
- `checks/` — verificación ejecutable del contrato de un modelo y del entorno.

## Reglas de oro (no negociables; detalle y citas en `lineamientos.md`)
1. **Contrato Kind 5:** `sol = f(pars, domain, opt)` — **siempre 3 entradas**; `sol.x` (1×Nt),
   `sol.y` (**salidas × tiempo**). Continuo: `ode45` + `deval` en rejilla fina, no pasos crudos.
2. **Orden de factores:** **condiciones iniciales primero, después parámetros** (convención). La regla
   dura es la **consistencia posicional**: fila i de rangos = `pars(i)` = `names{i}`.
3. **`margin`:** crudo (`0.1 = 10 %`) en `gsua_pe`/`gsua_costf`/`gsua_likelihood`; banda dual-mode en
   `gsua_covmetric`/`gsua_rcostf`/`gsua_costfMulti`. Default `margin > 0`.
4. **Flujo:** `gsua_dataprep → gsua_ua → gsua_sa → gsua_pe → gsua_noisefloor → gsua_bootstrap →
   gsua_ia → gsua_clustersummary`. `gsua_eval` es **posicional**.
5. **Identificabilidad práctica:** mínimo único del costo; se lee con la banda `PctWidthRel` de
   `gsua_clustersummary` (<0.25 identificable / >0.5 no / intermedio dudoso). El **índice es un
   ranking**, no un umbral. `|corr|` alto = combinación, no identificabilidad. La estructural es aparte.
6. **Confiabilidad del SA:** `sum(Si)/sum(|Si|) → 1` indica N suficiente; `sum(Si) → 1` es aditividad.
   **No** hay umbral publicado; no uses `sum(Si) > 0.65`. Un **STi estimado > 1** (o Si < 0) es ruido
   del estimador con N pequeño, no un error: sube N. Factores acoplados → `'SensMethod','Shapley'`
   (`metodos_sa.md`).
7. **Solver por defecto:** `fmincon`, con `TypicalX` automático (v1.18).

## Regla operativa
Haz cambios reversibles locales con normalidad. **Nada se publica (push a GitHub, release o sitio) sin
la aprobación explícita de Daniel.** El release sigue `reference_gsua_release_flow` (prj → buildToolbox
→ tag → gh release), bloqueado por el clasificador salvo que Daniel diga publicar.
