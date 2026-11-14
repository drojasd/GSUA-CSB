# Métodos de sensibilidad de `gsua_sa` — qué estima cada uno y cuándo usarlo

Estudio del 2026-11-14 (v1.19.0 / 0.5.0), verificado contra el código (`Functions/gsua_sa.m`) y la
función de Ishigami (a=7, b=0.1; índices analíticos conocidos). Complementa `lineamientos.md` §4.

## Qué estima cada método

| `SensMethod` | Estimando de `Si_vec`/`STi_vec` | Costo (simulaciones) | Nota |
|---|---|---|---|
| `'Xiao'` (default) | `STi = E‖YA−YABi‖^α / E‖YA−YB‖^α`, α = `'xAlpha'` | `N·(Np/2+1)` | Componentes de distancia / energy distance (Xiao, Lu & Wang 2018). **α=2 ≡ Jansen** (basado en varianza); **α=1 (default) es moment-independent**: más robusto, pero puede ordenar los factores distinto de Sobol' (en Ishigami pone x2 por encima de x1). |
| `'Sobol'` | Índices de Sobol' de `Y` | `N·(Np+1)` | Estimador clásico, insesgado. |
| `'Jansen'` | Índices de Sobol' de `Y` | `N·(Np/2+1)` | Insesgado; = Xiao con α=2. |
| `'Saltelli'` | Índices de Sobol' del **campo `(Y−ynom)²`** | `N·(Np/2+1)` | No mide `Y` sino el desvío cuadrático respecto a la trayectoria nominal. **No es un error**: converge exactamente a los índices de ese campo (Ishigami: `[0.585 0.440 0.423]` vs. verdad del campo `Y²` `[0.587 0.432 0.421]`). Otra pregunta, no documentada antes de v1.19. |
| `'Shapley'` | Efectos de Shapley de `Y` (`T.Sh_vec`) | ver abajo | v1.19. Escribe `T.Sh`/`T.Sh_vec`, no `Si`/`STi`. |
| `'brute-force'`, `'OAT'` | — | — | `'OAT'` está **roto** en la ruta de un objetivo (ver `errores.md` §13c). |

Para **todos** los métodos, los escalares `T.Si`/`T.STi`/`T.Sh` son índices del costo
`J = Σ(Y−ynom)²` (tiempo y series sumados), no de `Y`.

## STi > 1 (y Si < 0): ruido del estimador, no error

Los índices están en [0,1] por definición; sus **estimaciones** Monte Carlo no. El código **no recorta**
a propósito (recortar escondería que N es insuficiente). Frecuencia en Ishigami, 300 réplicas:

| N | Saltelli | Sobol | Jansen | Xiao α=1 | Xiao α=2 |
|---|---|---|---|---|---|
| 64 | P(algún STi>1)=0.25, máx 2.86 | 0.09, máx 1.60 | 0.06, máx 1.82 | **0.00**, máx 1.07 | 0.06, máx 1.85 |
| 256 | 0.04 | 0.01 | 0.00 | 0.00 | 0.00 |
| ≥1024 | 0 | 0 | 0 | 0 | 0 |

**Qué hacer:** subir N hasta que `sum(Si)/sum(|Si|)` se estabilice cerca de 1 (`lineamientos.md` §4).

## Efectos de Shapley (`'SensMethod','Shapley'`)

`Sh_i` reparte la varianza de la salida entre los factores dividiendo **cada interacción por partes
iguales** entre los factores que la forman (Owen 2014; Song, Nelson & Staum 2016). Propiedades:
`Σ Sh_i = 1` (exacto, por construcción), `Si ≤ Sh_i ≤ STi` con factores independientes, y — lo que
importa aquí — **siguen teniendo sentido con factores dependientes** (correlacionados o acoplados
estructuralmente, p. ej. `R0 = β·S0/γ`), donde los índices de Sobol' dejan de ser interpretables.

| `'ShapleyType'` | Muestra | Estimador | Costo | Úsalo cuando |
|---|---|---|---|---|
| `'Independent'` (default) | factores independientes de `T.Range` (como `gsua_dmatrix`) | MC anidado sobre permutaciones (Song, Nelson & Staum 2016) | `N + m·(Np−1)·No·Ni` (defaults `m=N`, `No=1`, `Ni=3`) | Quieres el reparto de interacciones sobre los rangos a priori. |
| `'Given'` | muestra **conjunta**: filas de `M` o `'Pool'` (`Np × nPool`, p. ej. `Tia.Est`) | vecinos más cercanos (Broto, Bachoc & Depecker 2020), `'ShapleyK'`=3 | una simulación por vector distinto | Los factores están acoplados: quieres saber cómo se reparte la varianza **dentro** de esa distribución conjunta. |

**Cuidado con `'Given'`:**
- Con un ensamble de estimaciones aceptadas (`Tia.Est`) la pregunta cambia a *sensibilidad dentro de la
  región identificada*, no sobre los rangos a priori.
- Los vectores repetidos se descartan. `gsua_dmatrix(...,'Method','Joint')` con el `'JointType'` por
  defecto (`'Bootstrap'`) solo repite vectores del pool → no aporta más que el pool. Usa
  `'JointType','SmoothBootstrap'` (o `'Gaussian'`) para tener puntos nuevos.
- Usa unos cientos de vectores distintos o más (avisa con menos de 100).

**Validación** (MATLAB `tests/gsuaShapleyTest.m`, Python `tests/test_sensitivity.py`):
Ishigami → `[0.4358 0.4424 0.1218]` exacto (ambos estimadores convergen a él); `y = x1 + 2·x2` con
`corr = 0.5` → `[0.3393 0.6607]` exacto con `'Given'`, frente a `[0.2 0.8]` con factores
independientes: **la dependencia, no solo el modelo, decide la atribución**.

**Ruido:** un `Sh_i` estimado puede salir ligeramente negativo con N pequeño (mismo fenómeno que
`Si < 0`). El `T.Sh` escalar, construido sobre el costo `J` (de cola pesada), necesita más N que
`T.Sh_vec` para estabilizarse.

## Cuál usar

- **Ranking robusto de factores en un modelo temporal:** `'Xiao'` (default). Si necesitas fracciones
  de varianza comparables con Sobol', `'xAlpha',2` (o `'Jansen'`).
- **Fracciones de varianza clásicas (Si, STi):** `'Jansen'` o `'Sobol'`.
- **¿Qué factores mueven el desvío cuadrático respecto al nominal?** `'Saltelli'`.
- **Factores acoplados o reparto justo de interacciones:** `'Shapley'` — `'Independent'` sobre los
  rangos, `'Given'` sobre un ensamble.

## Citas (verificadas en Crossref)

- Xiao, Lu & Wang (2018), *Risk Analysis* 38(12):2703–2721, [10.1111/risa.13133](https://doi.org/10.1111/risa.13133)
  (componentes de distancia; la que implementa `'Xiao'`). Compañera para salida dinámica: *Struct.
  Multidiscip. Optim.* 57(1):279–291 (en línea 2017), [10.1007/s00158-017-1753-x](https://doi.org/10.1007/s00158-017-1753-x).
- Saltelli *et al.* (2010), *Comput. Phys. Commun.* 181:259–270, [10.1016/j.cpc.2009.09.018](https://doi.org/10.1016/j.cpc.2009.09.018).
- Owen (2014), *SIAM/ASA JUQ* 2(1):245–251, [10.1137/130936233](https://doi.org/10.1137/130936233).
- Song, Nelson & Staum (2016), *SIAM/ASA JUQ* 4(1):1060–1083, [10.1137/15M1048070](https://doi.org/10.1137/15M1048070).
- Broto, Bachoc & Depecker (2020), *SIAM/ASA JUQ* 8(2):693–716, [10.1137/18M1234631](https://doi.org/10.1137/18M1234631).
