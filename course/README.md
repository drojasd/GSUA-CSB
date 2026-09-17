# Curso guiado de GSUA-CSB

Un curso práctico para aprender a usar la toolbox **GSUA-CSB** (Global Sensitivity and
Uncertainty Analysis — Confidence Sub-contour Box) de principio a fin: desde escribir un
modelo dinámico hasta analizar su sensibilidad, propagar incertidumbre, estimar parámetros
y diagnosticar identificabilidad práctica.

Cada sesión está escrita **dos veces** —una en **MATLAB** y otra en **Python**— resolviendo
lo mismo, para que sigas la que uses. Los scripts son ejecutables, están comentados paso a
paso e incluyen preguntas y ejercicios.

> Esta es la carpeta `course/` del repositorio [drojasd/GSUA-CSB](https://github.com/drojasd/GSUA-CSB).
> Material didáctico desarrollado en la Universidad EAFIT. Documentación de la toolbox:
> <https://drojasd.github.io/GSUA-CSB/>

---

## El flujo que enseña el curso

Los módulos siguen el mismo orden del flujo de trabajo de la toolbox:

```
  Modelo  ─►  Muestreo  ─►  Incertidumbre  ─►  Sensibilidad  ─►  Estimación + Identificabilidad
 (Módulo 1)  (Módulo 2)      (Módulo 3)         (Módulo 4)             (Módulo 5)
```

| # | Módulo | Tema | Funciones clave (MATLAB → Python) |
|---|--------|------|-----------------------------------|
| 01 | `01_Introduccion_Modelos` | Definir modelos discretos y continuos; simulación nominal y ensamble | `gsua_dataprep`, `gsua_eval` → `UserFunctionModel`, `.evaluate` |
| 02 | `02_Muestreo` | Diseño de experimentos: Uniforme, Hipercubo Latino, Sobol; escala log | `gsua_dmatrix` → `gc.design_matrix` |
| 03 | `03_Incertidumbre_MonteCarlo` | Propagación de incertidumbre, percentiles P5/P95, filtrado Monte Carlo | `gsua_ua`, `gsua_covmetric`, `gsua_MCF` → `gc.uncertainty_analysis`, `gc.coverage_metric`, `gc.monte_carlo_filter` |
| 04 | `04_Analisis_Sensibilidad` | Índices de Sobol (Sᵢ, STᵢ), Jansen/Saltelli/Xiao, OAT | `gsua_sa` → `gc.sensitivity_analysis` |
| 05 | `05_Estimacion_Identificabilidad` | Estimación multicomienzo, correlación, identificabilidad por bootstrap | `gsua_pe`, `gsua_ia` → `gc.parameter_estimation`, `gc.identifiability_analysis` |

Los modelos ejemplo son un **crecimiento logístico discreto** y un **modelo epidémico SIR
continuo** (resuelto con `ode45` / `solve_ivp`), que se reutilizan a lo largo de los módulos.

---

## Estructura

```
course/
├── README.md                     ← este archivo (índice del curso)
├── 01_Introduccion_Modelos/
│   ├── README.md
│   ├── matlab/                    ← demo_01_introduccion.m + modelos .m
│   └── python/                    ← demo_01_introduccion.py / .ipynb + modelos .py
├── 02_Muestreo/         { README.md, matlab/, python/ }
├── 03_Incertidumbre_MonteCarlo/  { README.md, matlab/, python/ }
├── 04_Analisis_Sensibilidad/     { README.md, matlab/, python/ }
└── 05_Estimacion_Identificabilidad/ { README.md, matlab/, python/ }
```

Cada módulo trae, dentro de `matlab/` y `python/`, un **demo ejecutable** (`demo_NN_*`) más
los **archivos de modelo** que ese demo utiliza. Los cinco módulos están implementados en
ambos lenguajes.

---

## Requisitos e instalación

### MATLAB

- **MATLAB R2019b o superior.**
- Toolboxes recomendados: **Statistics and Machine Learning** (muestreo LHS/Sobol),
  **Optimization** y **Global Optimization** (estimación del Módulo 05). *Symbolic Math* y/o
  *Simulink* solo si trabajas esos tipos de modelo (no se usan en este curso).
- La toolbox GSUA-CSB debe estar disponible. Tres opciones:
  - **Clonar el repositorio** y ejecutar los demos: cada demo localiza automáticamente la
    carpeta `Functions/` del propio repositorio (ruta relativa), así que funciona sin instalar nada.
  - **Add-On Explorer / File Exchange:** busca *"GSUA-CSB"*.
  - **Paquete:** doble clic sobre `GSUA-CSB.mltbx`.

> Los demos MATLAB **localizan la toolbox de forma portable**: si sus funciones ya están en
> el path no hacen nada; si no, buscan primero dentro del repositorio clonado y luego en
> ubicaciones habituales; si no la encuentran, se detienen con un mensaje claro. No hay rutas
> absolutas codificadas.

### Python

- **Python 3.10 o superior.**
- Instala el paquete `gsua-csb` en modo editable desde la carpeta `python/` de la toolbox:

```bash
cd python          # carpeta python/ del repositorio
pip install -e ".[all]"
```

Esto instala `numpy`, `scipy`, `scikit-learn`, `matplotlib` (y `sympy` para modelos
simbólicos). Luego:

```python
import gsua_csb as gc   # gc.UserFunctionModel, gc.design_matrix, gc.sensitivity_analysis, ...
```

---

## Cómo correr cada sesión

### En MATLAB

1. Abre MATLAB (con la toolbox disponible, ver arriba).
2. Navega a la carpeta `matlab/` del módulo, por ejemplo `course/01_Introduccion_Modelos/matlab`.
3. Ejecuta el demo, p. ej. `demo_01_introduccion`. El script agrega su propio modelo al path
   (ruta relativa) y genera las figuras.

### En Python

1. Activa el entorno donde instalaste `gsua-csb`.
2. Desde la carpeta `python/` del módulo:

   ```bash
   python demo_01_introduccion.py
   ```

   o abre el notebook `demo_01_introduccion.ipynb` en Jupyter. Las figuras se guardan como
   `.png` en la misma carpeta.

---

## Versiones en línea

Cada módulo se puede leer ejecutado (sin instalar nada) en el sitio de documentación:
**[drojasd.github.io/GSUA-CSB](https://drojasd.github.io/GSUA-CSB/)**.
