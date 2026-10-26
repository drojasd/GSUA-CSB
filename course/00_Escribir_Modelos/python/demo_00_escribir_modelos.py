# %% [markdown]
# # Módulo 00: Cómo se escribe un modelo para GSUA-CSB (Python)
# *Curso guiado de GSUA-CSB*
#
# Todo el curso (muestreo, incertidumbre, sensibilidad, estimación) trabaja sobre un
# **modelo de usuario**: una función que tú escribes. La toolbox no ve tus ecuaciones; solo
# llama a esa función con un vector de factores y lee lo que devuelve. Este módulo enseña ese
# contrato paso a paso, en espejo con la versión MATLAB (mismas secciones, mismos modelos,
# mismos resultados):
#
# 1. la firma del modelo y la forma de su salida;
# 2. el orden de los factores: primero los parámetros, después las condiciones iniciales;
# 3. un modelo discreto (un bucle) y uno continuo (un integrador de EDO evaluado en una rejilla);
# 4. cómo se declaran los rangos, los nombres de los factores y de las salidas, qué guarda la
#    toolbox del modelo y cómo se evalúa;
# 5. los errores típicos (firma, nombres de salidas, orientación, rejilla, pasos del
#    integrador, unidades);
# 6. un ejercicio: programar un SIR discreto desde cero.
#
# El código completo de cada modelo aparece en la página tal como está en su archivo.

# %% [markdown]
# ## 1. Configuración del entorno

# %%
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from gsua_csb import UserFunctionModel, design_matrix

# Los modelos de este módulo viven junto a este script (y los que tienen errores a
# propósito, en errores/).
aqui = os.path.dirname(os.path.abspath(__file__))
for carpeta in (aqui, os.path.join(aqui, "errores")):
    if carpeta not in sys.path:
        sys.path.insert(0, carpeta)

from mostrar_codigo import mostrar_codigo
from verificar_modelo import verificar_modelo
from logistica_discreta import logistica_discreta
from sir_continuo import sir_continuo


def ruta(nombre):
    return os.path.join(aqui, nombre)

# %% [markdown]
# ## 2. El contrato de un modelo de usuario
#
# En Python un modelo de usuario es una función (en su propio archivo `.py`, para poder
# importarla desde cualquier guion) con esta firma:
#
# ```python
# def mi_modelo(pars, domain, **opt) -> np.ndarray:
# ```
#
# **Entradas**
#
# * `pars`: arreglo con los factores. La posición manda: `pars[0]` es el factor de la fila 0 de
#   los rangos, `pars[1]` el de la fila 1, etc. Convención del curso: **primero los
#   parámetros, después las condiciones iniciales**.
# * `domain`: **la rejilla de puntos** donde quieres la salida. Es `model.domain` por defecto, o
#   el `xdata` que pases a `model.evaluate(pars, xdata)`. Como puede empezar en cualquier punto,
#   la condición inicial vale en un `t0` fijo que el modelo conoce por su cuenta.
# * `**opt`: opciones del modelo que no son factores; llega como argumentos con nombre lo que
#   pongas en `UserFunctionModel(..., opt={...})`. Si das `opt`, la función debe tener dónde
#   recibirlo (un argumento con ese nombre, o `**opt`).
#
# **Salida**: un arreglo de NumPy de forma **(n_salidas, len(domain))**: una fila por salida,
# en el mismo orden de los nombres de salida, evaluada *exactamente* en los puntos de `domain`.
# Con una sola salida también vale un vector `(len(domain),)`.
#
# **En MATLAB** el contrato es `sol = f(pars, domain, opt)`, pero `domain` es `[t0, tf]`, la
# función devuelve `sol.x` (tiempos) y `sol.y` (salidas x tiempo) y `gsua_eval` interpola
# linealmente en los tiempos pedidos.

# %% [markdown]
# ## 3. Paso 1 — un modelo discreto: la logística
#
# $N(t+1) = N(t) + r\,N(t)\,(1 - N(t)/K)$, con tres factores: `r`, `K` (parámetros) y `N0`
# (condición inicial). Un bucle avanza un paso de tiempo por iteración. Este es el archivo
# `logistica_discreta.py` completo:

# %%
mostrar_codigo(ruta("logistica_discreta.py"))

# %% [markdown]
# Antes de dárselo a la toolbox, llámalo tú a mano y mira las formas:

# %%
y = logistica_discreta([0.3, 500, 10], np.arange(0, 51))
print(f"salida: {y.shape}   (1 salida x 51 tiempos)")
print(f"N(0) = {y[0]:.1f},  N(50) = {y[-1]:.1f}")

# %% [markdown]
# ## 4. Paso 2 — declararlo ante la toolbox: rangos, nombres de factores y de salidas
#
# `UserFunctionModel` recibe la función misma (no su nombre) y una matriz de rangos `(Np, 2)`
# con `[mínimo, máximo]` por fila. La fila *i* corresponde a `pars[i]` y a `names[i]`: las tres
# listas van en el mismo orden. `output_names` nombra las filas de la salida, también en orden.
# `domain` es la rejilla por defecto.

# %%
rangos_log = np.array([
    [0.05, 0.50],    # r   (parámetro)
    [100., 1000.],   # K   (parámetro)
    [5.,   50.],     # N0  (condición inicial, al final)
])

model_log = UserFunctionModel(
    func=logistica_discreta,
    names=["r", "K", "N0"],
    range=rangos_log,
    domain=np.arange(0, 51),
    output_names=["Poblacion"],
)

# %% [markdown]
# El resultado `model` es un *objeto*: guarda cada factor con su rango y su valor nominal (por
# defecto, el punto medio del rango), la rejilla y las salidas; es el equivalente de la tabla
# `T` de MATLAB. Todas las funciones del curso (`design_matrix`, `uncertainty_analysis`,
# `sensitivity_analysis`, `parameter_estimation`...) leen el modelo desde aquí.

# %%
for nombre, (lo, hi), nom in zip(model_log.names, model_log.range, model_log.nominal):
    print(f"{nombre:>4}:  rango [{lo:g}, {hi:g}]   nominal {nom:g}")
print(f"domain : {model_log.domain[0]:g} ... {model_log.domain[-1]:g}  ({len(model_log.domain)} puntos)")
print(f"output_names: {', '.join(model_log.output_names)}   (nombres de salida)")
print(f"activas     : {', '.join(model_log.active_output_names)}   (salidas activas)")

# %% [markdown]
# ## 5. Paso 3 — evaluarlo
#
# `model.evaluate(pars, xdata)` hace **una** simulación (`pars` con los Np factores); `xdata` es
# opcional y por defecto es `model.domain`. `model.evaluate_batch(M)` hace una por **fila** de
# `M`, de forma `(N, Np)`, y devuelve `(N, salidas, tiempo)`.
#
# **En MATLAB** es al revés: `par` es Np x N (una columna por simulación) y la salida de
# `gsua_eval` es N x tiempo x salidas.

# %%
xdata = np.arange(0, 51)
y_nom = model_log.evaluate(model_log.nominal, xdata)
print(f"una simulación   : {y_nom.shape}")

M3 = np.array([model_log.nominal, [0.45, 800, 10], [0.10, 300, 40]])   # 3 filas = 3 simulaciones
y3 = model_log.evaluate_batch(M3, xdata)
print(f"tres simulaciones: {y3.shape}")

plt.figure(figsize=(8, 4))
for fila, etiqueta in zip(y3, ["nominal", "r=0.45, K=800, N0=10", "r=0.10, K=300, N0=40"]):
    plt.plot(xdata, fila, linewidth=1.5, label=etiqueta)
plt.xlabel("Paso de tiempo t"); plt.ylabel("Población N(t)"); plt.grid(True)
plt.legend(loc="lower right"); plt.title("Cada simulación es un vector de factores")
plt.tight_layout(); plt.savefig("demo_00_logistica.png")

# %% [markdown]
# ## 6. Paso 4 — un modelo continuo: el SIR
#
# Un modelo continuo resuelve sus EDO **dentro** de la función. El patrón es:
#
# 1. integrar con `solve_ivp` desde `t0` (donde vale la condición inicial) hasta el último punto pedido;
# 2. pedirle la solución exactamente en la rejilla con `t_eval=domain`;
# 3. devolver `sol.y`, que ya viene como (estados, tiempo).
#
# Factores: `beta`, `gamma` (parámetros), `S0`, `I0` (condiciones iniciales). Salidas: S, I y R,
# en ese orden de filas. Archivo `sir_continuo.py` completo:

# %%
mostrar_codigo(ruta("sir_continuo.py"))

# %%
y = sir_continuo([4.5e-4, 0.125, 950, 25], np.arange(0, 80.5, 0.5))
print(f"salida: {y.shape}  (S, I, R x tiempo)")

# %% [markdown]
# Al declararlo, `output_names` tiene tres nombres porque la salida tiene tres filas. Lo que
# pongas en `opt` llega como argumento con nombre: aquí, el tiempo `t0` de la condición inicial.

# %%
rangos_sir = np.array([
    [1e-4, 8e-4],    # beta  (1/(persona*día))
    [0.05, 0.20],    # gamma (1/día)
    [900., 999.],    # S0    (personas)
    [1.,   50.],     # I0    (personas)
])

model_sir = UserFunctionModel(
    func=sir_continuo,
    names=["beta", "gamma", "S0", "I0"],
    range=rangos_sir,
    domain=np.linspace(0, 80, 161),
    output_names=["S", "I", "R"],
    opt={"t0": 0.0},
)

xs = np.linspace(0, 80, 161)
y_sir = model_sir.evaluate(model_sir.nominal, xs)
print(f"tres salidas activas: {y_sir.shape}")

M = design_matrix(model_sir, n=5, method="latin_hypercube", seed=1)   # 5 x Np: una fila por muestra
Y5 = model_sir.evaluate_batch(M, xs)
print(f"design_matrix: {M.shape}  ->  evaluate_batch: {Y5.shape}")

# %% [markdown]
# Para quedarte con una sola salida (por ejemplo, la que tienes medida), usa `set_output`.
# Afecta a todas las funciones que reciban el modelo.

# %%
model_I = UserFunctionModel(func=sir_continuo, names=model_sir.names, range=rangos_sir,
                            domain=xs, output_names=["S", "I", "R"], opt={"t0": 0.0})
model_I.set_output("I")                            # fila 1 de la salida = I
y_I = model_I.evaluate(model_I.nominal, xs).ravel()
print(f"solo I: {model_I.evaluate(model_I.nominal, xs).shape},  pico de I = {y_I.max():.1f} personas")

plt.figure(figsize=(8, 4))
for fila, etiqueta in zip(y_sir, model_sir.output_names):
    plt.plot(xs, fila, linewidth=1.8, label=etiqueta)
plt.xlabel("Tiempo (días)"); plt.ylabel("Personas"); plt.grid(True)
plt.legend(loc="center right"); plt.title("Las tres salidas del SIR continuo")
plt.tight_layout(); plt.savefig("demo_00_sir.png")

# %% [markdown]
# ## 7. Fijar un factor: rango degenerado
#
# Si un factor no te interesa, dale un rango con mínimo = máximo (o usa `model.fix`). Queda fuera
# del muestreo y de la estimación, pero el modelo lo sigue recibiendo en su posición: `pars`
# siempre llega completo. En Python el factor sigue en el modelo y `pars` lleva siempre los Np
# factores; `design_matrix` mantiene constante la columna del fijo.
#
# **En MATLAB** el factor fijo desaparece de la tabla `T`, y `par` en `gsua_eval` lleva solo los
# factores *libres*.

# %%
rangos_fijo = rangos_sir.copy()
rangos_fijo[3] = [10., 10.]                        # I0 fijo en 10
model_fijo = UserFunctionModel(func=sir_continuo, names=["beta", "gamma", "S0", "I0"],
                               range=rangos_fijo, domain=xs, output_names=["S", "I", "R"])
print("fijos:", {n: bool(f) for n, f in zip(model_fijo.names, model_fijo.fixed)})
M_fijo = design_matrix(model_fijo, n=4, method="latin_hypercube", seed=1)
print(f"design_matrix: {M_fijo.shape} (todos los factores; columna I0 = {M_fijo[:, 3]})")
y_fijo = model_fijo.evaluate(model_fijo.nominal, xs)
print(f"el modelo recibió I0 = 10 en pars[3]: I(0) = {y_fijo[1, 0]:.1f}")

# %% [markdown]
# ## 8. Errores típicos
#
# **(a) Firma.** En Python, `UserFunctionModel` llama `func(pars, domain, **opt)`. Si das
# `opt={"t0": 0}` a una función que no tiene argumento para recibirlo, Python falla con
# "unexpected keyword argument". Declara `t0=...` o `**opt`. (En MATLAB el error equivalente es
# declarar el modelo con dos entradas: `gsua_eval` siempre lo llama con tres.)

# %%
mostrar_codigo(ruta("errores/logistica_sin_opt.py"))

# %%
from logistica_sin_opt import logistica_sin_opt

m_mal = UserFunctionModel(func=logistica_sin_opt, names=["r", "K", "N0"], range=rangos_log,
                          domain=xdata, opt={"t0": 0})
try:
    m_mal.evaluate(m_mal.nominal, xdata)
except TypeError as err:
    print(f"logistica_sin_opt -> error: {err}")

# %% [markdown]
# **(b) Los nombres de salida no coinciden con las filas.** Los nombres de salida no reordenan
# nada: solo nombran las filas, en orden. Si declaras una sola salida llamada `"I"` para un
# modelo cuya fila 0 es S, la toolbox analiza S con la etiqueta `"I"` y no da ningún error.

# %%
m_etiq = UserFunctionModel(func=sir_continuo, names=["beta", "gamma", "S0", "I0"],
                           range=rangos_sir, domain=xs, output_names=["I"], output="I")
print(f"salida etiquetada 'I': valor inicial {m_etiq.evaluate(m_etiq.nominal, xs)[0, 0]:.1f} (es S0, no I0)")

# %% [markdown]
# **(c) Orientación.** La convención es *salidas x tiempo*. En MATLAB, desde la versión 1.13, si
# devuelves `tiempo x salidas` la toolbox lo transpone; en Python nadie la corrige:
# `set_output` selecciona *filas* y te entrega un instante de tiempo en lugar de una salida. No
# dependas de la corrección: devuelve la salida bien orientada.

# %%
mostrar_codigo(ruta("errores/sir_transpuesto.py"))

# %%
from sir_transpuesto import sir_transpuesto

m_T = UserFunctionModel(func=sir_transpuesto, names=["beta", "gamma", "S0", "I0"],
                        range=rangos_sir, domain=xs, output_names=["S", "I", "R"], output="I")
y_T = m_T.evaluate(m_T.nominal, xs)
print(f"sir_transpuesto con set_output('I'): forma {y_T.shape}  <- debería ser (1, 161); Python no la corrige")

# %% [markdown]
# **(d) La rejilla de evaluación y t0.** La condición inicial pertenece a un `t0` fijo, no al
# primer punto que pidas. En MATLAB `t0 = domain(1)` queda fijado en `T`: puedes pedir `xdata`
# desde el día 20 sin problema, pero los puntos fuera de `[t0, tf]` salen `NaN` porque no hay
# simulación ahí. En Python no hay `NaN` (el modelo integra hasta donde le pidas), y el riesgo es
# el contrario: integrar desde `domain[0]`, que hace "arrancar" la epidemia en el primer punto
# pedido.

# %%
mostrar_codigo(ruta("errores/sir_arranca_en_domain0.py"))

# %%
from sir_arranca_en_domain0 import sir_arranca_en_domain0

x20 = np.arange(20, 101)
p = model_sir.nominal
I20_bien = sir_continuo(p, x20)[1]
I20_mal = sir_arranca_en_domain0(p, x20)[1]
print(f"xdata = 20:100: I(20) = {I20_bien[0]:.1f} con t0 fijo;  {I20_mal[0]:.1f} (= I0) si arranca en domain[0]")

plt.figure(figsize=(8, 4))
plt.plot(xs, y_I, "b-", linewidth=2, label="xdata = 0:80")
plt.plot(x20, I20_bien, "bo", markersize=3, label="xdata = 20:100, t0 fijo en 0")
plt.plot(x20, I20_mal, "r--", linewidth=1.5, label="xdata = 20:100, arranca en domain[0] (mal)")
plt.xlim(0, 100); plt.grid(True)
plt.xlabel("Tiempo (días)"); plt.ylabel("Infectados I(t)"); plt.legend(loc="upper right")
plt.title("La condición inicial pertenece a t0, no al primer punto pedido")
plt.tight_layout(); plt.savefig("demo_00_rejilla.png")

# %% [markdown]
# **(e) Devolver los pasos del integrador.** Sin `t_eval`, `solve_ivp` devuelve sus pasos
# adaptativos (pocos y desiguales) y la salida ya no tiene `len(domain)` columnas; la tentación
# es interpolar linealmente entre ellos, que es justo lo que hace `gsua_eval` en MATLAB con la
# estructura cruda de `ode45`. Este modelo lo hace así:

# %%
mostrar_codigo(ruta("errores/sir_continuo_crudo.py"))

# %%
from sir_continuo_crudo import sir_continuo_crudo, pasos_del_integrador

xf = np.arange(0, 80.05, 0.05)
I_bien = sir_continuo(p, xf)[1]
I_crudo = sir_continuo_crudo(p, xf)[1]
t_pasos, y_pasos = pasos_del_integrador(p, 80)
print(f"solve_ivp crudo: {len(t_pasos)} pasos en 80 días. Error máximo en I: "
      f"{np.max(np.abs(I_bien - I_crudo)):.2f} personas (pico {I_bien.max():.1f})")

plt.figure(figsize=(8, 4))
plt.plot(xf, I_bien, "b-", linewidth=2, label="evaluado en una rejilla fina")
plt.plot(xf, I_crudo, "r--", linewidth=1.2, label="pasos del integrador + interpolación lineal")
plt.plot(t_pasos, y_pasos[1], "ro", markersize=4, label="pasos del integrador")
plt.xlim(0, 40); plt.grid(True)
plt.xlabel("Tiempo (días)"); plt.ylabel("Infectados I(t)"); plt.legend(loc="upper left")
plt.title("Devuelve la solución evaluada en una rejilla, no los pasos del integrador")
plt.tight_layout(); plt.savefig("demo_00_pasos.png")

# %% [markdown]
# **(f) Unidades.** Rangos, condiciones iniciales y rejilla deben estar en las mismas unidades.
# En `sir_continuo`, `beta` va en 1/(persona·día) porque el término es `beta*S*I` con S e I en
# personas. Si escribes S0 e I0 como fracciones de la población (0.95, 0.025) sin cambiar
# `beta`, el número reproductivo $R_0 = \beta S_0/\gamma$ se hunde y la epidemia desaparece; si la
# rejilla está en semanas y `gamma` en 1/día, el tiempo de recuperación queda multiplicado por 7.

# %%
b, g = 4.5e-4, 0.125
print(f"R0 con S0 en personas (950) : {b * 950 / g:.2f}")
print(f"R0 con S0 en fracción (0.95): {b * 0.95 / g:.4f}   <- mismas ecuaciones, otra escala")

# %% [markdown]
# ## 9. Comprobar tu modelo antes de usarlo
#
# `verificar_modelo` (ayuda del curso, no de la toolbox) llama a tu modelo una vez y revisa la
# firma, la forma y orientación de la salida, y que no haya NaN ni Inf. Lo que no puede ver es
# la precisión: el modelo crudo de (e) cumple el contrato. Úsalo cada vez que escribas un modelo
# nuevo. (Esta versión revisa además que el valor en un tiempo no dependa de la rejilla pedida:
# el error (d) de Python. La de MATLAB revisa que la salida cubra el dominio `[t0, tf]`.)

# %%
grid = np.linspace(0, 80, 161)
verificar_modelo(logistica_discreta, [0.3, 500, 10], np.arange(0, 51), ["Poblacion"])
verificar_modelo(sir_continuo, model_sir.nominal, grid, ["S", "I", "R"])
verificar_modelo(logistica_sin_opt, [0.3, 500, 10], np.arange(0, 51), ["Poblacion"], t0=0)
verificar_modelo(sir_transpuesto, model_sir.nominal, grid, ["S", "I", "R"])
verificar_modelo(sir_continuo_crudo, model_sir.nominal, grid, ["S", "I", "R"])

# %% [markdown]
# ## 10. Ejercicio: programa un SIR discreto desde cero
#
# Con un paso de un día, el SIR discreto es
#
# $$S_{t+1} = S_t - \beta S_t I_t,\quad I_{t+1} = I_t + \beta S_t I_t - \gamma I_t,\quad R_{t+1} = R_t + \gamma I_t.$$
#
# 1. Copia `sir_discreto_plantilla.py` como `mi_sir_discreto.py` y cambia el nombre de la función.
# 2. Completa el `TODO` del bucle con las tres ecuaciones.
# 3. Pásalo por `verificar_modelo` con salidas S, I y R.
# 4. Decláralo con los mismos rangos, nombres de factores y de salidas que `sir_continuo`, y la
#    rejilla `np.arange(0, 81)`.
# 5. Evalúa el nominal y grafica I(t) junto a la de `sir_continuo`.
#
# **Preguntas.** ¿Por qué no coinciden las dos curvas de I aunque las ecuaciones "son las
# mismas"? ¿Qué pasa con S si `beta*S0 > 1`? ¿Se conserva S + I + R en tu modelo?
#
# La plantilla (ya cumple el contrato; solo le falta la dinámica):

# %%
mostrar_codigo(ruta("sir_discreto_plantilla.py"))

# %% [markdown]
# La solución está en `solucion/sir_discreto_solucion.py`. Ábrela cuando hayas terminado.

# %%
from sir_discreto_plantilla import sir_discreto_plantilla

verificar_modelo(sir_discreto_plantilla, model_sir.nominal, np.arange(0, 81), ["S", "I", "R"])
m_ej = UserFunctionModel(func=sir_discreto_plantilla, names=["beta", "gamma", "S0", "I0"],
                         range=rangos_sir, domain=np.arange(0, 81),
                         output_names=["S", "I", "R"], output="I")
I_ej = m_ej.evaluate(m_ej.nominal).ravel()
if np.all(np.abs(I_ej - I_ej[0]) < 1e-12):
    print("La plantilla cumple el contrato, pero I(t) sigue constante: completa el TODO.")
else:
    print(f"I(t) ya cambia: pico de {I_ej.max():.1f} personas en el día {int(np.argmax(I_ej))}.")

# %% [markdown]
# ## 11. Resumen: lista de chequeo
#
# * Función en su propio archivo `.py`: `func(pars, domain, **opt)`, con argumento para cada clave de `opt`.
# * `pars[i]` = fila *i* de los rangos = `names[i]`: parámetros primero, condiciones iniciales después.
# * Salida `(n_salidas, len(domain))` evaluada en los puntos de `domain`, filas en el orden de `output_names`; la condición inicial vale en un `t0` fijo.
# * Discreto: un bucle de `t0` al último punto pedido. Continuo: `solve_ivp(..., t_eval=domain)` desde `t0`.
# * `UserFunctionModel(func, names, range, domain, output_names, opt=...)` crea el modelo.
# * `model.evaluate_batch(M, xdata)` con `M` de N x Np devuelve N x salidas x len(xdata).
# * Una sola salida: `model.set_output("I")`. Un factor fijo: rango `[v, v]` o `model.fix("I0")`.
#
# En el módulo 1 estos mismos dos modelos (`discrete_logistic_userdefined.py` y
# `sir_continuous_userdefined.py`) se exploran con ensambles de trayectorias.
