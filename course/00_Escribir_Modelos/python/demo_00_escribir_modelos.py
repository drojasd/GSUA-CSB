# %% [markdown]
# # Módulo 00: Cómo se escribe un modelo para GSUA-CSB (Python)
# *Curso guiado de GSUA-CSB*
#
# Todo el curso (muestreo, incertidumbre, sensibilidad, estimación) trabaja sobre un
# **modelo de usuario**: una función de Python que tú escribes. El paquete `gsua_csb` no ve tus
# ecuaciones; solo llama a esa función con un vector de factores y lee lo que devuelve. Este
# módulo enseña ese contrato paso a paso, en espejo con la versión MATLAB:
#
# 1. la firma `func(pars, domain, **opt)`, `UserFunctionModel` y la forma de la salida;
# 2. el orden de los factores: primero los parámetros, después las condiciones iniciales;
# 3. un modelo discreto (un bucle) y uno continuo (`solve_ivp` evaluado en la rejilla);
# 4. cómo se declaran `range`, `names` y `output_names`, qué hace el objeto `model` (la tabla
#    `T` de MATLAB) y cómo se evalúa con `model.evaluate` (el `gsua_eval` de MATLAB);
# 5. los errores típicos (firma, orientación, rejilla, unidades);
# 6. un ejercicio: programar un SIR discreto desde cero.

# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from gsua_csb import UserFunctionModel, design_matrix
from verificar_modelo import verificar_modelo

# %% [markdown]
# ## 1. El contrato de un modelo de usuario
#
# ```python
# def mi_modelo(pars, domain, **opt) -> np.ndarray:   # forma (n_salidas, len(domain))
# ```
#
# **Entradas**
#
# * `pars`: arreglo con los factores. La posición manda: `pars[0]` es el factor de la fila 0 de
#   `range`, `pars[1]` el de la fila 1, etc. Convención del curso: **primero los parámetros,
#   después las condiciones iniciales**.
# * `domain`: **la rejilla de puntos** donde quieres la salida (no `[t0, tf]` como en MATLAB).
#   Es `model.domain` por defecto, o el `xdata` que pases a `model.evaluate(pars, xdata)`.
# * `**opt`: lo que pongas en `UserFunctionModel(..., opt={...})` llega como argumentos con
#   nombre. Si no usas `opt`, basta con `func(pars, domain)`.
#
# **Salida**: un arreglo de NumPy de forma **(n_salidas, len(domain))**: una fila por salida, en el
# orden de `output_names`, evaluada *exactamente* en los puntos de `domain`. Con una sola salida
# también vale un vector `(len(domain),)`.
#
# **Diferencia clave con MATLAB.** En MATLAB el modelo devuelve `sol.x`/`sol.y` y `gsua_eval`
# interpola. En Python nadie interpola por ti: la función recibe la rejilla y debe devolver los
# valores en esos puntos. Y como `domain` puede empezar en cualquier punto, el modelo tiene que
# saber por su cuenta en qué tiempo `t0` vale la condición inicial.

# %% [markdown]
# ## 2. Paso 1 — un modelo discreto: la logística
#
# $N(t+1) = N(t) + r\,N(t)\,(1 - N(t)/K)$, con tres factores: `r`, `K` (parámetros) y `N0`
# (condición inicial). El bucle avanza un paso por iteración desde `t0` hasta el último punto
# pedido; al final se devuelven los valores en los puntos de `domain`.

# %%
def logistica_discreta(pars, domain, t0=0):
    """N(t+1) = N(t) + r*N(t)*(1 - N(t)/K).  pars = [r, K, N0]; N(t0) = N0."""
    r, K, N0 = pars                                   # parámetros primero, condición inicial al final
    domain = np.asarray(domain, dtype=float)

    pasos = np.arange(t0, int(np.ceil(domain.max())) + 1)   # t0, t0+1, ..., último paso pedido
    N = np.zeros(len(pasos))
    N[0] = N0
    for k in range(len(pasos) - 1):
        N[k + 1] = N[k] + r * N[k] * (1 - N[k] / K)

    return np.interp(domain, pasos, N)               # (len(domain),): una sola salida

# %% [markdown]
# Antes de dárselo al paquete, llámalo tú a mano y mira la forma:

# %%
y = logistica_discreta([0.3, 500, 10], np.arange(0, 51))
print("forma de la salida:", y.shape, " (una salida x 51 tiempos)")
print(f"N(0) = {y[0]:.1f},  N(50) = {y[-1]:.1f}")

# %% [markdown]
# ## 3. Paso 2 — declararlo con `UserFunctionModel`: range, names y output_names
#
# `range` es una matriz `(Np, 2)` con `[mínimo, máximo]` por fila. La fila *i* corresponde a
# `pars[i]` y a `names[i]`: las tres listas van en el mismo orden. `output_names` nombra las filas
# de la salida, también en orden. `domain` es la rejilla por defecto.

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
# El objeto `model` hace el papel de la tabla `T` de MATLAB: guarda nombres, rangos, nominal
# (por defecto, el punto medio), rejilla y salidas. Todas las funciones del curso
# (`design_matrix`, `uncertainty_analysis`, `sensitivity_analysis`, `parameter_estimation`...)
# leen el modelo desde aquí.

# %%
for nombre, (lo, hi), nom in zip(model_log.names, model_log.range, model_log.nominal):
    print(f"{nombre:>4}:  rango [{lo:g}, {hi:g}]   nominal {nom:g}")
print("output_names:", model_log.output_names)
print("domain      :", model_log.domain[0], "...", model_log.domain[-1], f"({len(model_log.domain)} puntos)")

# %% [markdown]
# ## 4. Paso 3 — evaluar con `model.evaluate` y `model.evaluate_batch`
#
# `model.evaluate(pars, xdata)` hace **una** simulación (`pars` con los Np factores); `xdata` es
# opcional y por defecto es `model.domain`. `model.evaluate_batch(M)` hace una por **fila** de `M`
# (forma `(N, Np)`, como la que devuelve `design_matrix`).
#
# Ojo con el espejo: en MATLAB `par` es Np x N (una *columna* por simulación) y la salida de
# `gsua_eval` es N x tiempo x salidas; en Python `M` es N x Np (una *fila* por simulación) y la
# salida de `evaluate_batch` es N x salidas x tiempo.

# %%
y_nom = model_log.evaluate(model_log.nominal)
print("una simulación:", y_nom.shape)

M3 = np.array([model_log.nominal, [0.45, 800, 10], [0.10, 300, 40]])   # 3 filas = 3 simulaciones
Y3 = model_log.evaluate_batch(M3)
print("tres simulaciones:", Y3.shape)

plt.figure(figsize=(8, 4))
for fila, etiqueta in zip(Y3, ["nominal", "r=0.45, K=800, N0=10", "r=0.10, K=300, N0=40"]):
    plt.plot(model_log.domain, fila, linewidth=1.5, label=etiqueta)
plt.xlabel("Paso de tiempo t"); plt.ylabel("Población N(t)"); plt.grid(True)
plt.title("Cada fila de M es una simulación"); plt.legend(loc="lower right")
plt.tight_layout(); plt.savefig("demo_00_logistica.png")

# %% [markdown]
# ## 5. Paso 4 — un modelo continuo: el SIR con `solve_ivp`
#
# Un modelo continuo resuelve sus EDO **dentro** de la función. El patrón es:
#
# 1. integrar con `solve_ivp` desde `t0` (donde vale la condición inicial) hasta el último punto pedido;
# 2. pedirle la solución exactamente en la rejilla con `t_eval=domain`;
# 3. devolver `sol.y`, que ya viene como (estados, tiempo).
#
# Factores: `beta`, `gamma` (parámetros), `S0`, `I0` (condiciones iniciales). Salidas: S, I y R,
# en ese orden de filas.

# %%
def sir_continuo(pars, domain, t0=0.0):
    """SIR continuo. pars = [beta, gamma, S0, I0]; R(t0) = 0. Salida (3, len(domain)) = [S; I; R]."""
    beta, gamma, S0, I0 = pars
    domain = np.asarray(domain, dtype=float)

    def rhs(t, y):
        S, I, R = y
        return [-beta * S * I, beta * S * I - gamma * I, gamma * I]

    sol = solve_ivp(rhs, (t0, domain.max()), [S0, I0, 0.0],
                    t_eval=domain, rtol=1e-6, atol=1e-8)
    return sol.y                                       # (3, len(domain))

# %%
y = sir_continuo([4.5e-4, 0.125, 950, 25], np.linspace(0, 80, 161))
print("forma de la salida:", y.shape, " (S, I, R x tiempo)")

# %% [markdown]
# Al declararlo, `output_names` tiene tres nombres porque la salida tiene tres filas. Lo que pongas
# en `opt` llega como argumento con nombre: aquí fijamos `t0=0.0` de forma explícita.

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

y_sir = model_sir.evaluate(model_sir.nominal)
print("tres salidas activas:", y_sir.shape)

M = design_matrix(model_sir, n=5, method="latin_hypercube", seed=1)
print("design_matrix:", M.shape, "->  evaluate_batch:", model_sir.evaluate_batch(M).shape)

# %% [markdown]
# Para quedarte con una sola salida (la que tienes medida), usa `set_output`. Afecta a todas las
# funciones que reciban el modelo — es el `T.Properties.CustomProperties.output = 2` de MATLAB.

# %%
model_sir.set_output("I")
y_I = model_sir.evaluate(model_sir.nominal)
print("solo I:", y_I.shape, f" pico de I = {y_I.max():.1f} personas")
model_sir.set_output(None)                         # volver a las tres salidas

plt.figure(figsize=(8, 4))
for fila, etiqueta in zip(y_sir, model_sir.output_names):
    plt.plot(model_sir.domain, fila, linewidth=1.8, label=etiqueta)
plt.xlabel("Tiempo (días)"); plt.ylabel("Personas"); plt.grid(True)
plt.title("Las tres salidas: y[j] es la salida j"); plt.legend(loc="center right")
plt.tight_layout(); plt.savefig("demo_00_sir.png")

# %% [markdown]
# ## 6. Fijar un factor: rango degenerado
#
# Si un factor no te interesa, dale un rango con mínimo = máximo (o usa `model.fix`). Queda fuera
# del muestreo y la estimación. A diferencia de MATLAB, donde el factor fijo desaparece de `T` y
# `par` lleva solo los libres, en Python **`pars` siempre lleva los Np factores**: `design_matrix`
# devuelve todas las columnas y mantiene la del factor fijo constante.

# %%
model_fijo = UserFunctionModel(
    func=sir_continuo, names=["beta", "gamma", "S0", "I0"],
    range=np.vstack([rangos_sir[:3], [10., 10.]]),     # I0 fijo en 10
    domain=np.linspace(0, 80, 161), output_names=["S", "I", "R"],
)
print("fijos:", {n: bool(f) for n, f in zip(model_fijo.names, model_fijo.fixed)})
M_fijo = design_matrix(model_fijo, n=4, method="latin_hypercube", seed=1)
print("design_matrix:", M_fijo.shape, " columna I0:", M_fijo[:, 3])

# %% [markdown]
# ## 7. Errores típicos
#
# **(a) Firma: `opt` sin argumento que lo reciba.** Si das `opt={"t0": 0}` pero la función es
# `func(pars, domain)`, Python se queja del argumento inesperado. Declara `t0=...` o `**opt`.

# %%
def logistica_sin_opt(pars, domain):
    return logistica_discreta(pars, domain)

m_mal = UserFunctionModel(func=logistica_sin_opt, names=["r", "K", "N0"], range=rangos_log,
                          domain=np.arange(0, 51), opt={"t0": 0})
try:
    m_mal.evaluate(m_mal.nominal)
except TypeError as e:
    print("TypeError:", e)

# %% [markdown]
# **(b) `output_names` no coincide con las filas.** `output_names` no reordena nada: solo pone
# nombre a las filas, en orden. Si declaras una única salida `"I"` para un modelo cuya fila 0 es S,
# `set_output("I")` te entrega S con la etiqueta `"I"`, sin ningún error.

# %%
m_etiq = UserFunctionModel(func=sir_continuo, names=["beta", "gamma", "S0", "I0"],
                           range=rangos_sir, domain=np.linspace(0, 80, 161), output_names=["I"])
m_etiq.set_output("I")
print(f"salida etiquetada 'I': valor inicial {m_etiq.evaluate(m_etiq.nominal)[0, 0]:.1f} (es S0, no I0)")

# %% [markdown]
# **(c) Orientación.** Si devuelves `sol.y.T` (tiempo, salidas), Python no lo corrige: `set_output`
# selecciona *filas* y te entrega un instante de tiempo en lugar de una salida.

# %%
def sir_transpuesto(pars, domain, t0=0.0):
    return sir_continuo(pars, domain, t0).T          # (len(domain), 3): MAL

m_T = UserFunctionModel(func=sir_transpuesto, names=["beta", "gamma", "S0", "I0"],
                        range=rangos_sir, domain=np.linspace(0, 80, 161),
                        output_names=["S", "I", "R"], output="I")
print("forma con set_output('I'):", m_T.evaluate(m_T.nominal).shape, " <- deberían ser (1, 161)")

# %% [markdown]
# **(d) Rejilla: la condición inicial no vive en `domain[0]`.** Un error frecuente es integrar
# con `t_span=(domain[0], domain[-1])`. Funciona mientras la rejilla empiece en 0, pero si pides
# `xdata` desde el día 20 la epidemia "arranca" el día 20 con `I0` infectados.

# %%
def sir_arranca_en_domain0(pars, domain):
    beta, gamma, S0, I0 = pars
    rhs = lambda t, y: [-beta * y[0] * y[1], beta * y[0] * y[1] - gamma * y[1], gamma * y[1]]
    return solve_ivp(rhs, (domain[0], domain[-1]), [S0, I0, 0.0], t_eval=domain,
                     rtol=1e-6, atol=1e-8).y      # MAL: t0 = domain[0]

x20 = np.linspace(20, 80, 121)
p = model_sir.nominal
print(f"I(20) correcto           : {sir_continuo(p, x20)[1, 0]:.1f}")
print(f"I(20) arrancando en x[0] : {sir_arranca_en_domain0(p, x20)[1, 0]:.1f}  (= I0)")

plt.figure(figsize=(8, 4))
plt.plot(model_sir.domain, sir_continuo(p, model_sir.domain)[1], "b-", linewidth=2, label="t0 = 0 (correcto)")
plt.plot(x20, sir_arranca_en_domain0(p, x20)[1], "r--", linewidth=1.5, label="t0 = domain[0] = 20 (mal)")
plt.xlabel("Tiempo (días)"); plt.ylabel("Infectados I(t)"); plt.grid(True); plt.legend()
plt.title("La condición inicial pertenece a t0, no al primer punto pedido")
plt.tight_layout(); plt.savefig("demo_00_rejilla.png")

# %% [markdown]
# **(e) Rejilla: olvidar `t_eval`.** Sin `t_eval`, `solve_ivp` devuelve sus pasos adaptativos
# (pocos y desiguales) y la salida ya no tiene `len(domain)` columnas. `verificar_modelo`
# (sección 8) lo detecta.

# %%
sol = solve_ivp(lambda t, y: [-4.5e-4 * y[0] * y[1], 4.5e-4 * y[0] * y[1] - 0.125 * y[1], 0.125 * y[1]],
                (0, 80), [949.5, 25.5, 0.0])
print("sin t_eval:", sol.y.shape, " frente a len(domain) =", len(model_sir.domain))

# %% [markdown]
# **(f) Unidades.** Rangos, condiciones iniciales y rejilla deben estar en las mismas unidades. En
# `sir_continuo`, `beta` va en 1/(persona·día) porque el término es `beta*S*I` con S e I en
# personas. Si escribes S0 e I0 como fracciones de la población sin cambiar `beta`, el número
# reproductivo $R_0 = \beta S_0/\gamma$ se hunde; si la rejilla está en semanas y `gamma` en 1/día, el
# tiempo de recuperación queda multiplicado por 7.

# %%
b, g = 4.5e-4, 0.125
print(f"R0 con S0 en personas (950) : {b * 950 / g:.2f}")
print(f"R0 con S0 en fracción (0.95): {b * 0.95 / g:.4f}   <- mismas ecuaciones, otra escala")

# %% [markdown]
# ## 8. Comprobar tu modelo antes de usarlo
#
# `verificar_modelo.py` (ayuda del curso, no del paquete) llama a tu modelo y revisa la forma de
# la salida, que no haya NaN y que el valor en un tiempo no dependa de la rejilla que pidas (lo
# que delata el error (d)). Úsalo cada vez que escribas un modelo nuevo.

# %%
grid = np.linspace(0, 80, 161)
verificar_modelo(logistica_discreta, [0.3, 500, 10], np.arange(0, 51), ["Poblacion"])
verificar_modelo(sir_continuo, model_sir.nominal, grid, ["S", "I", "R"], t0=0.0)
verificar_modelo(sir_transpuesto, model_sir.nominal, grid, ["S", "I", "R"])
verificar_modelo(sir_arranca_en_domain0, model_sir.nominal, grid, ["S", "I", "R"])

# %% [markdown]
# ## 9. Ejercicio: programa un SIR discreto desde cero
#
# Con un paso de un día, el SIR discreto es
#
# $$S_{t+1} = S_t - \beta S_t I_t,\quad I_{t+1} = I_t + \beta S_t I_t - \gamma I_t,\quad R_{t+1} = R_t + \gamma I_t.$$
#
# 1. Copia `sir_discreto_plantilla.py` como `mi_sir_discreto.py`.
# 2. Completa el `TODO` del bucle con las tres ecuaciones.
# 3. Pásalo por `verificar_modelo` con `output_names = ["S", "I", "R"]`.
# 4. Envuélvelo en `UserFunctionModel` con los mismos `range`, `names` y `output_names` que
#    `sir_continuo` y la rejilla `np.arange(0, 81)`.
# 5. Evalúa el nominal y grafica I(t) junto a la de `sir_continuo`.
#
# **Preguntas.** ¿Por qué no coinciden las dos curvas de I aunque las ecuaciones "son las
# mismas"? ¿Qué pasa con S si `beta*S0 > 1`? ¿Se conserva S + I + R en tu modelo?
#
# La plantilla (ya cumple el contrato; solo le falta la dinámica):

# %%
with open("sir_discreto_plantilla.py", encoding="utf-8") as f:
    print(f.read())

# %% [markdown]
# La solución está en `solucion/sir_discreto_solucion.py`. Ábrela cuando hayas terminado.

# %%
from sir_discreto_plantilla import sir_discreto_plantilla

verificar_modelo(sir_discreto_plantilla, model_sir.nominal, np.arange(0, 81), ["S", "I", "R"])
m_ej = UserFunctionModel(func=sir_discreto_plantilla, names=["beta", "gamma", "S0", "I0"],
                         range=rangos_sir, domain=np.arange(0, 81),
                         output_names=["S", "I", "R"], output="I")
I_ej = m_ej.evaluate(m_ej.nominal).ravel()
if np.allclose(I_ej, I_ej[0]):
    print("La plantilla cumple el contrato, pero I(t) sigue constante: completa el TODO.")
else:
    print(f"I(t) ya cambia: pico de {I_ej.max():.1f} personas en el día {int(np.argmax(I_ej))}.")

# %% [markdown]
# ## 10. Resumen: lista de chequeo
#
# * `func(pars, domain, **opt)` devuelve un arreglo **(n_salidas, len(domain))**, evaluado en los puntos de `domain`.
# * `pars[i]` = fila *i* de `range` = `names[i]`: parámetros primero, condiciones iniciales después.
# * La condición inicial vale en un `t0` fijo que el modelo conoce (argumento `t0` u `opt`), no en `domain[0]`.
# * Discreto: un bucle de `t0` al último punto pedido. Continuo: `solve_ivp(..., t_eval=domain)`.
# * `UserFunctionModel(func, names, range, domain, output_names, opt=...)` es la tabla `T` de MATLAB.
# * `model.evaluate(pars)` → (salidas, tiempo); `model.evaluate_batch(M)` con `M` de N x Np → (N, salidas, tiempo).
# * Una sola salida: `model.set_output("I")`. Un factor fijo: rango `[v, v]` o `model.fix("I0")`.
#
# En el módulo 1 estos mismos dos modelos (`discrete_logistic_userdefined.py` y
# `sir_continuous_userdefined.py`) se exploran con ensambles de trayectorias.
