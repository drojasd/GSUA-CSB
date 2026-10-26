"""
sir_continuo.py — modelo epidémico SIR en tiempo continuo (modelo de usuario).

    y = sir_continuo(pars, domain, t0=0.0)

dS/dt = -beta*S*I
dI/dt =  beta*S*I - gamma*I
dR/dt =  gamma*I

Factores, en este orden (primero parámetros, después condiciones iniciales):
    pars[0] -> beta  : tasa de transmisión (1/(persona*día))
    pars[1] -> gamma : tasa de recuperación (1/día)
    pars[2] -> S0    : susceptibles en t0 (personas)
    pars[3] -> I0    : infectados en t0 (personas)
R(t0) = 0 no es un factor: es una constante del modelo.

domain = la rejilla de tiempos (días) donde se quiere la salida.
t0     = el día donde valen S0 e I0. Llega por opt={"t0": ...}; por defecto 0.

Salida: arreglo (3, len(domain)), filas = salidas [S; I; R].
"""
import numpy as np
from scipy.integrate import solve_ivp


def sir_continuo(pars, domain, t0=0.0):
    beta, gamma, S0, I0 = pars
    domain = np.asarray(domain, dtype=float)

    def rhs(t, y):
        S, I, R = y
        return [-beta * S * I, beta * S * I - gamma * I, gamma * I]

    # solve_ivp avanza con pasos adaptativos y desiguales. t_eval le pide la solución
    # exactamente en los puntos de domain (con su interpolante de alto orden).
    sol = solve_ivp(rhs, (t0, domain.max()), [S0, I0, 0.0],
                    t_eval=domain, rtol=1e-6, atol=1e-8)
    return sol.y                                        # (3, len(domain))
