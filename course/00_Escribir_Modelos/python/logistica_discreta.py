"""
logistica_discreta.py — crecimiento logístico en tiempo discreto (modelo de usuario).

    y = logistica_discreta(pars, domain, t0=0)

Ecuación en diferencias:  N(t+1) = N(t) + r*N(t)*(1 - N(t)/K)

Factores, en este orden (primero parámetros, después condición inicial):
    pars[0] -> r  : tasa de crecimiento por paso (adimensional)
    pars[1] -> K  : capacidad de carga (individuos)
    pars[2] -> N0 : población en t0 (individuos)

domain = la rejilla de tiempos donde se quiere la salida (model.domain o el xdata de evaluate).
t0     = el paso donde vale N0. Llega por opt={"t0": ...}; por defecto 0.

Salida: arreglo (len(domain),) con N en los puntos de domain: una sola salida.
"""
import numpy as np


def logistica_discreta(pars, domain, t0=0):
    r, K, N0 = pars
    domain = np.asarray(domain, dtype=float)

    t = np.arange(t0, int(np.ceil(domain.max())) + 1)   # un punto por paso de tiempo
    N = np.zeros(len(t))                                # 1 salida -> 1 fila
    N[0] = N0
    for k in range(len(t) - 1):
        N[k + 1] = N[k] + r * N[k] * (1 - N[k] / K)

    return np.interp(domain, t, N)                      # valores en los puntos pedidos
