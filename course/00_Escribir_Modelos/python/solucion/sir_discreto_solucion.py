"""
sir_discreto_solucion.py — solución del ejercicio del módulo 0 (SIR discreto).

Ábrela después de intentar sir_discreto_plantilla.py por tu cuenta.
"""
import numpy as np


def sir_discreto_solucion(pars, domain, t0=0):
    beta, gamma, S0, I0 = pars
    domain = np.asarray(domain, dtype=float)

    dias = np.arange(t0, int(np.ceil(domain.max())) + 1)
    Y = np.zeros((3, len(dias)))
    Y[:, 0] = [S0, I0, 0.0]

    for k in range(len(dias) - 1):
        S, I, R = Y[:, k]
        nuevos = beta * S * I          # contagios del día
        recuperados = gamma * I        # recuperaciones del día
        Y[:, k + 1] = [S - nuevos, I + nuevos - recuperados, R + recuperados]

    return np.vstack([np.interp(domain, dias, fila) for fila in Y])
