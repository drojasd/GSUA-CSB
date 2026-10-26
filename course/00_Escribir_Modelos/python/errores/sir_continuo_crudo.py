"""
sir_continuo_crudo.py — el mismo SIR que sir_continuo, pero con los pasos del integrador.

Solo existe para comparar en el módulo 0: solve_ivp sin t_eval devuelve sus pasos adaptativos
(pocos y desiguales, con tolerancias por defecto) y aquí se interpola linealmente entre ellos,
igual que hace gsua_eval en MATLAB con la estructura cruda de ode45.
"""
import numpy as np
from scipy.integrate import solve_ivp


def pasos_del_integrador(pars, tf, t0=0.0):
    """Los tiempos y estados que solve_ivp calcula por su cuenta (sin t_eval)."""
    beta, gamma, S0, I0 = pars
    rhs = lambda t, y: [-beta * y[0] * y[1], beta * y[0] * y[1] - gamma * y[1], gamma * y[1]]
    sol = solve_ivp(rhs, (t0, tf), [S0, I0, 0.0])                     # tolerancias por defecto
    return sol.t, sol.y


def sir_continuo_crudo(pars, domain, t0=0.0):
    domain = np.asarray(domain, dtype=float)
    t, y = pasos_del_integrador(pars, domain.max(), t0)
    return np.vstack([np.interp(domain, t, fila) for fila in y])   # interpolación lineal
