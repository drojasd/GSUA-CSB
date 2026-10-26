"""
sir_arranca_en_domain0.py — el SIR continuo integrando desde domain[0].

Solo existe para el módulo 0: funciona mientras la rejilla empiece en 0, pero si se pide
xdata desde el día 20 la epidemia "arranca" el día 20 con I0 infectados. La condición
inicial pertenece a un t0 fijo, no al primer punto pedido.
"""
import numpy as np
from scipy.integrate import solve_ivp


def sir_arranca_en_domain0(pars, domain):
    beta, gamma, S0, I0 = pars
    domain = np.asarray(domain, dtype=float)
    rhs = lambda t, y: [-beta * y[0] * y[1], beta * y[0] * y[1] - gamma * y[1], gamma * y[1]]
    return solve_ivp(rhs, (domain[0], domain[-1]), [S0, I0, 0.0], t_eval=domain,
                     rtol=1e-6, atol=1e-8).y      # MAL: t0 = domain[0]
