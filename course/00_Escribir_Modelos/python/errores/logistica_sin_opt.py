"""
logistica_sin_opt.py — la logística discreta sin argumento para opt.

Solo existe para el módulo 0: si el modelo se declara con opt={"t0": 0}, UserFunctionModel
llama logistica_sin_opt(pars, domain, t0=0) y Python falla con "unexpected keyword argument".
La versión correcta es logistica_discreta (declara t0=0, o **opt).
"""
import numpy as np


def logistica_sin_opt(pars, domain):
    r, K, N0 = pars
    domain = np.asarray(domain, dtype=float)
    t = np.arange(0, int(np.ceil(domain.max())) + 1)
    N = np.zeros(len(t))
    N[0] = N0
    for k in range(len(t) - 1):
        N[k + 1] = N[k] + r * N[k] * (1 - N[k] / K)
    return np.interp(domain, t, N)
