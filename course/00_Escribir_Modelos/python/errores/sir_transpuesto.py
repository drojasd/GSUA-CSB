"""
sir_transpuesto.py — el SIR continuo devolviendo la salida transpuesta (tiempo, salidas).

Solo existe para el módulo 0: en Python nadie corrige la orientación, así que
set_output("I") selecciona una fila que es un instante de tiempo, no la salida I.
"""
from sir_continuo import sir_continuo


def sir_transpuesto(pars, domain, t0=0.0):
    return sir_continuo(pars, domain, t0).T          # (len(domain), 3): MAL
