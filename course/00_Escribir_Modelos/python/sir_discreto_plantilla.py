"""
sir_discreto_plantilla.py — ejercicio del módulo 0: un SIR en tiempo discreto.

Copia este archivo como mi_sir_discreto.py y completa el TODO.

Con un paso de un día, el SIR discreto es:
    S(t+1) = S(t) - beta*S(t)*I(t)
    I(t+1) = I(t) + beta*S(t)*I(t) - gamma*I(t)
    R(t+1) = R(t) + gamma*I(t)

Factores (primero parámetros, después condiciones iniciales):
    pars[0] -> beta  (1/(persona*día))
    pars[1] -> gamma (1/día)
    pars[2] -> S0    (personas)
    pars[3] -> I0    (personas)
R0 = 0 al inicio, en el día t0.

Salida: arreglo (3, len(domain)) = [S; I; R] en los días de domain.
"""
import numpy as np


def sir_discreto_plantilla(pars, domain, t0=0):
    beta, gamma, S0, I0 = pars
    domain = np.asarray(domain, dtype=float)

    dias = np.arange(t0, int(np.ceil(domain.max())) + 1)   # t0, t0+1, ..., último día pedido
    Y = np.zeros((3, len(dias)))                            # 3 salidas -> 3 filas
    Y[:, 0] = [S0, I0, 0.0]

    for k in range(len(dias) - 1):
        S, I, R = Y[:, k]
        # TODO: escribe aquí las tres ecuaciones de actualización.
        #       Mientras tanto la plantilla deja S, I y R constantes.
        Y[:, k + 1] = [S, I, R]

    # Devolver exactamente los puntos pedidos (interpolación lineal entre días)
    return np.vstack([np.interp(domain, dias, fila) for fila in Y])
