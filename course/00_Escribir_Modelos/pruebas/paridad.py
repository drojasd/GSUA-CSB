"""Paridad MATLAB vs Python de los modelos del módulo 0.

Para cada modelo corre en Python los 20 vectores de pruebas/pars_<modelo>.csv sobre la rejilla
pruebas/x_<modelo>.csv y compara con pruebas/salidas_<modelo>.csv (generados por
generar_paridad.m).

Error relativo: |py - matlab| / max|fila de matlab| (escala por salida y por vector).
Tolerancias: 1e-8 en los discretos, 1e-4 en los continuos (ode45 frente a solve_ivp con las
mismas tolerancias).

Uso:  python pruebas/paridad.py     (también se puede correr con pytest)
"""
import importlib
import pathlib
import sys

import numpy as np

PRUEBAS = pathlib.Path(__file__).resolve().parent
MOD0 = PRUEBAS.parent
sys.path[:0] = [str(MOD0 / "python"), str(MOD0 / "python" / "solucion")]

MODELOS = {                     # nombre: tipo
    "logistica_discreta": "discreto",
    "sir_continuo": "continuo",
    "sir_discreto_plantilla": "discreto",
    "sir_discreto_solucion": "discreto",
}
TOL = {"discreto": 1e-8, "continuo": 1e-4}


def comparar(nombre):
    f = getattr(importlib.import_module(nombre), nombre)
    P = np.loadtxt(PRUEBAS / f"pars_{nombre}.csv", delimiter=",", ndmin=2)
    x = np.loadtxt(PRUEBAS / f"x_{nombre}.csv", delimiter=",", ndmin=1)
    M = np.loadtxt(PRUEBAS / f"salidas_{nombre}.csv", delimiter=",", ndmin=2)
    err = 0.0
    for p, m in zip(P, M):
        y = np.atleast_2d(f(p, x))                   # (salidas, tiempo)
        ym = m.reshape(y.shape)
        escala = np.maximum(np.abs(ym).max(axis=1, keepdims=True), 1e-12)
        err = max(err, float(np.max(np.abs(y - ym) / escala)))
    return err, len(P)


def test_paridad():
    fallos = []
    for nombre, tipo in MODELOS.items():
        err, n = comparar(nombre)
        tol = TOL[tipo]
        estado = "OK   " if err < tol else "FALLA"
        print(f"  {estado} {nombre:<24s} {n} vectores, error relativo máximo {err:.2e} (tol {tol:g})")
        if err >= tol:
            fallos.append(nombre)
    assert not fallos, f"paridad fallida en: {fallos}"


if __name__ == "__main__":
    test_paridad()
    print("PARIDAD OK")
