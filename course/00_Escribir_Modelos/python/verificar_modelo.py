"""
verificar_modelo.py — revisa que un modelo de usuario cumpla el contrato de gsua_csb.

Es una ayuda del curso, no una función del paquete. Llama al modelo a mano, antes de
envolverlo en UserFunctionModel, y comprueba lo que el paquete da por hecho:

1. devuelve un arreglo de forma (n_salidas, len(domain)) — o (len(domain),) si hay una salida;
2. no hay NaN ni Inf;
3. el valor en un tiempo no depende de qué rejilla le pidas (la condición inicial vive en t0,
   no en el primer punto de domain).
"""
import numpy as np


def verificar_modelo(func, pars, domain, output_names, **opt):
    nombre = getattr(func, "__name__", str(func))
    print(f"Revisando {nombre} ...")
    domain = np.asarray(domain, dtype=float)
    n_out, n_t = len(output_names), len(domain)
    ok = True

    try:
        y = np.asarray(func(np.asarray(pars, dtype=float), domain, **opt), dtype=float)
    except TypeError as e:
        print(f"  FALLA  no se pudo llamar func(pars, domain, **opt): {e}")
        print("  => revisa los puntos marcados\n")
        return False
    print(f"  ok     se llama como {nombre}(pars, domain{', **opt' if opt else ''})")

    print(f"  info   salida {y.shape};  len(domain) = {n_t};  output_names: {n_out}")
    if y.shape == (n_out, n_t) or (n_out == 1 and y.shape == (n_t,)):
        print("  ok     forma (salidas, tiempo)")
    elif y.shape == (n_t, n_out):
        print("  FALLA  forma (tiempo, salidas): devuelve la transpuesta. En Python nadie la corrige")
        ok = False
    else:
        print(f"  FALLA  debería ser ({n_out}, {n_t}): una fila por salida, una columna por punto de domain")
        ok = False

    if not np.all(np.isfinite(y)):
        print("  FALLA  hay NaN o Inf")
        ok = False
    else:
        print("  ok     sin NaN ni Inf")

    if ok and n_t >= 4:
        mitad = domain[n_t // 2:]
        y_mitad = np.asarray(func(np.asarray(pars, dtype=float), mitad, **opt), dtype=float)
        y_cola = y[..., n_t // 2:]
        if np.allclose(y_mitad, y_cola, rtol=1e-4, atol=1e-6):
            print("  ok     pedir solo la segunda mitad de la rejilla da los mismos valores")
        else:
            print("  FALLA  con la rejilla recortada cambian los valores: ¿arranca en domain[0] en vez de t0?")
            ok = False

    print("  => cumple el contrato\n" if ok else "  => revisa los puntos marcados\n")
    return ok
