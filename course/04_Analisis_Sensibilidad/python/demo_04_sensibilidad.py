"""
DEMO 04: Análisis de Sensibilidad Global (Python / gsua_csb)
Curso guiado de GSUA-CSB

Cuantifica cuánto contribuye cada factor a la varianza de la salida:
  - Índices de Sobol de primer orden  S_i,
  - Índices de orden total            ST_i,
  - Índices dependientes del tiempo   S_i(t), ST_i(t),
  - Comparación entre métodos estimadores (xiao vs saltelli).

Con esto decidimos qué parámetros ESTIMAR (los influyentes) y cuáles FIJAR.
Modelo: SIR continuo del Módulo 01, salida Infectados.
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from gsua_csb import (
    design_matrix,
    sensitivity_analysis,
    plot_sensitivity_bar,
    plot_sensitivity_area,
)

module1_py_path = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "01_Introduccion_Modelos", "python")
)
if module1_py_path not in sys.path:
    sys.path.append(module1_py_path)
from sir_continuous_userdefined import create_sir_continuous_model


def run_demo():
    print("=========================================================================")
    print("  MÓDULO 04: ANÁLISIS DE SENSIBILIDAD GLOBAL (PYTHON)")
    print("=========================================================================\n")

    # 1. Modelo SIR sobre TODO el espacio de factores (rangos amplios por defecto).
    model = create_sir_continuous_model(output_selection="Infectados")
    free_names = [n for n, fx in zip(model.names, model.fixed) if not fx]

    # 2. Muestreo base para el estimador de índices.
    M = design_matrix(model, n=500, method="latin_hypercube", seed=0)

    # 3. Análisis de sensibilidad con el método Xiao (estimador propio de la toolbox).
    res = sensitivity_analysis(model, M, method="xiao")
    print("Índices de sensibilidad (método Xiao):")
    for name, si, sti in zip(free_names, res.Si, res.STi):
        print(f"   {name:6s}  S_i = {si: .3f}   ST_i = {sti: .3f}")

    i_max = int(np.argmax(res.STi))
    i_min = int(np.argmin(res.STi))
    print(f"\nFactor MÁS influyente : {free_names[i_max]} (ST_i = {res.STi[i_max]:.3f})")
    print(f"Factor MENOS influyente: {free_names[i_min]} (ST_i = {res.STi[i_min]:.3f})")
    print("=> Estima los factores con ST_i alto; considera FIJAR los de ST_i ~ 0.\n")

    # 4. Gráfica de barras S_i y ST_i (dos paneles con las funciones de la toolbox).
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    plot_sensitivity_bar(res, index="Si", ax=axes[0]); axes[0].set_title("Primer orden S_i")
    plot_sensitivity_bar(res, index="STi", ax=axes[1]); axes[1].set_title("Orden total ST_i")
    fig.tight_layout(); fig.savefig("demo_04_indices_barras.png")

    # 5. Sensibilidad dependiente del tiempo (¿cuándo importa cada factor?).
    fig2, ax2 = plt.subplots(figsize=(9, 4.5))
    plot_sensitivity_area(res, xdata=model.domain, index="STi_vec", ax=ax2)
    ax2.set_title("Fracción de ST_i(t): factor dominante en cada instante")
    fig2.tight_layout(); fig2.savefig("demo_04_sensibilidad_tiempo.png")

    # 6. Comparación entre métodos estimadores (robustez del resultado).
    res_salt = sensitivity_analysis(model, M, method="saltelli")
    print("Comparación ST_i (xiao vs saltelli):")
    for name, a, b in zip(free_names, res.STi, res_salt.STi):
        print(f"   {name:6s}  xiao = {a: .3f}   saltelli = {b: .3f}")

    print("\nFiguras guardadas: demo_04_indices_barras.png, demo_04_sensibilidad_tiempo.png")

    # 7. PREGUNTAS Y EJERCICIOS PARA EL ESTUDIANTE
    print("\n=========================================================================")
    print("  PREGUNTAS DE ANÁLISIS Y EJERCICIOS DOCENTES (MÓDULO 04)  ")
    print("=========================================================================")
    print("1. ¿Qué factor domina la varianza de Infectados? ¿Coincide con tu intuición?")
    print("2. ¿Para algún factor ST_i >> S_i? ¿Con quién interactúa? (pista: beta, S0, gamma via R0).")
    print("3. Cambia n de 500 a 2000. ¿Se estabilizan los índices? (convergencia del estimador).")
    print("4. Si fijas el factor menos sensible, ¿cambia mucho la banda de incertidumbre del Módulo 03?")
    print("=========================================================================\n")


if __name__ == "__main__":
    run_demo()
