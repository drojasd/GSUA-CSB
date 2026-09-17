"""
DEMO 03: Análisis de Incertidumbre y Filtrado Monte Carlo (Python / gsua_csb)
Curso guiado de GSUA-CSB

Propaga la incertidumbre de los factores hacia la salida del modelo y evalúa:
  - el ENSAMBLE de trayectorias (uncertainty_analysis),
  - la BANDA de percentiles P5/P50/P95 y su métrica de cobertura (coverage_metric),
  - el FILTRADO Monte Carlo que separa combinaciones aceptables/no aceptables
    (monte_carlo_filter).

Modelo de trabajo: SIR continuo del Módulo 01, observando la curva de Infectados.
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from gsua_csb import (
    design_matrix,
    uncertainty_analysis,
    coverage_metric,
    monte_carlo_filter,
    plot_uncertainty,
    plot_mcf,
)

# Reutilizamos el modelo SIR del Módulo 01 (ruta relativa a este script).
module1_py_path = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "01_Introduccion_Modelos", "python")
)
if module1_py_path not in sys.path:
    sys.path.append(module1_py_path)
from sir_continuous_userdefined import create_sir_continuous_model


def run_demo():
    print("=========================================================================")
    print("  MÓDULO 03: ANÁLISIS DE INCERTIDUMBRE Y FILTRADO MONTE CARLO (PYTHON)")
    print("=========================================================================\n")

    # 1. Modelo SIR con la salida de Infectados seleccionada
    #    (equivalente a fijar T.Properties.CustomProperties.output = 2 en MATLAB).
    #    Rangos "previos" acotados (± ~20% del valor real): a propósito más estrechos
    #    que en el Módulo 02, para propagar una incertidumbre plausible.
    model = create_sir_continuous_model(
        beta_range=(4.0e-4, 6.0e-4),
        gamma_range=(0.10, 0.15),
        S0_range=(920.0, 980.0),
        I0_range=(3.0, 8.0),
        output_selection="Infectados",
    )

    # 2. Datos "observados" sintéticos: verdad conocida + ruido de medición del 8%.
    true_par = np.array([0.0005, 0.12, 950.0, 5.0])          # beta, gamma, S0, I0
    y_true = np.squeeze(model.evaluate(true_par, model.domain))
    rng = np.random.default_rng(7)
    ydata = y_true * (1 + 0.08 * rng.standard_normal(y_true.shape))

    # 3. Muestreo del espacio de incertidumbre (Hipercubo Latino).
    M = design_matrix(model, n=300, method="latin_hypercube", seed=0)

    # 4. Análisis de incertidumbre: ensamble de trayectorias sobre M.
    #    A diferencia de MATLAB, en Python NADA se grafica implícitamente:
    #    uncertainty_analysis solo calcula; las figuras las pedimos con plot_*.
    ua = uncertainty_analysis(model, M, y_exp=ydata)          # ua.Y: (N, Nd); ua.y_nom: referencia
    print(f"Ensamble simulado: {ua.Y.shape[0]} trayectorias de {ua.Y.shape[1]} puntos.")

    # 5. Banda de percentiles + métricas de cobertura.
    #    coverage_metric devuelve (cost_data, cost_band, P5, P50, P95); < 1 = dentro de tolerancia.
    cost_data, cost_band, P5, P50, P95 = coverage_metric(ua.Y, ydata, margin=0.1)
    print("\n--- Métricas de cobertura (margen 10%) ---")
    print(f"cost_data = {cost_data:.3f}  (exactitud: mediana vs datos)"
          + ("  -> datos ALCANZABLES" if cost_data < 1 else "  -> datos NO alcanzables"))
    print(f"cost_band = {cost_band:.3g}  (precisión: ancho de la banda P5-P95)"
          + ("  -> banda angosta" if cost_band < 1
             else "  -> banda ancha (normal >>1 sin ajustar; el Módulo 05 la reduce)"))

    # P5/P50/P95 vienen como (nout, Nd); con una sola salida tomamos la fila 0.
    p5, p50, p95 = np.atleast_2d(P5)[0], np.atleast_2d(P50)[0], np.atleast_2d(P95)[0]

    # 6. Filtrado Monte Carlo: subconjuntos low/high por parámetro.
    mcf = monte_carlo_filter(model, M, ua.Y, ua.y_nom)
    print(f"\nFiltrado Monte Carlo: {mcf.low.shape[0]} corridas 'low' y "
          f"{mcf.high.shape[0]} 'high' sobre {len(mcf.names)} factores libres {list(mcf.names)}.")

    # 7. Figuras: (a) banda de incertidumbre propia, (b) plot_uncertainty de la toolbox,
    #    (c) paneles ECDF del filtrado Monte Carlo.
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.fill_between(model.domain, p5, p95, color=(0.80, 0.88, 0.98), label="Banda P5-P95")
    ax.plot(model.domain, p50, color="tab:blue", lw=2, label="Mediana P50")
    ax.plot(model.domain, ydata, ".", color="tab:red", ms=4, label="Datos observados")
    ax.set_xlabel("Tiempo (Días)"); ax.set_ylabel("Infectados I(t)")
    ax.set_title(f"Incertidumbre propagada (cost_data={cost_data:.2f}, cost_band={cost_band:.2f})")
    ax.grid(True); ax.legend()
    fig.tight_layout(); fig.savefig("demo_03_banda_incertidumbre.png")

    fig2, ax2 = plt.subplots(figsize=(9, 4.5))
    plot_uncertainty(ua, ax=ax2)
    fig2.tight_layout(); fig2.savefig("demo_03_uncertainty_toolbox.png")

    # plot_mcf necesita la matriz de diseño restringida a columnas libres; como aquí
    # los 4 factores son libres, pasamos M completa.
    plot_mcf(mcf, M)                    # abre su propia figura (un panel por parámetro)
    plt.savefig("demo_03_mcf.png")
    print("\nFiguras guardadas: demo_03_banda_incertidumbre.png, demo_03_uncertainty_toolbox.png, demo_03_mcf.png")

    # 8. PREGUNTAS Y EJERCICIOS PARA EL ESTUDIANTE
    print("\n=========================================================================")
    print("  PREGUNTAS DE ANÁLISIS Y EJERCICIOS DOCENTES (MÓDULO 03)  ")
    print("=========================================================================")
    print("1. Sube el ruido de 8% a 25%. ¿Cómo cambian cost_data y cost_band?")
    print("2. Reduce el rango de beta a la mitad. ¿La banda P5-P95 se vuelve más angosta?")
    print("3. En plot_mcf, ¿qué parámetro separa más sus curvas low/high? ¿Por qué?")
    print("4. ¿Puede cost_band ser bajo (banda angosta) y cost_data alto a la vez? ¿Qué significaría?")
    print("=========================================================================\n")


if __name__ == "__main__":
    run_demo()
