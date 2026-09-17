"""
DEMO 05: Estimación de Parámetros e Identificabilidad Práctica (Python / gsua_csb)
Curso guiado de GSUA-CSB

Cierra el ciclo de calibración y llega a la lección central del curso:
UN BUEN AJUSTE NO GARANTIZA PARÁMETROS CONFIABLES.
  - ESTIMACIÓN (parameter_estimation): ajuste por optimización multicomienzo.
  - IDENTIFICABILIDAD (identifiability_analysis): generamos un ENSAMBLE re-ajustando
    varias realizaciones de ruido (bootstrap); su dispersión y su matriz de
    correlación revelan qué parámetros están entrelazados.
  - CANTIDAD IDENTIFICABLE: aunque beta, gamma y S0 se dispersen, su combinación
    R0 = beta*S0/gamma queda determinada.

Modelo: SIR continuo. Observamos SOLO Infectados; I0 fijo; estimamos beta, gamma, S0.
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from gsua_csb import (
    parameter_estimation,
    identifiability_analysis,
    plot_identifiability_correlation,
    plot_identifiability_index,
)

module1_py_path = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "01_Introduccion_Modelos", "python")
)
if module1_py_path not in sys.path:
    sys.path.append(module1_py_path)
from sir_continuous_userdefined import create_sir_continuous_model


def run_demo():
    print("=========================================================================")
    print("  MÓDULO 05: ESTIMACIÓN DE PARÁMETROS E IDENTIFICABILIDAD (PYTHON)")
    print("=========================================================================\n")

    # 1. Modelo SIR con I0 FIJO (range degenerado); estimamos beta, gamma, S0.
    model = create_sir_continuous_model(I0_range=(5.0, 5.0), output_selection="Infectados")
    free_idx = [i for i, fx in enumerate(model.fixed) if not fx]
    free_names = [model.names[i] for i in free_idx]
    print(f"Factores a estimar: {free_names} (I0 fijo en 5); se observa solo Infectados.\n")

    # 2. Datos "observados" sintéticos (verdad + 5% de ruido).
    true_par = np.array([5e-4, 0.12, 950.0, 5.0])          # beta, gamma, S0, I0 (fijo)
    true_free = true_par[free_idx]
    y_clean = np.squeeze(model.evaluate(true_par, model.domain))
    rng = np.random.default_rng(3)
    ydata = y_clean * (1 + 0.05 * rng.standard_normal(y_clean.shape))

    # 3. Estimación multicomienzo sobre los datos observados.
    pe = parameter_estimation(model, model.domain, ydata, n=10, solver="least_squares", seed=0)
    best = pe.x[0]
    print("--- Mejor ajuste ---")
    print(f"   beta  = {best[0]:.3e}   (verdad 5.000e-04)")
    print(f"   gamma = {best[1]:.4f}    (verdad 0.1200)")
    print(f"   S0    = {best[2]:.1f}     (verdad 950.0)")
    print(f"   El ajuste es EXCELENTE (costo = {pe.cost[0]:.4g}). ¿Son confiables los parámetros?\n")

    y_best = np.squeeze(model.evaluate(best, model.domain))
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(model.domain, ydata, ".", color="tab:red", ms=5, label="Datos observados")
    ax.plot(model.domain, y_best, "-", color="tab:blue", lw=2, label="Mejor ajuste")
    ax.set_xlabel("Tiempo (Días)"); ax.set_ylabel("Infectados I(t)")
    ax.set_title("Estimación de parámetros: el modelo ajusta muy bien los datos")
    ax.grid(True); ax.legend()
    fig.tight_layout(); fig.savefig("demo_05_mejor_ajuste.png")

    # 4. Identificabilidad por bootstrap: ensamble de estimaciones sobre B realizaciones
    #    de ruido. La dispersión = incertidumbre; su estructura = entrelazamientos.
    B = 20
    print(f"Generando ensamble bootstrap ({B} re-ajustes)...")
    est = np.zeros((B, model.range.shape[0]))
    for b in range(B):
        yb = y_clean * (1 + 0.05 * rng.standard_normal(y_clean.shape))
        pe_b = parameter_estimation(model, model.domain, yb, n=3, solver="least_squares", seed=b)
        est[b] = pe_b.x[0]

    ia = identifiability_analysis(model, est, correction=False)
    est_free = est[:, free_idx]
    cv = lambda v: np.std(v) / np.mean(v) * 100
    print("\n--- Dispersión del ensamble (coeficiente de variación) ---")
    for j, name in enumerate(free_names):
        print(f"   {name:6s} CV = {cv(est_free[:, j]):5.1f}%   índice identif = {np.asarray(ia.index)[free_idx][j]:.2f}")

    Rcorr = np.corrcoef(est_free, rowvar=False)
    print("\n--- Correlación entre parámetros ---")
    print("        " + "  ".join(f"{n:>6s}" for n in free_names))
    for i, n in enumerate(free_names):
        print(f"   {n:6s} " + "  ".join(f"{Rcorr[i, j]:6.2f}" for j in range(len(free_names))))
    Rabs = np.abs(Rcorr.copy()); np.fill_diagonal(Rabs, 0.0)
    a, b2 = np.unravel_index(int(np.argmax(Rabs)), Rabs.shape)
    print(f"\nPar más entrelazado: {free_names[a]} <-> {free_names[b2]}  (|corr| = {Rabs[a, b2]:.2f})")
    if Rabs[a, b2] > 0.9:
        print("=> |corr| > 0.9: aunque el ajuste sea perfecto, los datos identifican su")
        print("   COMBINACIÓN, no cada parámetro por separado (no-identificabilidad práctica).\n")

    # Figuras diagnósticas de identificabilidad.
    fig2, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    plot_identifiability_correlation(ia, ax=axes[0]); axes[0].set_title("Correlación entre parámetros")
    plot_identifiability_index(ia, ax=axes[1]); axes[1].set_title("Índice de identificabilidad")
    fig2.tight_layout(); fig2.savefig("demo_05_identificabilidad.png")

    # 5. La cantidad identificable: R0 = beta*S0/gamma.
    R0_ens = est[:, 0] * est[:, 2] / est[:, 1]
    R0_true = true_par[0] * true_par[2] / true_par[1]
    print("--- Número reproductivo básico R0 = beta*S0/gamma ---")
    cv_ind = np.mean([cv(est_free[:, j]) for j in range(len(free_names))])
    print(f"   CV de los parámetros individuales ~ {cv_ind:.1f}%, pero CV de R0 = {cv(R0_ens):.2f}%")
    print(f"   R0 estimado = {R0_ens.mean():.2f} ± {R0_ens.std():.2f}   (verdad = {R0_true:.2f})")
    print("=> Moraleja: la COMBINACIÓN identificable (R0) es confiable aunque sus")
    print("   componentes por separado no lo sean.\n")
    print("Figuras guardadas: demo_05_mejor_ajuste.png, demo_05_identificabilidad.png")

    # 6. PREGUNTAS Y EJERCICIOS PARA EL ESTUDIANTE
    print("\n=========================================================================")
    print("  PREGUNTAS DE ANÁLISIS Y EJERCICIOS DOCENTES (MÓDULO 05)  ")
    print("=========================================================================")
    print("1. Explica la correlación beta<->S0 usando R0 = beta*S0/gamma.")
    print("2. Fija S0 en su valor real (S0_range=(950,950)) y repite. ¿Mejora la identificabilidad?")
    print("3. Observa los 3 compartimentos (output_selection=None) en vez de solo Infectados.")
    print("   ¿Se rompen los entrelazamientos? ¿Bajan los CV individuales?")
    print("4. Sube el ruido de 5% a 15%. ¿Cómo cambian los CV y el CV de R0?")
    print("=========================================================================\n")


if __name__ == "__main__":
    run_demo()
