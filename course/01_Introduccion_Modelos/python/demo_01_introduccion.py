"""
DEMO 01: Introducción Avanzada a Modelos Definidos por el Usuario en gsua_csb
Curso guiado de GSUA-CSB

Script ejecutable demostrativo para evaluar modelos en tiempo discreto y continuo,
mostrando ensamble de trayectorias, plano de fase epidemiológico y cálculo de R0.
"""
import numpy as np
import matplotlib.pyplot as plt
from gsua_csb import design_matrix
from discrete_logistic_userdefined import create_discrete_logistic_model
from sir_continuous_userdefined import create_sir_continuous_model

# ## Los modelos de este módulo
# Este módulo usa dos modelos de usuario ya escritos. Así están, completos (el módulo 0
# explica paso a paso cómo se escribe uno: la firma `func(pars, domain)`, el orden de los
# factores y la forma (salidas, tiempo) de la salida). Lo que ve gsua_csb son las funciones
# `discrete_logistic_func` y `sir_continuous_func`; `create_*_model` las envuelve en
# `UserFunctionModel`.
#
# Nota: ambos modelos arrancan la simulación en `domain[0]`. Es correcto en este curso, donde la
# rejilla siempre empieza en 0; si pides una rejilla que empiece después, la condición inicial
# se desplazaría con ella. El módulo 0 muestra cómo fijar `t0` para evitarlo.

import inspect
import discrete_logistic_userdefined as _dl
import sir_continuous_userdefined as _sir
for _f in (_dl.discrete_logistic_func, _dl.create_discrete_logistic_model,
           _sir._sir_rhs, _sir.sir_continuous_func, _sir.create_sir_continuous_model):
    print(inspect.getsource(_f))

def run_demo():
    print("=========================================================================")
    print("  MÓDULO 01: INTRODUCCIÓN A MODELOS EN TIEMPO DISCRETO Y CONTINUO (PYTHON)")
    print("=========================================================================\n")
    
    # --- 1. MODELO LOGÍSTICO DISCRETO ---
    print("--- 1. MODELO LOGÍSTICO DISCRETO ---")
    model_discrete = create_discrete_logistic_model()
    
    print(f"Factores (Parámetros + C.I.): {model_discrete.names}")
    print(f"Límites de los factores:\n{model_discrete.range}")
    print(f"Nombres de salidas: {model_discrete.output_names}")
    print(f"Punto nominal: {model_discrete.nominal}\n")
    
    # A. Simulación Nominal
    y_discrete_nom = np.squeeze(model_discrete.evaluate(model_discrete.nominal, model_discrete.domain))
    
    # B. Muestreo de Ensamble de Incertidumbre (design_matrix)
    M_discrete = design_matrix(model_discrete, n=25, method="latin_hypercube", seed=42)
    Y_discrete_ens = model_discrete.evaluate_batch(M_discrete, model_discrete.domain) # shape (25, len(domain))
    
    plt.figure(figsize=(9, 4.5))
    for i in range(25):
        plt.plot(model_discrete.domain, Y_discrete_ens[i], color='lightblue', alpha=0.6, linewidth=0.8)
    plt.plot(model_discrete.domain, y_discrete_nom, color='tab:blue', linewidth=2.5, label='Nominal')
    plt.xlabel('Paso de Tiempo (t)')
    plt.ylabel('Población N(t)')
    plt.title('Crecimiento Logístico Discreto: Ensamble de 25 Trayectorias')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("demo_logistic_discrete_ensemble.png")
    print("Gráfica del ensamble discreto guardada como 'demo_logistic_discrete_ensemble.png'")
    
    # --- 2. MODELO SIR CONTINUO ---
    print("\n--- 2. MODELO EPIDÉMICO SIR CONTINUO ---")
    model_sir = create_sir_continuous_model(output_selection=None) # Todas las salidas (3)
    
    print(f"Factores (Parámetros + C.I.): {model_sir.names}")
    print(f"Límites de los factores:\n{model_sir.range}")
    print(f"Nombres de salidas: {model_sir.output_names}")
    print(f"Punto nominal: {model_sir.nominal}\n")
    
    # A. Simulación Nominal de todos los compartimentos (S, I, R)
    y_sir_all = model_sir.evaluate(model_sir.nominal, model_sir.domain) # shape (3, len(domain))
    S_nom, I_nom, R_nom = y_sir_all[0], y_sir_all[1], y_sir_all[2]
    
    beta_nom, gamma_nom, S0_nom, I0_nom = model_sir.nominal
    R0_basic = (beta_nom * S0_nom) / gamma_nom
    
    print(f"Parámetros Nominales: beta = {beta_nom:.5f}, gamma = {gamma_nom:.2f}, S0 = {S0_nom:.0f}")
    print(f"Número Reproductivo Básico Calculado: R0 = {R0_basic:.2f}")
    if R0_basic > 1:
        print("-> R0 > 1: Existe brote epidémico (El número de infectados crece inicialmente).\n")
    else:
        print("-> R0 <= 1: No hay brote epidémico.\n")
        
    # B. Visualización Dinámica y Plano de Fase (S vs I)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    
    axes[0].plot(model_sir.domain, S_nom, color='tab:blue', linewidth=2, label='S(t) Susceptibles')
    axes[0].plot(model_sir.domain, I_nom, color='tab:red', linewidth=2, label='I(t) Infectados')
    axes[0].plot(model_sir.domain, R_nom, color='tab:green', linewidth=2, label='R(t) Recuperados')
    axes[0].set_xlabel('Tiempo (Días)')
    axes[0].set_ylabel('Población (Personas)')
    axes[0].set_title(f'Evolución Temporal SIR (R0 = {R0_basic:.2f})')
    axes[0].grid(True)
    axes[0].legend()
    
    axes[1].plot(S_nom, I_nom, color='purple', linewidth=2.5)
    axes[1].plot(S_nom[0], I_nom[0], 'ko', markersize=8, label=f'Inicio ({S0_nom:.0f}, {I0_nom:.0f})')
    axes[1].set_xlabel('Susceptibles S(t)')
    axes[1].set_ylabel('Infectados I(t)')
    axes[1].set_title('Plano de Fase Epidemiológico (S vs I)')
    axes[1].grid(True)
    axes[1].legend()
    
    plt.tight_layout()
    plt.savefig("demo_sir_continuous_phase.png")
    print("Gráfica de dinámica y plano de fase guardada como 'demo_sir_continuous_phase.png'")
    
    # C. Muestreo de Ensamble para la Curva de Infectados (Salida "Infectados")
    model_sir.set_output("Infectados")
    M_sir = design_matrix(model_sir, n=30, method="latin_hypercube", seed=42)
    Y_sir_ens = model_sir.evaluate_batch(M_sir, model_sir.domain) # shape (30, 1, len(domain)) o (30, len(domain))
    if Y_sir_ens.ndim == 3:
        Y_sir_ens = Y_sir_ens.squeeze(axis=1)
        
    plt.figure(figsize=(9, 4.5))
    for i in range(30):
        plt.plot(model_sir.domain, Y_sir_ens[i], color='pink', alpha=0.6, linewidth=0.8)
    plt.plot(model_sir.domain, I_nom, color='tab:red', linewidth=2.5, label='Nominal')
    plt.xlabel('Tiempo (Días)')
    plt.ylabel('Número de Infectados I(t)')
    plt.title('Variabilidad de la Curva de Infectados frente a Incertidumbre en Factores')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("demo_sir_continuous_ensemble.png")
    print("Gráfica del ensamble SIR guardada como 'demo_sir_continuous_ensemble.png'")
    
    # --- 3. PREGUNTAS Y EJERCICIOS PARA EL ESTUDIANTE ---
    print("\n=========================================================================")
    print("  PREGUNTAS DE ANÁLISIS Y EJERCICIOS DOCENTES  ")
    print("=========================================================================")
    print("1. ¿Cómo cambia el tiempo de pico epidémico (t_max) cuando beta aumenta un 20%?")
    print("2. Si la condición inicial I0 cambia de 1 a 50, ¿afecta el valor máximo del pico I_max?")
    print("3. En el modelo logístico discreto, ¿qué sucede con la trayectoria si r = 2.2?")
    print("=========================================================================\n")

if __name__ == "__main__":
    run_demo()
