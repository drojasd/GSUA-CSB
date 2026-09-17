"""
DEMO 02: Muestreo del Espacio de Parámetros (Design of Experiments - DoE) en Python
Curso guiado de GSUA-CSB

Compara métodos de muestreo:
- Uniforme Aleatorio ("uniform")
- Hipercubo Latino ("latin_hypercube")
- Secuencias de Sobol ("sobol")
- Demostración de escala logarítmica (log_scale)
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from gsua_csb import design_matrix, UserFunctionModel

# Añadir la ruta del Módulo 01 para reutilizar el modelo SIR continuo
module1_py_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "01_Introduccion_Modelos", "python"))
if module1_py_path not in sys.path:
    sys.path.append(module1_py_path)

from sir_continuous_userdefined import create_sir_continuous_model

def run_demo():
    print("=========================================================================")
    print("  MÓDULO 02: MUESTREO DEL ESPACIO DE PARÁMETROS (PYTHON gsua_csb)")
    print("=========================================================================\n")
    
    # 1. Crear Modelo SIR
    model_sir = create_sir_continuous_model()
    print(f"Factores del Modelo: {model_sir.names}")
    print(f"Límites:\n{model_sir.range}\n")
    
    N_samples = 250
    
    # 2. Generar Matrices de Diseño con diferentes métodos
    M_uni   = design_matrix(model_sir, n=N_samples, method="uniform", seed=42)
    M_lhs   = design_matrix(model_sir, n=N_samples, method="latin_hypercube", seed=42)
    M_sobol = design_matrix(model_sir, n=N_samples, method="sobol", seed=42)
    
    # 3. Visualizar y comparar proyecciones 2D (beta vs gamma)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    
    axes[0].scatter(M_uni[:, 0], M_uni[:, 1], color='tab:red', alpha=0.7, s=15)
    axes[0].set_xlabel('beta')
    axes[0].set_ylabel('gamma')
    axes[0].set_title('Muestreo Uniforme Aleatorio')
    axes[0].grid(True)
    
    axes[1].scatter(M_lhs[:, 0], M_lhs[:, 1], color='tab:green', alpha=0.7, s=15)
    axes[1].set_xlabel('beta')
    axes[1].set_ylabel('gamma')
    axes[1].set_title('Hipercubo Latino (LHS)')
    axes[1].grid(True)
    
    axes[2].scatter(M_sobol[:, 0], M_sobol[:, 1], color='tab:blue', alpha=0.7, s=15)
    axes[2].set_xlabel('beta')
    axes[2].set_ylabel('gamma')
    axes[2].set_title('Secuencia de Sobol')
    axes[2].grid(True)
    
    plt.tight_layout()
    plt.savefig("demo_sampling_comparison.png")
    print("Gráfica comparativa de muestreo guardada como 'demo_sampling_comparison.png'")
    
    # 4. Demostración de Escala Logarítmica (log_scale)
    print("\n--- Demostración de Escala Logarítmica (log_scale) ---")
    
    # Modelo con log_scale activado para el primer parámetro
    def dummy_func(params):
        return params[0] + params[1]
        
    model_log = UserFunctionModel(
        func=dummy_func,
        names=["beta_multidecadal", "gamma_std"],
        range=np.array([[1e-6, 1e-1], [0.05, 0.20]]),
        log_scale=[True, False] # beta varía en escala logarítmica
    )
    
    M_log = design_matrix(model_log, n=500, method="latin_hypercube", seed=42)
    
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].hist(M_log[:, 0], bins=25, color='tab:orange', edgecolor='black')
    axes[0].set_xlabel('beta_multidecadal (Escala Lineal)')
    axes[0].set_ylabel('Frecuencia')
    axes[0].set_title('Distribución en Valores Reales (Concentrada en década superior)')
    axes[0].grid(True)
    
    axes[1].hist(np.log10(M_log[:, 0]), bins=25, color='tab:purple', edgecolor='black')
    axes[1].set_xlabel('log10(beta_multidecadal) (Escala Logarítmica)')
    axes[1].set_ylabel('Frecuencia')
    axes[1].set_title('Distribución en log10 (Uniformemente distribuida)')
    axes[1].grid(True)
    
    plt.tight_layout()
    plt.savefig("demo_log_scale_sampling.png")
    print("Gráfica de muestreo logarítmico guardada como 'demo_log_scale_sampling.png'")
    
    # 5. PREGUNTAS Y EJERCICIOS PARA EL ESTUDIANTE
    print("\n=========================================================================")
    print("  PREGUNTAS DE ANÁLISIS Y EJERCICIOS DOCENTES (MÓDULO 02)  ")
    print("=========================================================================")
    print("1. Al observar la proyección 2D (beta vs gamma), ¿cuál método muestra mayores huecos o aglomeraciones?")
    print("2. ¿Por qué el Hipercubo Latino garantiza que no haya dos muestras en el mismo intervalo unidimensional?")
    print("3. Si un parámetro varía entre 1e-5 y 100, ¿qué ocurre si usamos muestreo lineal vs logarítmico?")
    print("=========================================================================\n")

if __name__ == "__main__":
    run_demo()
