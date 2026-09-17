"""
Modelo Logístico en Tiempo Discreto utilizando gsua_csb.UserFunctionModel
Este script muestra el flujo completo:
1. Definición de la función en tiempo discreto.
2. Creación del objeto de control `Model` (equivalente a la tabla T de MATLAB).
3. Simulación y evaluación con `model.evaluate()` (equivalente a gsua_eval).
"""
import numpy as np
import matplotlib.pyplot as plt
from gsua_csb import UserFunctionModel

def discrete_logistic_func(params, domain):
    """
    Función del modelo de crecimiento logístico discreto.
    N(t+1) = N(t) + r * N(t) * (1 - N(t) / K)
    """
    r = params[0]
    K = params[1]
    N0 = params[2] if len(params) >= 3 else 10.0
    
    num_steps = len(domain)
    N = np.zeros(num_steps, dtype=float)
    N[0] = N0
    
    for t in range(num_steps - 1):
        N[t + 1] = N[t] + r * N[t] * (1.0 - N[t] / K)
        if N[t + 1] < 0:
            N[t + 1] = 0.0
            
    return N

def create_discrete_logistic_model(
    r_range=(0.05, 0.50),
    K_range=(100.0, 1000.0),
    N0_range=(5.0, 50.0),
    t_max=50
):
    """
    Construye y retorna una instancia de UserFunctionModel (Objeto de Control del Sistema).
    """
    domain = np.arange(0, t_max + 1)
    names = ["r", "K", "N0"]
    param_ranges = np.array([r_range, K_range, N0_range])
    
    model = UserFunctionModel(
        func=discrete_logistic_func,
        names=names,
        range=param_ranges,
        domain=domain,
        output_names=["Poblacion"]
    )
    return model

if __name__ == "__main__":
    print("=== MODELO LOGÍSTICO DISCRETO EN PYTHON (gsua_csb) ===\n")
    
    # 1. Crear el objeto de control del modelo (equivalente a la tabla T)
    model = create_discrete_logistic_model()
    print("1. Objeto de control Model creado:")
    print(f"   - Nombres de factores: {model.names}")
    print(f"   - Rangos de variación:\n{model.range}")
    print(f"   - Punto nominal por defecto: {model.nominal}\n")
    
    # 2. Evaluar el modelo con model.evaluate (equivalente a gsua_eval)
    y_nominal = np.squeeze(model.evaluate(model.nominal, model.domain))
    
    # Evaluación con parámetros alternativos
    test_params = np.array([0.40, 800.0, 10.0])
    y_test = np.squeeze(model.evaluate(test_params, model.domain))
    
    # 3. Graficar simulaciones
    plt.figure(figsize=(8, 4))
    plt.plot(model.domain, y_nominal, 'o-', color='tab:blue', label='Nominal (r=0.275, K=550)')
    plt.plot(model.domain, y_test, 's-', color='tab:red', label='Prueba (r=0.40, K=800)')
    plt.xlabel('Paso de Tiempo (t)')
    plt.ylabel('Población N(t)')
    plt.title('Simulación con model.evaluate() [Equivalente a gsua_eval]')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("standalone_discrete_logistic.png")
    print("Gráfica guardada como 'standalone_discrete_logistic.png'")
