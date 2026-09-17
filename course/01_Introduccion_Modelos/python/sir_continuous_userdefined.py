"""
Modelo Epidémico SIR Continuo utilizando scipy.integrate.solve_ivp y gsua_csb.UserFunctionModel
Este script muestra el flujo completo:
1. Definición del sistema de EDOs.
2. Creación del objeto de control `Model` (equivalente a la tabla T de MATLAB).
3. Simulación y evaluación con `model.evaluate()` (equivalente a gsua_eval).
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from gsua_csb import UserFunctionModel

def _sir_rhs(t, y, beta, gamma):
    S, I, R = y
    dSdt = -beta * S * I
    dIdt =  beta * S * I - gamma * I
    dRdt =  gamma * I
    return [dSdt, dIdt, dRdt]

def sir_continuous_func(params, domain):
    """
    Función de evaluación para el modelo SIR continuo resuelto con solve_ivp.
    """
    beta = params[0]
    gamma = params[1]
    S0 = params[2] if len(params) >= 3 else 990.0
    I0 = params[3] if len(params) >= 4 else 10.0
    R0_pop = params[4] if len(params) >= 5 else 0.0
        
    y0 = [S0, I0, R0_pop]
    
    sol = solve_ivp(
        _sir_rhs,
        t_span=(domain[0], domain[-1]),
        y0=y0,
        args=(beta, gamma),
        t_eval=domain,
        method="RK45",
        rtol=1e-6,
        atol=1e-8
    )
    return sol.y

def create_sir_continuous_model(
    beta_range=(0.0001, 0.0008),
    gamma_range=(0.05, 0.20),
    S0_range=(900.0, 999.0),
    I0_range=(1.0, 50.0),
    t_max=80,
    n_points=200,
    output_selection=None
):
    """
    Construye y retorna una instancia de UserFunctionModel (Objeto de Control del Sistema).
    """
    domain = np.linspace(0, t_max, n_points)
    names = ["beta", "gamma", "S0", "I0"]
    param_ranges = np.array([beta_range, gamma_range, S0_range, I0_range])
    output_names = ["Susceptibles", "Infectados", "Recuperados"]
    
    model = UserFunctionModel(
        func=sir_continuous_func,
        names=names,
        range=param_ranges,
        domain=domain,
        output_names=output_names
    )
    
    if output_selection is not None:
        model.set_output(output_selection)
        
    return model

if __name__ == "__main__":
    print("=== MODELO SIR CONTINUO EN PYTHON (gsua_csb) ===\n")
    
    # 1. Crear el objeto de control del modelo (equivalente a la tabla T)
    model = create_sir_continuous_model()
    print("1. Objeto de control Model creado:")
    print(f"   - Nombres de factores: {model.names}")
    print(f"   - Rangos de variación:\n{model.range}")
    print(f"   - Salidas del modelo: {model.output_names}\n")
    
    # 2. Evaluar el modelo con model.evaluate (equivalente a gsua_eval) para TODAS las salidas
    y_all = model.evaluate(model.nominal, model.domain) # shape (3, 200)
    S, I, R = y_all[0], y_all[1], y_all[2]
    
    # 3. Filtrar únicamente la salida de Infectados
    model.set_output("Infectados")
    y_infectados_only = np.squeeze(model.evaluate(model.nominal, model.domain))
    
    # 4. Graficar simulaciones
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    
    axes[0].plot(model.domain, S, label='Susceptibles S(t)')
    axes[0].plot(model.domain, I, label='Infectados I(t)')
    axes[0].plot(model.domain, R, label='Recuperados R(t)')
    axes[0].set_xlabel('Tiempo (Días)')
    axes[0].set_ylabel('Población')
    axes[0].set_title('Todas las Salidas (model.evaluate())')
    axes[0].grid(True); axes[0].legend()
    
    axes[1].plot(model.domain, y_infectados_only, color='tab:red', label='Infectados I(t)')
    axes[1].set_xlabel('Tiempo (Días)')
    axes[1].set_ylabel('Infectados I(t)')
    axes[1].set_title('Salida Filtrada (set_output("Infectados"))')
    axes[1].grid(True); axes[1].legend()
    
    plt.tight_layout()
    plt.savefig("standalone_sir_continuous.png")
    print("Gráfica guardada como 'standalone_sir_continuous.png'")
