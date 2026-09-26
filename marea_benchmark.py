"""
MAREA - Metaheuristic Algorithm for Regulated Exploration and Adaptation
Formulación genérica basada en ALGO.pdf (formulación AGEO)

Ecuación de actualización:
  X^(t+1)_i = X^t_i + α_t(X^t_g - X^t_i) + β_t(μ^t - X^t_i) + γ_t·ξ^t_i

Verificación sobre funciones de referencia estándar.

Parámetros calibrados (grid search normalizado, 80 configs × 15 semillas):
  N=80, alpha_min=0.1, alpha_max=0.99, beta_max=0.3, gamma_max=0.2
  Mejora sobre configuración inicial: 34.9% (score normalizado 0.65 vs 1.00)
  Fecha de calibración: 2026-05-12
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Callable


# ─── Parámetros del algoritmo ────────────────────────────────────────────────

@dataclass
class MAREAParams:
    N: int   = 80       # tamaño de la población  [calibrado: 40→80]
    T: int   = 500      # iteraciones máximas
    alpha_min: float = 0.1
    alpha_max: float = 0.99   # [calibrado: 0.9→0.99]
    beta_max:  float = 0.3    # [calibrado: 0.5→0.3]
    gamma_max: float = 0.2    # [calibrado: 0.3→0.2]
    # ── Parada anticipada ────────────────────────────────────────────────────
    # DESACTIVADA por defecto. El protocolo CEC exige agotar MaxFES: si el
    # algoritmo abandona antes, compite con una fracción del presupuesto de sus
    # competidores y la comparación deja de ser válida.
    # Poner early_stop=True sólo para reproducir la variante con parada
    # anticipada (ablación).
    early_stop:  bool  = False
    tol_fitness: float = 1e-6   # variación mínima del mejor fitness
    tol_iters:   int   = 30     # iteraciones sin mejora para parar
    tol_var:     float = 1e-8   # varianza poblacional mínima
    # ── Escala del término de diversificación ────────────────────────────────
    # gamma_max multiplica xi ~ U(-1,1) en unidades ABSOLUTAS y se calibró en
    # dominios de ancho ~10.24 (Sphere/Rastrigin), donde 0.2 es el 2% del ancho.
    # En CEC el dominio es [-100,100] y ese mismo 0.2 pasa a ser el 0.1%: el
    # termino de diversificacion queda inerte. gamma_rel=True lo reescala por
    # (U-L)/CAL_WIDTH, preservando la razon gamma/ancho de la calibracion.
    # Es el equivalente de sigma_rel en marea_bio_benchmark.py; sirve para el
    # 2x2 que separa el efecto de la escala del efecto de la formulacion.
    gamma_rel:   bool  = False


# ─── Coeficientes de fase ────────────────────────────────────────────────────

CAL_WIDTH = 10.24   # ancho del dominio donde se calibro gamma_max


def alpha(t: int, T: int, p: MAREAParams) -> float:
    """Intensificación: crece linealmente de α_min a α_max."""
    return p.alpha_min + (p.alpha_max - p.alpha_min) * (t / T)

def beta(t: int, T: int, p: MAREAParams) -> float:
    """Exploración estructurada: decrece de β_max a 0."""
    return p.beta_max * (1 - t / T)

def gamma(t: int, T: int, p: MAREAParams) -> float:
    """Diversificación: decrece de γ_max a 0."""
    return p.gamma_max * (1 - t / T)


# ─── Algoritmo principal ─────────────────────────────────────────────────────

def marea(f: Callable, bounds: list[tuple], p: MAREAParams = MAREAParams(),
          maximize: bool = False, seed: int = 42) -> dict:
    """
    Ejecuta MAREA sobre la función f.

    Args:
        f:        función objetivo f: R^n -> R (minimización por defecto)
        bounds:   lista de (L_j, U_j) por dimensión
        p:        parámetros del algoritmo
        maximize: si True, maximiza f
        seed:     semilla para reproducibilidad

    Returns:
        dict con best_x, best_f, history (mejor fitness por iteración)
    """
    rng = np.random.default_rng(seed)
    n   = len(bounds)
    L   = np.array([b[0] for b in bounds])
    U   = np.array([b[1] for b in bounds])
    sign = -1 if maximize else 1
    g_scale = float(np.mean(U - L) / CAL_WIDTH) if p.gamma_rel else 1.0

    # Inicialización uniforme
    X = rng.uniform(L, U, size=(p.N, n))
    F = np.array([sign * f(X[i]) for i in range(p.N)])

    best_idx = np.argmin(F)
    X_g = X[best_idx].copy()
    f_g = F[best_idx]

    history = [sign * f_g]
    no_improve = 0

    for t in range(1, p.T + 1):
        at = alpha(t, p.T, p)
        bt = beta(t,  p.T, p)
        gt = gamma(t, p.T, p) * g_scale

        mu = X.mean(axis=0)                              # centro de masa
        xi = rng.uniform(-1, 1, size=(p.N, n))          # perturbación

        # Ecuación de actualización (Ec. 1 de ALGO.pdf)
        X = (X
             + at * (X_g - X)       # intensificación
             + bt * (mu - X)        # exploración estructurada
             + gt * xi)             # diversificación

        # Restricción de dominio (Ec. 8)
        X = np.clip(X, L, U)

        # Evaluación
        F = np.array([sign * f(X[i]) for i in range(p.N)])

        # Actualización del mejor global
        best_idx = np.argmin(F)
        if F[best_idx] < f_g:
            improvement = abs(f_g - F[best_idx])
            f_g = F[best_idx]
            X_g = X[best_idx].copy()
            no_improve = 0 if improvement > p.tol_fitness else no_improve + 1
        else:
            no_improve += 1

        history.append(sign * f_g)

        # Criterio de parada anticipada (opt-in — ver MAREAParams.early_stop)
        if p.early_stop:
            pop_var = np.var(F)
            if no_improve >= p.tol_iters or pop_var < p.tol_var:
                break

    return {"best_x": X_g, "best_f": sign * f_g, "history": history, "iters": t}


# ─── Funciones de referencia ─────────────────────────────────────────────────

def sphere(x):
    """Óptimo: f(0,...,0) = 0"""
    return np.sum(x**2)

def rastrigin(x):
    """Óptimo: f(0,...,0) = 0. Altamente multimodal."""
    n = len(x)
    return 10*n + np.sum(x**2 - 10*np.cos(2*np.pi*x))

def rosenbrock(x):
    """Óptimo: f(1,...,1) = 0. Valle estrecho."""
    return np.sum(100*(x[1:] - x[:-1]**2)**2 + (1 - x[:-1])**2)

def schwefel(x):
    """Óptimo: f(420.97,...) ≈ 0. Engañosa: óptimo global lejos del centro."""
    n = len(x)
    return 418.9829*n - np.sum(x * np.sin(np.sqrt(np.abs(x))))


# ─── Ejecución y visualización ───────────────────────────────────────────────

def run_benchmark():
    benchmarks = [
        ("Sphere",     sphere,     [(-5.12,  5.12)] * 10, 0.0),
        ("Rastrigin",  rastrigin,  [(-5.12,  5.12)] * 10, 0.0),
        ("Rosenbrock", rosenbrock, [(-2.048, 2.048)] * 10, 0.0),
        ("Schwefel",   schwefel,   [(-500,   500)]   * 10, 0.0),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle("MAREA — Convergencia en funciones de referencia (n=10)", fontsize=13)

    print(f"\n{'Función':<12} {'Mejor f':<14} {'Óptimo conocido':<16} {'Iteraciones'}")
    print("-" * 58)

    for ax, (name, f, bounds, optimum) in zip(axes.flat, benchmarks):
        result = marea(f, bounds)
        history = result["history"]

        ax.plot(history, linewidth=1.5, color="steelblue")
        ax.axhline(optimum, color="crimson", linestyle="--", linewidth=1, label=f"Óptimo = {optimum}")
        ax.set_title(name)
        ax.set_xlabel("Iteración")
        ax.set_ylabel("Mejor f(x)")
        ax.legend(fontsize=8)
        ax.set_yscale("symlog")
        ax.grid(True, alpha=0.3)

        print(f"{name:<12} {result['best_f']:<14.6f} {optimum:<16} {result['iters']}")

    plt.tight_layout()
    plt.savefig("/Users/benjaminserranobarba/Desktop/Trabajo/Candidatura/marea_convergencia.png",
                dpi=150, bbox_inches="tight")
    plt.show()
    print("\nGráfico guardado en marea_convergencia.png")


if __name__ == "__main__":
    run_benchmark()
