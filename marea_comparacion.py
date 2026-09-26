"""
Comparación entre dos versiones de MAREA:
  - MAREA-base:  formulación genérica de ALGO.pdf (tres fuerzas fijas)
  - MAREA-bio:   formulación de tres fases biológicas (página web)

Los indicadores biológicos μ~, B~, V~ se derivan del estado del algoritmo
(sin datos de Cobetia sp.) para poder probar sobre funciones de referencia.

Parámetros calibrados (grid search normalizado, 2026-05-12):
  N=80, alpha_min=0.1, alpha_max=0.99, beta_max=0.3, gamma_max=0.2
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from scipy.special import gamma as gamma_fn


# ─── Parámetros ──────────────────────────────────────────────────────────────

@dataclass
class Params:
    N:           int   = 80      # tamaño de la población  [calibrado: 40→80]
    T:           int   = 500     # iteraciones máximas
    # MAREA-base
    alpha_min:   float = 0.1
    alpha_max:   float = 0.99    # [calibrado: 0.9→0.99]
    beta_max:    float = 0.3     # [calibrado: 0.5→0.3]
    gamma_max:   float = 0.2     # [calibrado: 0.3→0.2]
    # MAREA-bio — Fase I (planctónica)
    a1:          float = 0.7
    sigma1:      float = 0.4
    # MAREA-bio — Fase II (biofilm)
    a2:          float = 0.5
    c2:          float = 0.4
    sigma2:      float = 0.15
    # MAREA-bio — Fase III (desprendimiento)
    a3:          float = 0.8
    sigma3:      float = 0.1
    eta:         float = 0.08    # escala del salto de Lévy
    beta_levy:   float = 1.5     # exponente de Lévy
    p_det_base:  float = 0.20    # prob. base de desprendimiento
    B_threshold: float = 0.45    # umbral para activar Fase II
    tau:         int   = 25      # constante de decaimiento para μ~
    # Parada
    tol_iters:   int   = 50
    tol_var:     float = 1e-10


# ─── Utilidades ──────────────────────────────────────────────────────────────

def levy_flight(n: int, beta: float, rng) -> np.ndarray:
    """Salto de Lévy mediante el algoritmo de Mantegna."""
    num = gamma_fn(1 + beta) * np.sin(np.pi * beta / 2)
    den = gamma_fn((1 + beta) / 2) * beta * 2**((beta - 1) / 2)
    sigma_u = (num / den) ** (1 / beta)
    u = rng.normal(0, sigma_u, n)
    v = rng.normal(0, 1, n)
    return u / (np.abs(v) ** (1 / beta))


def indicadores(X, F, stagnation, initial_spread, initial_std_F, p):
    """
    Calcula los tres indicadores biológicos a partir del estado del algoritmo.

    μ~ (growth): calidad de explotación — decae cuando hay estancamiento.
    B~ (biofilm): clustering poblacional — sube cuando la población converge.
    V~ (entropy): varianza del fitness — alta cuando hay diversidad.
    """
    spread = np.mean(np.std(X, axis=0))
    mu_tilde = np.exp(-stagnation / p.tau)
    B_tilde  = 1.0 - min(spread / (initial_spread + 1e-12), 1.0)
    std_F    = np.std(F)
    V_tilde  = min(std_F / (initial_std_F + 1e-12), 1.0)
    return mu_tilde, B_tilde, V_tilde


# ─── MAREA-base (ALGO.pdf) ───────────────────────────────────────────────────

def marea_base(f, bounds, p=Params(), seed=42):
    rng = np.random.default_rng(seed)
    n = len(bounds)
    L, U = np.array([b[0] for b in bounds]), np.array([b[1] for b in bounds])

    X  = rng.uniform(L, U, size=(p.N, n))
    F  = np.array([f(X[i]) for i in range(p.N)])
    gi = np.argmin(F)
    Xg, fg = X[gi].copy(), F[gi]
    history, no_imp = [fg], 0

    for t in range(1, p.T + 1):
        at = p.alpha_min + (p.alpha_max - p.alpha_min) * t / p.T
        bt = p.beta_max  * (1 - t / p.T)
        gt = p.gamma_max * (1 - t / p.T)
        mu = X.mean(axis=0)
        xi = rng.uniform(-1, 1, size=(p.N, n))

        X  = np.clip(X + at*(Xg - X) + bt*(mu - X) + gt*xi, L, U)
        F  = np.array([f(X[i]) for i in range(p.N)])
        gi = np.argmin(F)

        no_imp = 0 if F[gi] < fg else no_imp + 1
        if F[gi] < fg:
            fg, Xg = F[gi], X[gi].copy()
        history.append(fg)

        if no_imp >= p.tol_iters or np.var(F) < p.tol_var:
            break

    return {"best_f": fg, "best_x": Xg, "history": history, "iters": t}


# ─── MAREA-bio (tres fases biológicas) ───────────────────────────────────────

def marea_bio(f, bounds, p=Params(), seed=42):
    rng = np.random.default_rng(seed)
    n = len(bounds)
    L, U = np.array([b[0] for b in bounds]), np.array([b[1] for b in bounds])

    X  = rng.uniform(L, U, size=(p.N, n))
    F  = np.array([f(X[i]) for i in range(p.N)])
    gi = np.argmin(F)
    Xg, fg = X[gi].copy(), F[gi]
    history, no_imp = [fg], 0

    initial_spread = np.mean(np.std(X, axis=0)) + 1e-12
    initial_std_F  = np.std(F) + 1e-12
    phase_counts   = [0, 0, 0]

    for t in range(1, p.T + 1):
        mu_t, B_t, V_t = indicadores(X, F, no_imp, initial_spread, initial_std_F, p)
        p_det = p.p_det_base * V_t * (1 + B_t)
        X_new = np.empty_like(X)

        for i in range(p.N):
            r = rng.random()
            xi = rng.uniform(-1, 1, n)

            if r < p_det:
                # Fase III: desprendimiento — salto de Lévy
                scale = p.eta * (U - L)
                lv    = levy_flight(n, p.beta_levy, rng)
                X_new[i] = np.clip(Xg + scale * lv, L, U)
                phase_counts[2] += 1

            elif B_t > p.B_threshold:
                # Fase II: biofilm — exploración estructurada con vecindad
                mu_pop   = X.mean(axis=0)
                X_new[i] = np.clip(
                    X[i]
                    + p.a2 * B_t  * (Xg - X[i])
                    + p.c2 * V_t  * (mu_pop - X[i])
                    + p.sigma2 * xi,
                    L, U)
                phase_counts[1] += 1

            else:
                # Fase I: planctónica — explotación intensa
                X_new[i] = np.clip(
                    X[i]
                    + p.a1 * mu_t   * (Xg - X[i])
                    + p.sigma1 * (1 - B_t) * xi,
                    L, U)
                phase_counts[0] += 1

        X = X_new
        F = np.array([f(X[i]) for i in range(p.N)])
        gi = np.argmin(F)

        no_imp = 0 if F[gi] < fg else no_imp + 1
        if F[gi] < fg:
            fg, Xg = F[gi], X[gi].copy()
        history.append(fg)

        if no_imp >= p.tol_iters or np.var(F) < p.tol_var:
            break

    return {"best_f": fg, "best_x": Xg, "history": history,
            "iters": t, "phases": phase_counts}


# ─── Funciones de referencia ─────────────────────────────────────────────────

def sphere(x):     return np.sum(x**2)
def rastrigin(x):  return 10*len(x) + np.sum(x**2 - 10*np.cos(2*np.pi*x))
def rosenbrock(x): return np.sum(100*(x[1:]-x[:-1]**2)**2 + (1-x[:-1])**2)
def schwefel(x):   return 418.9829*len(x) - np.sum(x*np.sin(np.sqrt(np.abs(x))))


# ─── Ejecución ───────────────────────────────────────────────────────────────

def run():
    benchmarks = [
        ("Sphere",     sphere,     [(-5.12,  5.12)]*10, 0.0),
        ("Rastrigin",  rastrigin,  [(-5.12,  5.12)]*10, 0.0),
        ("Rosenbrock", rosenbrock, [(-2.048, 2.048)]*10, 0.0),
        ("Schwefel",   schwefel,   [(-500,    500)]*10,  0.0),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle("MAREA-base vs MAREA-bio — Convergencia (n=10)", fontsize=13)

    header = f"{'Función':<12} {'MAREA-base':>14} {'MAREA-bio':>14} {'Óptimo':>10} {'Iters-base':>12} {'Iters-bio':>10}"
    print("\n" + header)
    print("-" * 76)

    for ax, (name, f, bounds, opt) in zip(axes.flat, benchmarks):
        r_base = marea_base(f, bounds)
        r_bio  = marea_bio(f, bounds)

        ax.plot(r_base["history"], lw=1.5, color="steelblue",  label="MAREA-base")
        ax.plot(r_bio["history"],  lw=1.5, color="darkorange", label="MAREA-bio",
                linestyle="--")
        ax.axhline(opt, color="crimson", linestyle=":", lw=1, label=f"Óptimo={opt}")
        ax.set_title(name)
        ax.set_xlabel("Iteración")
        ax.set_ylabel("Mejor f(x)")
        ax.legend(fontsize=8)
        ax.set_yscale("symlog")
        ax.grid(True, alpha=0.3)

        phases = r_bio["phases"]
        total  = sum(phases) or 1
        print(f"{name:<12} {r_base['best_f']:>14.4f} {r_bio['best_f']:>14.4f} "
              f"{opt:>10} {r_base['iters']:>12} {r_bio['iters']:>10}")
        print(f"{'':12}  Fases bio → I:{phases[0]/total:.0%}  "
              f"II:{phases[1]/total:.0%}  III:{phases[2]/total:.0%}")

    plt.tight_layout()
    out = "/Users/benjaminserranobarba/Desktop/Trabajo/Candidatura/marea_comparacion.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.show()
    print(f"\nGráfico guardado en marea_comparacion.png")


if __name__ == "__main__":
    run()
