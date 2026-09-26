"""
MAREA-bio — formulación de tres fases biológicas (ciclo de Cobetia sp.)

Port fiel de `marea_comparacion.marea_bio` con la firma que espera el runner
de `experimento_config`, sin matplotlib (se importa en cada worker).

Fase I  (planctónica)   : explotación intensa, guiada por mu~ (growth)
Fase II (biofilm)       : exploración estructurada con vecindad, si B~ > umbral
Fase III(desprendimiento): salto de Levy desde el mejor global, con prob. p_det

Diferencias deliberadas respecto de marea_comparacion.py:
  - early_stop DESACTIVADO por defecto. El protocolo CEC exige agotar MaxFES
    (mismo criterio que marea_benchmark.MAREAParams). Ver la nota en ese archivo.
  - sigma_rel: opcional. Ver la nota de escala mas abajo.

Coste por iteracion: N evaluaciones (identico a MAREA-base y a PSO/DE/GWO/WOA).
"""

import numpy as np
from dataclasses import dataclass
from typing import Callable
from scipy.special import gamma as gamma_fn


# ─── Parámetros ──────────────────────────────────────────────────────────────

@dataclass
class MAREABioParams:
    N: int = 80                  # tamaño de la población
    T: int = 500                 # iteraciones máximas
    # Fase I — planctónica
    a1:          float = 0.7
    sigma1:      float = 0.4
    # Fase II — biofilm
    a2:          float = 0.5
    c2:          float = 0.4
    sigma2:      float = 0.15
    # Fase III — desprendimiento
    a3:          float = 0.8     # (declarado en la formulación; no usado aún)
    sigma3:      float = 0.1     # (idem)
    eta:         float = 0.08    # escala del salto de Lévy (× ancho de dominio)
    beta_levy:   float = 1.5     # exponente de Lévy
    p_det_base:  float = 0.20    # prob. base de desprendimiento
    B_threshold: float = 0.45    # umbral de B~ para activar Fase II
    tau:         int   = 25      # constante de decaimiento de mu~

    # ── Escala de los términos de ruido ──────────────────────────────────────
    # sigma1/sigma2 multiplican xi ~ U(-1,1) en unidades ABSOLUTAS. Se calibraron
    # en dominios de ancho ~10.24 (Sphere/Rastrigin), donde sigma1=0.4 es ~3.9%
    # del ancho. En CEC el dominio es [-100,100] (ancho 200) y ese mismo 0.4 pasa
    # a ser el 0.2%: el término de ruido queda inerte.
    # sigma_rel=True reescala sigma1/sigma2 por (U-L)/CAL_WIDTH, preservando la
    # razón sigma/ancho-de-dominio de la calibración. False = port fiel.
    sigma_rel:   bool  = False

    # ── Parada anticipada ────────────────────────────────────────────────────
    # DESACTIVADA por defecto: el protocolo CEC exige agotar MaxFES.
    early_stop:  bool  = False
    tol_iters:   int   = 50      # iteraciones sin mejora para parar
    tol_var:     float = 1e-10   # varianza mínima del fitness poblacional


CAL_WIDTH = 10.24   # ancho del dominio donde se calibraron sigma1/sigma2


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
    Indicadores biológicos derivados del estado del algoritmo.

    mu~ (growth) : calidad de explotación — decae con el estancamiento.
    B~  (biofilm): clustering poblacional — sube cuando la población converge.
    V~  (entropy): dispersión del fitness — alta cuando hay diversidad.
    """
    spread   = np.mean(np.std(X, axis=0))
    mu_tilde = np.exp(-stagnation / p.tau)
    B_tilde  = 1.0 - min(spread / (initial_spread + 1e-12), 1.0)
    std_F    = np.std(F)
    V_tilde  = min(std_F / (initial_std_F + 1e-12), 1.0)
    return mu_tilde, B_tilde, V_tilde


# ─── Algoritmo principal ─────────────────────────────────────────────────────

def marea_bio(f: Callable, bounds: list, p: MAREABioParams = MAREABioParams(),
              seed: int = 42) -> dict:
    rng = np.random.default_rng(seed)
    n = len(bounds)
    L, U = np.array([b[0] for b in bounds]), np.array([b[1] for b in bounds])

    s1, s2 = p.sigma1, p.sigma2
    if p.sigma_rel:
        k = np.mean(U - L) / CAL_WIDTH
        s1, s2 = s1 * k, s2 * k

    X  = rng.uniform(L, U, size=(p.N, n))
    F  = np.array([f(X[i]) for i in range(p.N)])
    gi = np.argmin(F)
    Xg, fg = X[gi].copy(), F[gi]
    history, no_imp = [fg], 0

    initial_spread = np.mean(np.std(X, axis=0)) + 1e-12
    initial_std_F  = np.std(F) + 1e-12
    phase_counts   = [0, 0, 0]

    t = 0
    for t in range(1, p.T + 1):
        mu_t, B_t, V_t = indicadores(X, F, no_imp, initial_spread, initial_std_F, p)
        p_det = p.p_det_base * V_t * (1 + B_t)
        X_new = np.empty_like(X)

        for i in range(p.N):
            r  = rng.random()
            xi = rng.uniform(-1, 1, n)

            if r < p_det:
                # Fase III: desprendimiento — salto de Lévy desde el mejor global
                scale = p.eta * (U - L)
                lv    = levy_flight(n, p.beta_levy, rng)
                X_new[i] = np.clip(Xg + scale * lv, L, U)
                phase_counts[2] += 1

            elif B_t > p.B_threshold:
                # Fase II: biofilm — exploración estructurada con vecindad
                mu_pop   = X.mean(axis=0)
                X_new[i] = np.clip(
                    X[i]
                    + p.a2 * B_t * (Xg - X[i])
                    + p.c2 * V_t * (mu_pop - X[i])
                    + s2 * xi,
                    L, U)
                phase_counts[1] += 1

            else:
                # Fase I: planctónica — explotación intensa
                X_new[i] = np.clip(
                    X[i]
                    + p.a1 * mu_t * (Xg - X[i])
                    + s1 * (1 - B_t) * xi,
                    L, U)
                phase_counts[0] += 1

        X  = X_new
        F  = np.array([f(X[i]) for i in range(p.N)])
        gi = np.argmin(F)

        no_imp = 0 if F[gi] < fg else no_imp + 1
        if F[gi] < fg:
            fg, Xg = F[gi], X[gi].copy()
        history.append(fg)

        if p.early_stop and (no_imp >= p.tol_iters or np.var(F) < p.tol_var):
            break

    return {"best_f": fg, "best_x": Xg, "history": history,
            "iters": t, "phases": phase_counts}
