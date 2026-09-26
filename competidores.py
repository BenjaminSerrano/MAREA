"""
Algoritmos competidores para comparación con MAREA.

Todos comparten la misma interfaz:
    resultado = algoritmo(f, bounds, seed=0)
    resultado = {"best_f": float, "best_x": np.ndarray, "history": list, "iters": int}

Paradigmas representados:
    Swarm        : PSO (Kennedy & Eberhart, 1995)
    Evolutivo    : DE/rand/1/bin (Storn & Price, 1997)
    Inmunológico : CLONALG (de Castro & Von Zuben, 2002)
    Lobo         : GWO (Mirjalili et al., 2014)
    Distribución : CMA-ES (Hansen & Ostermeier, 2001) vía scipy
    Recocido     : SA — Simulated Annealing clásico
    Firefly      : FA (Yang, 2010)

Parámetros: valores recomendados por los autores originales (sin re-tuning).
Presupuesto: MaxFES controlado externamente; aquí T=500 iteraciones × N agentes.
"""

import numpy as np
from scipy.optimize import minimize
from scipy.special import gamma as gamma_fn


# ══════════════════════════════════════════════════════════════════════════════
# INTERFAZ COMÚN
# ══════════════════════════════════════════════════════════════════════════════

def _parse_bounds(bounds):
    L = np.array([b[0] for b in bounds])
    U = np.array([b[1] for b in bounds])
    return L, U


# ══════════════════════════════════════════════════════════════════════════════
# CONTROL DE PRESUPUESTO (MaxFES)
#
# El protocolo CEC exige que todos los algoritmos gasten el MISMO numero de
# evaluaciones. `run_experiment` fija T = MaxFES // N asumiendo N evaluaciones
# por iteracion — cierto para PSO/DE/GWO/WOA/ABC (1.00-1.01x el presupuesto),
# FALSO para FA (~N^2/2 por iteracion, medido 39.6x) y CLONALG (~3.5x).
#
# `_presupuesto` envuelve la funcion objetivo, cuenta evaluaciones y corta la
# corrida al agotar el presupuesto, devolviendo el mejor valor hasta ese punto
# — que es exactamente la semantica de "best_f a MaxFES" del protocolo.
#
# El presupuesto es (T+1)*N, no T*N: es lo que consumen de hecho los otros
# algoritmos (poblacion inicial + T generaciones), asi que FA y CLONALG quedan
# en pie de igualdad con ellos y no con una definicion mas estricta.
# ══════════════════════════════════════════════════════════════════════════════

class _PresupuestoAgotado(Exception):
    """Se lanza al alcanzar MaxFES; la atrapa el algoritmo y devuelve su mejor."""


def _presupuesto(f, maxfes):
    """Envuelve f contando evaluaciones y registrando el mejor punto visto."""
    est = {"n": 0, "best_f": np.inf, "best_x": None}

    def g(x):
        if est["n"] >= maxfes:
            raise _PresupuestoAgotado
        est["n"] += 1
        v = f(x)
        if v < est["best_f"]:
            est["best_f"] = v
            est["best_x"] = np.asarray(x).copy()
        return v

    return g, est



# ══════════════════════════════════════════════════════════════════════════════
# 1. PSO — Particle Swarm Optimization
#    Kennedy & Eberhart (1995). "Particle swarm optimization."
#    Parámetros: w=0.729, c1=1.494, c2=1.494  (Clerc & Kennedy, 2002)
# ══════════════════════════════════════════════════════════════════════════════

def pso(f, bounds, N=80, T=500, w=0.729, c1=1.494, c2=1.494, seed=0):
    """
    PSO con factor de constricción de Clerc & Kennedy (2002).
    w=0.729, c1=c2=1.494 son los valores canónicos más citados.
    """
    rng = np.random.default_rng(seed)
    L, U = _parse_bounds(bounds)
    n = len(bounds)

    # Inicialización
    X  = rng.uniform(L, U, size=(N, n))
    V  = rng.uniform(-(U-L), (U-L), size=(N, n)) * 0.1
    F  = np.array([f(X[i]) for i in range(N)])

    pbest_x = X.copy()
    pbest_f = F.copy()
    gi      = np.argmin(F)
    gbest_x = X[gi].copy()
    gbest_f = F[gi]
    history = [gbest_f]

    for t in range(1, T + 1):
        r1 = rng.uniform(0, 1, size=(N, n))
        r2 = rng.uniform(0, 1, size=(N, n))

        V = (w * V
             + c1 * r1 * (pbest_x - X)
             + c2 * r2 * (gbest_x - X))
        X = np.clip(X + V, L, U)
        F = np.array([f(X[i]) for i in range(N)])

        # Actualizar best personal
        improved = F < pbest_f
        pbest_x[improved] = X[improved].copy()
        pbest_f[improved] = F[improved]

        # Actualizar best global
        gi = np.argmin(pbest_f)
        if pbest_f[gi] < gbest_f:
            gbest_f = pbest_f[gi]
            gbest_x = pbest_x[gi].copy()
        history.append(gbest_f)

    return {"best_f": gbest_f, "best_x": gbest_x, "history": history, "iters": T}


# ══════════════════════════════════════════════════════════════════════════════
# 2. DE/rand/1/bin — Differential Evolution
#    Storn & Price (1997). "Differential evolution — a simple and efficient
#    heuristic for global optimization over continuous spaces."
#    Parámetros: F=0.8, CR=0.9  (valores originales del paper)
# ══════════════════════════════════════════════════════════════════════════════

def de(f, bounds, N=80, T=500, F_scale=0.8, CR=0.9, seed=0):
    """
    DE/rand/1/bin — variante estándar más usada en literatura.
    F_scale=0.8, CR=0.9: valores recomendados por Storn & Price (1997).
    """
    rng = np.random.default_rng(seed)
    L, U = _parse_bounds(bounds)
    n = len(bounds)

    X  = rng.uniform(L, U, size=(N, n))
    F_val = np.array([f(X[i]) for i in range(N)])
    gi    = np.argmin(F_val)
    gbest_f = F_val[gi]
    gbest_x = X[gi].copy()
    history = [gbest_f]

    for t in range(1, T + 1):
        for i in range(N):
            # Seleccionar 3 índices distintos de i
            idxs = rng.choice([j for j in range(N) if j != i], 3, replace=False)
            a, b, c = idxs

            # Mutación: V = X_a + F*(X_b - X_c)
            mutant = np.clip(X[a] + F_scale * (X[b] - X[c]), L, U)

            # Crossover binomial
            mask  = rng.uniform(0, 1, n) < CR
            j_rand = rng.integers(0, n)
            mask[j_rand] = True
            trial = np.where(mask, mutant, X[i])

            # Selección
            f_trial = f(trial)
            if f_trial <= F_val[i]:
                X[i]     = trial
                F_val[i] = f_trial
                if f_trial < gbest_f:
                    gbest_f = f_trial
                    gbest_x = trial.copy()

        history.append(gbest_f)

    return {"best_f": gbest_f, "best_x": gbest_x, "history": history, "iters": T}


# ══════════════════════════════════════════════════════════════════════════════
# 3. CLONALG — Artificial Immune System
#    de Castro & Von Zuben (2002). "Learning and optimization using the
#    clonal selection principle."
#    Parámetros: n_clones=10, beta=100, d=0.1  (valores del paper original)
# ══════════════════════════════════════════════════════════════════════════════

def clonalg(f, bounds, N=80, T=500, n_clones=10, beta=100, d=0.1, seed=0):
    """
    CLONALG: selección clonal + hipermutación inversamente proporcional al fitness.
    beta controla el número de clones (n_clones_i = floor(beta/i)).
    d: fracción de la población que se reemplaza aleatoriamente cada generación.

    PRESUPUESTO (corregido 2026-09-04): cada generación evalúa ~278 clones más
    los reemplazos, no N — medido en 3.48x MaxFES. Ahora la corrida se corta al
    agotar (T+1)*N evaluaciones. Ver `_presupuesto`.
    """
    rng = np.random.default_rng(seed)
    L, U = _parse_bounds(bounds)
    n = len(bounds)

    fc, est = _presupuesto(f, (T + 1) * N)
    history = []
    t = 0
    try:
        X     = rng.uniform(L, U, size=(N, n))
        F_val = np.array([fc(X[i]) for i in range(N)])
        order = np.argsort(F_val)
        X     = X[order]; F_val = F_val[order]
        history.append(est["best_f"])

        d_replace = max(1, int(d * N))  # cuántos se reemplazan por generación

        for t in range(1, T + 1):
            clones_X = []
            clones_F = []

            for i in range(N):
                # Número de clones inversamente proporcional al rango
                nc = max(1, int(np.floor(beta * N / (i + 1) / N)))
                nc = min(nc, n_clones)  # tope máximo

                # Tasa de hipermutación: inversamente proporcional al fitness relativo
                rank_norm = (i + 1) / N   # 0→mejor, 1→peor
                sigma = np.exp(-rank_norm)  # alta mutación para peores soluciones

                for _ in range(nc):
                    clone = X[i] + rng.normal(0, sigma * (U - L), n)
                    clone = np.clip(clone, L, U)
                    clones_X.append(clone)
                    clones_F.append(fc(clone))

            # Selección de los mejores clones
            clones_X = np.array(clones_X)
            clones_F = np.array(clones_F)
            all_X    = np.vstack([X, clones_X])
            all_F    = np.concatenate([F_val, clones_F])
            order    = np.argsort(all_F)[:N]
            X        = all_X[order]
            F_val    = all_F[order]

            # Reemplazo aleatorio (diversidad)
            replace_idx = rng.choice(N, d_replace, replace=False)
            X[replace_idx]     = rng.uniform(L, U, size=(d_replace, n))
            F_val[replace_idx] = np.array([fc(X[i]) for i in replace_idx])

            # Reordenar
            order  = np.argsort(F_val)
            X      = X[order]; F_val = F_val[order]
            history.append(est["best_f"])
    except _PresupuestoAgotado:
        pass

    return {"best_f": est["best_f"], "best_x": est["best_x"],
            "history": history, "iters": t, "fes": est["n"]}


# ══════════════════════════════════════════════════════════════════════════════
# 4. GWO — Grey Wolf Optimizer
#    Mirjalili et al. (2014). "Grey wolf optimizer."
#    Sin parámetros libres de usuario (a decrece linealmente de 2 a 0).
# ══════════════════════════════════════════════════════════════════════════════

def gwo(f, bounds, N=80, T=500, seed=0):
    """
    GWO: jerarquía alfa-beta-delta-omega. a decrece de 2 a 0 linealmente.
    Implementación directa de Mirjalili et al. (2014) — sin parámetros extra.
    """
    rng = np.random.default_rng(seed)
    L, U = _parse_bounds(bounds)
    n = len(bounds)

    X     = rng.uniform(L, U, size=(N, n))
    F_val = np.array([f(X[i]) for i in range(N)])
    order = np.argsort(F_val)

    alpha_x, alpha_f = X[order[0]].copy(), F_val[order[0]]
    beta_x,  beta_f  = X[order[1]].copy(), F_val[order[1]]
    delta_x, delta_f = X[order[2]].copy(), F_val[order[2]]
    history = [alpha_f]

    for t in range(1, T + 1):
        a = 2 - 2 * t / T   # decrece linealmente de 2 a 0

        for i in range(N):
            r1, r2 = rng.uniform(0,1,n), rng.uniform(0,1,n)
            A1 = 2*a*r1 - a;  C1 = 2*r2
            D_alpha = np.abs(C1*alpha_x - X[i])
            X1 = alpha_x - A1*D_alpha

            r1, r2 = rng.uniform(0,1,n), rng.uniform(0,1,n)
            A2 = 2*a*r1 - a;  C2 = 2*r2
            D_beta = np.abs(C2*beta_x - X[i])
            X2 = beta_x - A2*D_beta

            r1, r2 = rng.uniform(0,1,n), rng.uniform(0,1,n)
            A3 = 2*a*r1 - a;  C3 = 2*r2
            D_delta = np.abs(C3*delta_x - X[i])
            X3 = delta_x - A3*D_delta

            X[i] = np.clip((X1 + X2 + X3) / 3, L, U)

        F_val = np.array([f(X[i]) for i in range(N)])
        order = np.argsort(F_val)

        if F_val[order[0]] < alpha_f:
            alpha_x, alpha_f = X[order[0]].copy(), F_val[order[0]]
        if F_val[order[1]] < beta_f:
            beta_x,  beta_f  = X[order[1]].copy(), F_val[order[1]]
        if F_val[order[2]] < delta_f:
            delta_x, delta_f = X[order[2]].copy(), F_val[order[2]]

        history.append(alpha_f)

    return {"best_f": alpha_f, "best_x": alpha_x, "history": history, "iters": T}


# ══════════════════════════════════════════════════════════════════════════════
# 5. WOA — Whale Optimization Algorithm
#    Mirjalili & Lewis (2016). "The whale optimization algorithm."
#    Sin parámetros libres — b=1 (espiral logarítmica).
# ══════════════════════════════════════════════════════════════════════════════

def woa(f, bounds, N=80, T=500, b=1.0, seed=0):
    """
    WOA: burbuja de red + ataque espiral. b=1 (parámetro de espiral logarítmica).
    """
    rng = np.random.default_rng(seed)
    L, U = _parse_bounds(bounds)
    n = len(bounds)

    X     = rng.uniform(L, U, size=(N, n))
    F_val = np.array([f(X[i]) for i in range(N)])
    gi    = np.argmin(F_val)
    gbest_x, gbest_f = X[gi].copy(), F_val[gi]
    history = [gbest_f]

    for t in range(1, T + 1):
        a = 2 - 2 * t / T   # decrece de 2 a 0
        a2 = -1 - t / T     # decrece de -1 a -2

        for i in range(N):
            p  = rng.random()
            r1 = rng.uniform(0, 1, n)
            r2 = rng.uniform(0, 1, n)
            A  = 2*a*r1 - a
            C  = 2*r2

            if p < 0.5:
                if np.linalg.norm(A) < 1:
                    # Shrinking encircling
                    D = np.abs(C*gbest_x - X[i])
                    X[i] = np.clip(gbest_x - A*D, L, U)
                else:
                    # Random search
                    rand_idx = rng.integers(0, N)
                    X_rand   = X[rand_idx]
                    D = np.abs(C*X_rand - X[i])
                    X[i] = np.clip(X_rand - A*D, L, U)
            else:
                # Spiral update
                l = rng.uniform(-1, 1, n)
                D = np.abs(gbest_x - X[i])
                X[i] = np.clip(D * np.exp(b*l) * np.cos(2*np.pi*l) + gbest_x, L, U)

        F_val = np.array([f(X[i]) for i in range(N)])
        gi    = np.argmin(F_val)
        if F_val[gi] < gbest_f:
            gbest_f = F_val[gi]; gbest_x = X[gi].copy()
        history.append(gbest_f)

    return {"best_f": gbest_f, "best_x": gbest_x, "history": history, "iters": T}


# ══════════════════════════════════════════════════════════════════════════════
# 6. ABC — Artificial Bee Colony
#    Karaboga & Basturk (2007). "A powerful and efficient algorithm for
#    numerical function optimization: artificial bee colony (ABC) algorithm."
#    Parámetros: limit = N*D/2  (recomendado por Karaboga)
# ══════════════════════════════════════════════════════════════════════════════

def abc(f, bounds, N=80, T=500, seed=0):
    """
    ABC: empleadas + observadoras + exploradoras. limit = N*D/2 (Karaboga, 2007).
    N aquí es el número de fuentes de alimento (= mitad de la colonia total).
    """
    rng = np.random.default_rng(seed)
    L, U = _parse_bounds(bounds)
    n  = len(bounds)
    NF = N // 2   # fuentes de alimento = empleadas = observadoras/2
    limit = NF * n // 2

    X       = rng.uniform(L, U, size=(NF, n))
    F_val   = np.array([f(X[i]) for i in range(NF)])
    trial   = np.zeros(NF, dtype=int)
    gi      = np.argmin(F_val)
    gbest_f = F_val[gi]; gbest_x = X[gi].copy()
    history = [gbest_f]

    for t in range(1, T + 1):
        # ── Empleadas ─────────────────────────────────────────────────────────
        for i in range(NF):
            k   = rng.choice([j for j in range(NF) if j != i])
            d   = rng.integers(0, n)
            phi = rng.uniform(-1, 1)
            candidate    = X[i].copy()
            candidate[d] = np.clip(X[i][d] + phi*(X[i][d] - X[k][d]), L[d], U[d])
            fc = f(candidate)
            if fc < F_val[i]:
                X[i] = candidate; F_val[i] = fc; trial[i] = 0
            else:
                trial[i] += 1

        # ── Observadoras (selección por ruleta) ───────────────────────────────
        fitness  = 1 / (1 + F_val - F_val.min() + 1e-12)
        probs    = fitness / fitness.sum()
        selected = rng.choice(NF, size=NF, p=probs)
        for i in selected:
            k   = rng.choice([j for j in range(NF) if j != i])
            d   = rng.integers(0, n)
            phi = rng.uniform(-1, 1)
            candidate    = X[i].copy()
            candidate[d] = np.clip(X[i][d] + phi*(X[i][d] - X[k][d]), L[d], U[d])
            fc = f(candidate)
            if fc < F_val[i]:
                X[i] = candidate; F_val[i] = fc; trial[i] = 0
            else:
                trial[i] += 1

        # ── Exploradoras (reemplazar agotadas) ────────────────────────────────
        for i in range(NF):
            if trial[i] >= limit:
                X[i]     = rng.uniform(L, U, n)
                F_val[i] = f(X[i])
                trial[i] = 0

        gi = np.argmin(F_val)
        if F_val[gi] < gbest_f:
            gbest_f = F_val[gi]; gbest_x = X[gi].copy()
        history.append(gbest_f)

    return {"best_f": gbest_f, "best_x": gbest_x, "history": history, "iters": T}


# ══════════════════════════════════════════════════════════════════════════════
# 7. FA — Firefly Algorithm
#    Yang (2010). "Firefly algorithm, stochastic test functions and design
#    optimisation."
#    Parámetros: alpha=0.5, beta0=1.0, gamma=1.0  (Yang, 2010)
# ══════════════════════════════════════════════════════════════════════════════

def fa(f, bounds, N=80, T=500, alpha=0.5, beta0=1.0, gamma=1.0, seed=0):
    """
    FA: atracción inversamente proporcional a la distancia y al brillo.
    alpha=0.5, beta0=1.0, gamma=1.0: valores del paper original (Yang, 2010).

    DOS CORRECCIONES (2026-09-04), ambas necesarias para que la comparación
    sea legítima:

    1. ESCALA DE gamma. La atracción es beta0*exp(-gamma*r^2) con r^2 la
       distancia AL CUADRADO. Con gamma=1.0 en el dominio CEC [-100,100], r^2
       tiene mediana 6.6e4 (D=10) y 3.4e5 (D=50), así que exp(-gamma*r^2) hace
       underflow a EXACTAMENTE 0.0 en float64 — incluso en el par más cercano.
       El término de atracción desaparecía y FA quedaba como un paseo aleatorio
       puro. Yang define gamma = 1/Gamma^2 con Gamma la "longitud característica
       del sistema"; aquí Gamma^2 = n*(U-L)^2 (diagonal del dominio al
       cuadrado), lo que deja gamma*r^2 ~ O(1) en cualquier dimensión.
       Es el mismo tipo de defecto que gamma_max/sigma en MAREA: una constante
       calibrada en dominio unitario, inerte en [-100,100].

    2. PRESUPUESTO. El doble bucle evalúa una vez por par (i,j) con j más
       brillante, ~N^2/2 por iteración en vez de N — medido en 39.6x MaxFES.
       Ahora la corrida se corta al agotar (T+1)*N evaluaciones.
    """
    rng = np.random.default_rng(seed)
    L, U = _parse_bounds(bounds)
    n = len(bounds)

    # gamma reescalado a la longitud característica del dominio (Yang, 2010)
    gamma_ef = gamma / (n * float(np.mean(U - L)) ** 2)

    fc, est = _presupuesto(f, (T + 1) * N)
    history = []
    t = 0
    try:
        X     = rng.uniform(L, U, size=(N, n))
        F_val = np.array([fc(X[i]) for i in range(N)])
        history.append(est["best_f"])

        for t in range(1, T + 1):
            # alpha decrece para afinar en iteraciones tardías
            alpha_t = alpha * (0.97 ** t)

            for i in range(N):
                for j in range(N):
                    if F_val[j] < F_val[i]:   # j más brillante que i
                        r2   = np.sum((X[i] - X[j])**2)
                        beta = beta0 * np.exp(-gamma_ef * r2)
                        X[i] = np.clip(
                            X[i]
                            + beta * (X[j] - X[i])
                            + alpha_t * rng.uniform(-0.5, 0.5, n) * (U - L),
                            L, U)
                        F_val[i] = fc(X[i])

            history.append(est["best_f"])
    except _PresupuestoAgotado:
        pass

    return {"best_f": est["best_f"], "best_x": est["best_x"],
            "history": history, "iters": t, "fes": est["n"]}


# ══════════════════════════════════════════════════════════════════════════════
# REGISTRO DE ALGORITMOS
# Permite iterar sobre todos con un bucle simple.
# ══════════════════════════════════════════════════════════════════════════════

COMPETIDORES = {
    "PSO":     {"fn": pso,     "paradigma": "Swarm",       "ref": "Kennedy & Eberhart (1995)"},
    "DE":      {"fn": de,      "paradigma": "Evolutivo",   "ref": "Storn & Price (1997)"},
    "CLONALG": {"fn": clonalg, "paradigma": "Inmunológico","ref": "de Castro & Von Zuben (2002)"},
    "GWO":     {"fn": gwo,     "paradigma": "Swarm",       "ref": "Mirjalili et al. (2014)"},
    "WOA":     {"fn": woa,     "paradigma": "Swarm",       "ref": "Mirjalili & Lewis (2016)"},
    "ABC":     {"fn": abc,     "paradigma": "Swarm",       "ref": "Karaboga & Basturk (2007)"},
    "FA":      {"fn": fa,      "paradigma": "Swarm",       "ref": "Yang (2010)"},
}
