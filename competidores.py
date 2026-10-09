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
    DE adaptativo: L-SHADE (Tanabe & Fukunaga, 2014), jSO (Brest et al., 2017)

Parámetros: valores recomendados por los autores originales (sin re-tuning).
Presupuesto: MaxFES controlado externamente; aquí T=500 iteraciones × N agentes.
L-SHADE y jSO usan su propia población y reciben el presupuesto (T+1)*N.
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
#    Parámetros: n_clones=10, beta=100, d=0.1, rho=5, sigma0=0.1
# ══════════════════════════════════════════════════════════════════════════════

def clonalg(f, bounds, N=80, T=500, n_clones=10, beta=100, d=0.1, rho=5.0, sigma0=0.1,
            seed=0):
    """
    CLONALG: selección clonal + hipermutación inversamente proporcional a la afinidad.
    beta controla el número de clones (n_clones_i = floor(beta/i)).
    d: fracción de la población (las de peor afinidad) que se reemplaza al azar.
    Hipermutación: alpha = exp(-rho * f*), con f* la afinidad normalizada por rango
    (1 = la mejor), y el clon es x + N(0, alpha * sigma0 * (U - L)).

    PRESUPUESTO (corregido 2026-09-04): cada generación evalúa ~278 clones más
    los reemplazos, no N — medido en 3.48x MaxFES. Ahora la corrida se corta al
    agotar (T+1)*N evaluaciones. Ver `_presupuesto`.

    MUTACIÓN (corregida 2026-10-09): era sigma = exp(-rank) * (U - L), con rank 0
    para la mejor. Iba al revés (la mejor mutaba con sigma ~ 0.99 * (U - L)) y la
    escala hacía de cada clon un punto casi al azar: perdía contra la búsqueda
    aleatoria en 53-60 de 71 funciones. Además el reemplazo elegía índices al azar
    (podía tocar a la mejor) en vez de las d peores. rho = 5 queda entre 2.5 (no
    alcanza precisión) y 10 (a 5k no mejora a la búsqueda aleatoria) en esfera,
    Rastrigin y Rosenbrock desplazadas, D10.
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

                # Hipermutación inversamente proporcional a la afinidad
                afin = 1.0 - i / (N - 1)            # 1 = mejor, 0 = peor
                alpha = np.exp(-rho * afin)         # de Castro & Von Zuben (2002)

                for _ in range(nc):
                    clone = X[i] + rng.normal(0, alpha * sigma0 * (U - L), n)
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

            # Reemplazo de las d peores por soluciones al azar (diversidad)
            replace_idx = np.arange(N - d_replace, N)
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
# 8-9. L-SHADE y jSO — DE adaptativo de estado del arte
#    L-SHADE: Tanabe & Fukunaga (2014). "Improving the search performance of
#             SHADE using linear population size reduction." CEC 2014 (1º).
#    jSO:     Brest, Maučec & Bošković (2017). "Single objective real-parameter
#             optimization: algorithm jSO." CEC 2017 (2º; 1º entre los DE).
#
# A diferencia de los 7 clásicos, NO usan N=80 fijo: derivan su población de D
# y la reducen linealmente hasta 4 según las FES consumidas (LPSR). Por eso
# reciben un presupuesto, no (N, T). Con la firma común (N, T) del runner el
# presupuesto es (T+1)*N — el mismo que gastan los clásicos —; `maxfes`
# explícito lo sustituye (validación contra CEC2017 a 10 000·D exactos).
#
# Parámetros de los autores, sin re-tuning:
#                 L-SHADE           jSO
#   N_init        18·D              round(25·ln(D)·√D)
#   N_min         4                 4
#   H             6                 5 (la última celda fija en 0.9)
#   M_F, M_CR     0.5, 0.5          0.3, 0.8
#   archivo       2.6·N             1.0·N
#   p (p-best)    0.11              0.125 → 0.25 lineal en FES
#   mutación      cur-to-pbest/1    cur-to-pbest-w/1 (F_w = 0.7F/0.8F/1.2F)
#   memoria       Lehmer ponderada  (Lehmer ponderada + valor previo) / 2
#   jSO además:   CR >= 0.7 (<25% FES), >= 0.6 (<50%); F <= 0.7 (<60%).
# Límites: si v_j sale del dominio, v_j = (límite + x_ij)/2 (ambos papers).
# El valor terminal de M_CR (⊥ → CR = 0) se representa con NaN.
# ══════════════════════════════════════════════════════════════════════════════

def _shade(f, bounds, maxfes, seed, jso):
    rng = np.random.default_rng(seed)
    L, U = _parse_bounds(bounds)
    n = len(bounds)

    if jso:
        N_init = int(round(25 * np.log(n) * np.sqrt(n)))
        H, r_arc, mf0, mcr0 = 5, 1.0, 0.3, 0.8
    else:
        N_init = 18 * n
        H, r_arc, mf0, mcr0 = 6, 2.6, 0.5, 0.5
    N_min = 4
    M_F = np.full(H, mf0); M_CR = np.full(H, mcr0)
    if jso:
        M_F[-1] = M_CR[-1] = 0.9
    k_max = H - 1 if jso else H          # celdas actualizables de la memoria
    k = 0

    fc, est = _presupuesto(f, maxfes)
    history = []
    g = 0
    try:
        N = N_init
        X = rng.uniform(L, U, size=(N, n))
        Fx = np.array([fc(x) for x in X])
        A = np.empty((0, n))
        history.append(est["best_f"])

        while True:
            g += 1
            frac = est["n"] / maxfes
            idx = np.arange(N)

            # ── parámetros de control desde la memoria ────────────────────────
            r = rng.integers(0, H, N)
            mcr = M_CR[r]
            CR = np.where(np.isnan(mcr), 0.0,
                          np.clip(rng.normal(np.nan_to_num(mcr), 0.1), 0, 1))
            F = M_F[r] + 0.1 * rng.standard_cauchy(N)
            while np.any(F <= 0):
                bad = F <= 0
                F[bad] = M_F[r[bad]] + 0.1 * rng.standard_cauchy(bad.sum())
            F = np.minimum(F, 1.0)
            if jso:
                if frac < 0.25:  CR = np.maximum(CR, 0.7)
                elif frac < 0.5: CR = np.maximum(CR, 0.6)
                if frac < 0.6:   F = np.minimum(F, 0.7)
                Fw = F * (0.7 if frac < 0.2 else 0.8 if frac < 0.4 else 1.2)
                p = 0.125 + 0.125 * frac
            else:
                Fw = F
                p = 0.11

            # ── current-to-pbest(-w)/1 con archivo ────────────────────────────
            n_p = max(2, int(round(p * N)))
            pb = np.argsort(Fx)[rng.integers(0, n_p, N)]
            r1 = rng.integers(0, N, N)
            while np.any(m := r1 == idx):
                r1[m] = rng.integers(0, N, m.sum())
            PA = np.vstack([X, A])
            r2 = rng.integers(0, len(PA), N)
            while np.any(m := (r2 == idx) | (r2 == r1)):
                r2[m] = rng.integers(0, len(PA), m.sum())
            V = X + Fw[:, None] * (X[pb] - X) + F[:, None] * (X[r1] - PA[r2])
            V = np.where(V < L, (L + X) / 2, V)
            V = np.where(V > U, (U + X) / 2, V)

            # ── cruce binomial ────────────────────────────────────────────────
            mask = rng.random((N, n)) < CR[:, None]
            mask[idx, rng.integers(0, n, N)] = True
            Ut = np.where(mask, V, X)
            Fu = np.array([fc(u) for u in Ut])

            # ── selección, archivo y memoria ──────────────────────────────────
            better = Fu < Fx
            if better.any():
                S_CR, S_F = CR[better], F[better]
                w = Fx[better] - Fu[better]; w = w / w.sum()
                A = np.vstack([A, X[better]])
                if np.isnan(M_CR[k]) or S_CR.max() == 0:
                    cr_new = np.nan
                else:
                    cr_new = np.sum(w * S_CR**2) / np.sum(w * S_CR)
                f_new = np.sum(w * S_F**2) / np.sum(w * S_F)
                if jso:
                    cr_new = (cr_new + M_CR[k]) / 2
                    f_new = (f_new + M_F[k]) / 2
                M_CR[k], M_F[k] = cr_new, f_new
                k = (k + 1) % k_max
            keep = Fu <= Fx
            X[keep], Fx[keep] = Ut[keep], Fu[keep]

            # ── reducción lineal de la población (LPSR) ───────────────────────
            N_new = int(round((N_min - N_init) * est["n"] / maxfes + N_init))
            if N_new < N:
                o = np.argsort(Fx)[:N_new]
                X, Fx, N = X[o], Fx[o], N_new
            A_max = int(round(r_arc * N))
            if len(A) > A_max:
                A = A[rng.choice(len(A), A_max, replace=False)]
            history.append(est["best_f"])
    except _PresupuestoAgotado:
        pass

    return {"best_f": est["best_f"], "best_x": est["best_x"],
            "history": history, "iters": g, "fes": est["n"]}


def lshade(f, bounds, N=80, T=500, seed=0, maxfes=None):
    """L-SHADE (Tanabe & Fukunaga, 2014). N y T sólo fijan el presupuesto (T+1)*N."""
    return _shade(f, bounds, maxfes or (T + 1) * N, seed, jso=False)


def jso(f, bounds, N=80, T=500, seed=0, maxfes=None):
    """jSO (Brest et al., 2017). N y T sólo fijan el presupuesto (T+1)*N."""
    return _shade(f, bounds, maxfes or (T + 1) * N, seed, jso=True)


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
    "LSHADE":  {"fn": lshade,  "paradigma": "Evolutivo",   "ref": "Tanabe & Fukunaga (2014)"},
    "JSO":     {"fn": jso,     "paradigma": "Evolutivo",   "ref": "Brest et al. (2017)"},
}
