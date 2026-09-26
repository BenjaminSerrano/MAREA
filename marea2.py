"""
MAREA-2 — Metaheurística bacteriana por microcolonias
=====================================================

Implementación de referencia. Cada componente tiene un flag en MAREA2Params
para poder apagarlo en la ablación. Con TODOS los flags apagados el algoritmo
se reduce a AGEO/MAREA_REL (ecuación (1) con gamma en unidades relativas).

Interfaz compatible con experimento_config:

    from marea2 import marea2, MAREA2Params
    p = MAREA2Params()
    def alg_fn(f, bounds, N, T, seed):
        p.N0 = N
        return marea2(f, bounds, p, maxfes=N*T, seed=seed)

El presupuesto se controla en EVALUACIONES (maxfes), no en iteraciones:
con reducción de población el número de iteraciones es > T. El smoke test
debe comprobar r['fes'] <= maxfes, no r['iters'] == T.

Mapeo biológico (ver MAREA2_especificacion.md):
  microcolonias          -> nichos con atractor local (X_c*, mu_c)
  inanición -> dispersión-> reinicio parcial disparado por estancamiento
  quorum sensing         -> adaptación del paso por tasa de éxito (regla 1/5)
  logística de Cobetia   -> schedule de tamaño poblacional (capacidad de carga);
                            forma y parámetros de Stehlik et al. 2025, Tabla 1
  trade-off mu_max/biofilm -> regla 1/5 (éxito alto = planctónico, paso mayor;
                            éxito bajo = sésil, paso menor -> desprendimiento)
  Student-t nu=3.05      -> forma de la perturbación (familia TF del paper)
  fase de latencia       -> fracción inicial de FES de exploración garantizada
  memoria de nutrientes  -> archivo de cuencas agotadas (tabú suave)
  células persistentes   -> elitismo por colonia
"""
from __future__ import annotations
from dataclasses import dataclass, field
import math
import numpy as np


# --------------------------------------------------------------------------
# Parámetros
# --------------------------------------------------------------------------
@dataclass
class MAREA2Params:
    # ---- población -------------------------------------------------------
    N0: int = 80                   # lo fija el runner (POP_SIZE)
    N_min: int = 20                # capacidad de carga final (> D+1 no hace falta: ver diff_axis_mix)
    pop_schedule: str = 'cobetia'  # 'cobetia' | 'gompertz' | 'linear' | 'constant'
    lag: float = 0.05              # fase de latencia: fracción de FES
    # 'cobetia': logística en log-tiempo del paper (Stehlik et al. 2025, Tabla 1)
    #   m(t) = A / (1 + exp((xmid - log t)/s)),  t en horas, ventana del ensayo 63 h
    cob_xmid: float = 3.429
    cob_scal: float = 0.784
    cob_tmax: float = 63.0
    mu_g: float = 3.0              # solo 'gompertz': pendiente normalizada (G(1)~0.999)

    # ---- microcolonias ---------------------------------------------------
    colonies: bool = True
    cells_per_colony: int = 20
    recluster_every: int = 10      # iteraciones entre re-clusterings

    # ---- ecuación de movimiento (AGEO, ec. 1) ----------------------------
    alpha_min: float = 0.10
    alpha_max: float = 0.50        # 0.99 en AGEO colapsa la colonia sobre X_c*
    beta_max: float = 0.10         # 0.3 en AGEO: la atracción al centroide era el mayor freno
    gamma_max: float = 0.80        # multiplica el vector de perturbación
    noise: str = 'diff'            # 'diff' (dif. de vecinas, invariante a rotación) | 'normal' | 'uniform'
    diff_axis_mix: float = 0.3     # solo 'diff': peso de la componente por eje
    axis_dist: str = 'student'     # 'student' | 'normal': forma de la componente por eje
    student_nu: float = 3.05       # nu de la familia TF ajustada en el paper (Tabla 1)
    relative_step: bool = True     # True: gamma*sigma_c ; False: gamma*0.02*(U-L)
    sigma_floor: float = 1e-8      # relativo al ancho del dominio
    greedy: bool = True            # aceptar el movimiento solo si mejora

    # ---- quorum sensing (regla 1/5 por colonia) --------------------------
    quorum: bool = True
    qs_target: float = 0.20
    qs_rate: float = 0.10          # EMA de la tasa de éxito
    qs_gain: float = 2.0
    qs_clip: tuple = (0.25, 4.0)

    # ---- inanición -> dispersión -----------------------------------------
    starvation: bool = True
    limit_base: int = 10           # limit = limit_base + 2*D iteraciones sin mejora
    limit_mode: str = 'dim'        # 'dim': 10+2D | 'budget': min(10+2D, max(10, 0.2*T))
                                   # T = maxfes/N0. Solo cambia las celdas de 5k FES (T=62 ->
                                   # limit 12). Con 'dim' la inanicion nunca se activa en D50/5k
                                   # (110 > 62 iteraciones) y no cambia nada en 50k/500k.
    p_disperse: float = 0.5        # fracción (peores) de la colonia que se dispersa
    levy_beta: float = 1.5
    disp_scale: float = 0.10       # escala del salto de Lévy, fracción de (U-L)

    # ---- memoria de nutrientes (cuencas agotadas) ------------------------
    nutrient_memory: bool = True
    depletion_radius: float = 0.05 # fracción de (U-L), distancia normalizada
    memory_size: int = 50

    # ---- otros -----------------------------------------------------------
    elitism: bool = True
    bound_mode: str = 'reflect'    # 'reflect' | 'clip'
    record_history: bool = False


# --------------------------------------------------------------------------
# Utilidades
# --------------------------------------------------------------------------
def _parse_bounds(bounds, D):
    b = np.asarray(bounds, dtype=float)
    if b.ndim == 1 and b.size == 2:
        lo = np.full(D, b[0]); hi = np.full(D, b[1])
    elif b.ndim == 2 and b.shape == (2, D):
        lo, hi = b[0], b[1]
    elif b.ndim == 2 and b.shape == (D, 2):
        lo, hi = b[:, 0], b[:, 1]
    elif b.ndim == 2 and b.shape[0] == 2:
        lo = np.full(D, b[0, 0]); hi = np.full(D, b[1, 0])
    else:
        raise ValueError(f'bounds no reconocido: shape={b.shape}')
    return lo, hi


def _gompertz(u, lag, mu):
    """Gompertz modificada (Zwietering 1990), normalizada A=1. u in [0,1]."""
    return math.exp(-math.exp(mu * math.e * (lag - u) + 1.0))


def _cobetia(u, xmid, scal, tmax):
    """Logística en log-tiempo ajustada a Cobetia sp. (Stehlik et al. 2025, ec. 1,
    Tabla 1), con t = u*tmax horas y normalizada para que G(1) = 1."""
    def L(t):
        return 1.0 / (1.0 + math.exp((xmid - math.log(max(t, 1e-9))) / scal))
    return L(u * tmax) / L(tmax)


def _levy(rng, size, beta):
    """Vuelo de Lévy, algoritmo de Mantegna."""
    num = math.gamma(1 + beta) * math.sin(math.pi * beta / 2)
    den = math.gamma((1 + beta) / 2) * beta * 2 ** ((beta - 1) / 2)
    sigma_u = (num / den) ** (1 / beta)
    u = rng.normal(0, sigma_u, size)
    v = rng.normal(0, 1, size)
    return u / np.abs(v) ** (1 / beta)


def _bound(X, lo, hi, mode):
    if mode == 'clip':
        return np.clip(X, lo, hi)
    # reflexión
    w = hi - lo
    Y = (X - lo) % (2 * w)
    Y = np.where(Y > w, 2 * w - Y, Y)
    return lo + Y


def _kmeans(X, k, rng, init=None, iters=8):
    """k-means ligero. Devuelve etiquetas (n,)."""
    n = X.shape[0]
    k = max(1, min(k, n))
    if k == 1:
        return np.zeros(n, dtype=int)
    if init is not None and init.shape[0] >= k:
        C = init[:k].copy()
    else:
        C = X[rng.choice(n, k, replace=False)].copy()
    lab = np.zeros(n, dtype=int)
    for _ in range(iters):
        d = ((X[:, None, :] - C[None, :, :]) ** 2).sum(-1)
        new = d.argmin(1)
        if np.array_equal(new, lab) and _ > 0:
            break
        lab = new
        for j in range(k):
            m = lab == j
            if m.any():
                C[j] = X[m].mean(0)
            else:                       # colonia vacía: re-sembrar
                C[j] = X[rng.integers(n)]
    # compactar etiquetas (quitar colonias vacías)
    uniq, lab = np.unique(lab, return_inverse=True)
    return lab


# --------------------------------------------------------------------------
# Algoritmo
# --------------------------------------------------------------------------
def marea2(f, bounds, p: MAREA2Params, maxfes: int, seed: int = 0):
    rng = np.random.default_rng(seed)

    # -- dimensión: se infiere de bounds o de una evaluación de prueba -----
    b = np.asarray(bounds, dtype=float)
    if b.ndim == 2 and b.shape[0] not in (2,):
        D = b.shape[0]
    elif b.ndim == 2 and b.shape[1] > 2:
        D = b.shape[1]
    else:
        D = getattr(f, 'dim', None) or getattr(f, 'D', None)
        if D is None:
            raise ValueError('No puedo inferir la dimensión: pasa bounds con shape (D,2) '
                             'o (2,D), o un f con atributo .dim')
    lo, hi = _parse_bounds(bounds, D)
    width = hi - lo
    sig_floor = p.sigma_floor * width

    def evaluate(X):
        return np.array([float(f(x)) for x in X])

    def axis_noise(shape):
        """Forma de la perturbación por eje. Student-t con la nu del paper
        (colas pesadas), reescalada a varianza unitaria para nu > 2."""
        if p.axis_dist == 'student':
            nu = p.student_nu
            z = rng.standard_t(nu, shape)
            return z * math.sqrt((nu - 2.0) / nu) if nu > 2 else z
        return rng.normal(0, 1, shape)

    # -- inicialización -----------------------------------------------------
    N = p.N0
    X = lo + rng.random((N, D)) * width
    F = evaluate(X); fes = N
    gi = int(F.argmin()); gbest_x = X[gi].copy(); gbest_f = F[gi]

    k = max(1, math.ceil(N / p.cells_per_colony)) if p.colonies else 1
    lab = _kmeans(X, k, rng)
    k = lab.max() + 1
    tau = np.zeros(k, dtype=int)                 # iteraciones sin mejora por colonia
    succ = np.full(k, p.qs_target)               # EMA tasa de éxito por colonia
    cbest_f = np.array([F[lab == c].min() for c in range(k)])
    memory = []                                  # centros de cuencas agotadas
    limit = p.limit_base + 2 * D
    if p.limit_mode == 'budget':
        limit = int(min(limit, max(10, 0.2 * maxfes / p.N0)))
    history = []
    t = 0

    def schedule(u):
        if p.pop_schedule == 'constant':
            return p.N0
        if p.pop_schedule == 'linear':
            return int(round(p.N0 - (p.N0 - p.N_min) * u))
        if p.pop_schedule == 'gompertz':
            return int(round(p.N0 - (p.N0 - p.N_min) * _gompertz(u, p.lag, p.mu_g)))
        return int(round(p.N0 - (p.N0 - p.N_min) * _cobetia(u, p.cob_xmid, p.cob_scal, p.cob_tmax)))

    def recluster(X, F):
        nonlocal lab, k, tau, succ, cbest_f
        n = X.shape[0]
        k_new = max(1, math.ceil(n / p.cells_per_colony)) if p.colonies else 1
        # sembrar con los mejores actuales de cada colonia
        seeds = np.array([X[np.where(lab == c)[0][F[lab == c].argmin()]]
                          for c in range(k) if (lab == c).any()])
        new_lab = _kmeans(X, k_new, rng, init=seeds)
        k_new = new_lab.max() + 1
        # heredar tau/succ de la colonia dominante en cada nuevo cluster
        tau_new = np.zeros(k_new, dtype=int); succ_new = np.full(k_new, p.qs_target)
        for c in range(k_new):
            m = new_lab == c
            if m.any():
                dom = np.bincount(lab[m]).argmax()
                tau_new[c] = tau[dom] if dom < len(tau) else 0
                succ_new[c] = succ[dom] if dom < len(succ) else p.qs_target
        lab, k, tau, succ = new_lab, k_new, tau_new, succ_new
        cbest_f = np.array([F[lab == c].min() for c in range(k)])

    # -- bucle principal ---------------------------------------------------
    while fes < maxfes:
        u = fes / maxfes
        t += 1

        # (1) capacidad de carga: eliminar peores si la población debe reducirse
        N_t = max(p.N_min, schedule(u))
        if N_t < X.shape[0]:
            order = np.argsort(F)
            keep = np.sort(order[:N_t])
            if p.elitism and gi not in keep:
                keep[-1] = gi
            X, F, lab = X[keep], F[keep], lab[keep]
            recluster(X, F)
        elif t % p.recluster_every == 0:
            recluster(X, F)
        n = X.shape[0]

        # (2) coeficientes AGEO tras la fase de latencia
        v = 0.0 if u < p.lag else (u - p.lag) / max(1e-12, 1.0 - p.lag)
        alpha = p.alpha_min + (p.alpha_max - p.alpha_min) * v
        beta = p.beta_max * (1.0 - v)
        gamma_base = p.gamma_max * (1.0 - v) + 0.05 * p.gamma_max  # nunca exactamente 0

        # (3) movimiento por colonia
        Xn = np.empty_like(X)
        for c in range(k):
            idx = np.where(lab == c)[0]
            if idx.size == 0:
                continue
            Xc = X[idx]; Fc = F[idx]
            xstar = Xc[Fc.argmin()]
            mu_c = Xc.mean(0)
            if p.relative_step:
                sigma = Xc.std(0) if idx.size > 1 else 0.02 * width
                sigma = np.maximum(sigma, sig_floor)
            else:
                sigma = 0.02 * width           # equivalente a MAREA_REL (0.2/10.24)
            g = gamma_base
            if p.quorum:
                g *= float(np.clip(math.exp(p.qs_gain * (succ[c] - p.qs_target)),
                                   *p.qs_clip))
            if p.noise == 'diff' and idx.size >= 3:
                # transferencia horizontal: diferencia de dos vecinas distintas
                # (sigue la covarianza de la colonia; invariante a rotación)
                a = rng.integers(0, idx.size, idx.size)
                b = rng.integers(0, idx.size, idx.size)
                same = a == b
                b[same] = (b[same] + 1) % idx.size
                pert = (Xc[a] - Xc[b]) * rng.normal(1.0, 0.3, (idx.size, 1))
                # componente por eje: cubre el subespacio que las diferencias
                # no abarcan cuando la colonia tiene <= D células
                pert += p.diff_axis_mix * sigma * axis_noise(Xc.shape)
                Xn[idx] = Xc + alpha * (xstar - Xc) + beta * (mu_c - Xc) + g * pert
                continue
            if p.noise == 'uniform':
                xi = rng.uniform(-1, 1, Xc.shape) * math.sqrt(3)
            else:
                xi = axis_noise(Xc.shape)
            Xn[idx] = Xc + alpha * (xstar - Xc) + beta * (mu_c - Xc) + g * sigma * xi
        Xn = _bound(Xn, lo, hi, p.bound_mode)

        # (4) evaluación (respetando el presupuesto)
        room = maxfes - fes
        if room < n:
            sel = rng.choice(n, room, replace=False)
        else:
            sel = np.arange(n)
        Fn = np.full(n, np.inf)
        Fn[sel] = evaluate(Xn[sel]); fes += sel.size

        improved = Fn < F
        if p.greedy:
            X = np.where(improved[:, None], Xn, X)
            F = np.where(improved, Fn, F)
        else:
            X[sel] = Xn[sel]; F[sel] = Fn[sel]
            if p.elitism:                    # célula persistente: el élite no se pierde
                if F.min() > gbest_f:
                    w = int(F.argmax()); X[w] = gbest_x; F[w] = gbest_f

        # (5) estado por colonia: éxito, estancamiento
        for c in range(k):
            idx = np.where(lab == c)[0]
            if idx.size == 0:
                continue
            rate = improved[idx].mean()
            succ[c] = (1 - p.qs_rate) * succ[c] + p.qs_rate * rate
            fmin = F[idx].min()
            if fmin < cbest_f[c] - 1e-12 * max(1.0, abs(cbest_f[c])):
                cbest_f[c] = fmin; tau[c] = 0
            else:
                tau[c] += 1

        # (6) inanición -> dispersión
        if p.starvation and fes < maxfes:
            for c in range(k):
                if tau[c] <= limit:
                    continue
                idx = np.where(lab == c)[0]
                if idx.size < 2:
                    tau[c] = 0; continue
                order = idx[np.argsort(F[idx])]
                n_keep = max(1, math.ceil((1 - p.p_disperse) * idx.size))
                disp = order[n_keep:]
                if p.nutrient_memory:
                    memory.append(X[order[0]].copy())
                    if len(memory) > p.memory_size:
                        memory.pop(0)
                # salto de Lévy desde la posición PROPIA (no desde X_g)
                step = _levy(rng, (disp.size, D), p.levy_beta) * p.disp_scale * width
                Xd = _bound(X[disp] + step, lo, hi, p.bound_mode)
                if p.nutrient_memory and memory:
                    M = np.asarray(memory)
                    dn = np.sqrt((((Xd[:, None, :] - M[None, :, :]) / width) ** 2).mean(-1))
                    near = (dn < p.depletion_radius).any(1)
                    if near.any():           # cuenca agotada: re-siembra uniforme
                        Xd[near] = lo + rng.random((near.sum(), D)) * width
                room = maxfes - fes
                if room <= 0:
                    break
                if room < disp.size:
                    disp = disp[:room]; Xd = Xd[:room]
                X[disp] = Xd; F[disp] = evaluate(Xd); fes += disp.size
                tau[c] = 0; succ[c] = p.qs_target
                cbest_f[c] = F[lab == c].min()

        # (7) élite global
        gi = int(F.argmin())
        if F[gi] < gbest_f:
            gbest_f = F[gi]; gbest_x = X[gi].copy()
        if p.record_history:
            history.append((fes, gbest_f, k, n))

    return dict(best_f=float(gbest_f), best_x=gbest_x, iters=t, fes=fes,
                history=history)


# --------------------------------------------------------------------------
# Configuraciones de ablación
# --------------------------------------------------------------------------
def config(name: str) -> MAREA2Params:
    """Devuelve la configuración por nombre. FULL = todo encendido."""
    p = MAREA2Params()
    if name == 'FULL':
        return p
    if name == 'sin_colonias':
        p.colonies = False
    elif name == 'sin_inanicion':
        p.starvation = False
    elif name == 'sin_quorum':
        p.quorum = False
    elif name == 'sin_schedule':      # población constante
        p.pop_schedule = 'constant'
    elif name == 'lineal':            # LPSR tipo L-SHADE
        p.pop_schedule = 'linear'
    elif name == 'gompertz':          # Zwietering en vez de la logística del paper
        p.pop_schedule = 'gompertz'
    elif name == 'ruido_gauss':       # componente por eje gaussiana en vez de Student-t
        p.axis_dist = 'normal'
    elif name == 'sin_paso_relativo':
        p.relative_step = False
    elif name == 'sin_memoria':
        p.nutrient_memory = False
    elif name == 'sin_greedy':
        p.greedy = False
    elif name == 'sin_lag':
        p.lag = 0.0
    elif name == 'ruido_normal':      # perturbacion POR EJE (Student-t, axis_dist por defecto)
        p.noise = 'normal'; p.gamma_max = 1.0   # sin diferencia de vecinas. NO es gaussiana.
    elif name == 'eje_student':       # alias explicito de 'ruido_normal'
        p.noise = 'normal'; p.gamma_max = 1.0
    elif name == 'eje_gauss':         # perturbacion por eje GAUSSIANA (cambia noise y axis_dist)
        p.noise = 'normal'; p.gamma_max = 1.0; p.axis_dist = 'normal'
    # ---- Fase A2: cadena de un cambio a la vez desde 'sin_colonias' ----------
    # Resultado de la ablacion one-at-a-time: colonias perjudican (sobre todo en
    # D50); memoria y lag inertes; perturbacion por eje gana a la diferencia de
    # vecinas en D50; la inanicion no se activa en 5k FES con limit=10+2D.
    # Cada brazo cambia UNA cosa respecto a 'sin_colonias', salvo A2_todo.
    elif name == 'A2_eje_student':            # sin colonias + por eje Student-t (sin dif. vecinas)
        p.colonies = False; p.noise = 'normal'; p.gamma_max = 1.0
    elif name == 'A2_eje_gauss':              # sin colonias + por eje gaussiana
        p.colonies = False; p.noise = 'normal'; p.gamma_max = 1.0; p.axis_dist = 'normal'
    elif name == 'A2_limit':                  # sin colonias + inanicion que se activa en 5k
        p.colonies = False; p.limit_mode = 'budget'
    elif name == 'A2_lite':                   # sin colonias, sin memoria, sin lag (los inertes)
        p.colonies = False; p.nutrient_memory = False; p.lag = 0.0
    elif name == 'A2_gompertz':               # sin colonias + Gompertz en vez de la logistica
        p.colonies = False; p.pop_schedule = 'gompertz'
    elif name == 'A2_todo':                   # todos los cambios juntos (Student-t se conserva)
        p.colonies = False; p.noise = 'normal'; p.gamma_max = 1.0
        p.limit_mode = 'budget'; p.nutrient_memory = False; p.lag = 0.0
    elif name == 'AGEO_REL':          # ~ MAREA_REL: todo apagado
        p.colonies = False; p.starvation = False; p.quorum = False
        p.pop_schedule = 'constant'; p.relative_step = False
        p.nutrient_memory = False; p.greedy = False; p.lag = 0.0
        p.axis_dist = 'normal'
        p.alpha_max = 0.99; p.beta_max = 0.3; p.gamma_max = 0.977  # 0.2*(200/10.24)=3.906 abs = 0.977*0.02*200
        p.noise = 'uniform'; p.bound_mode = 'clip'
    else:
        raise ValueError(f'config desconocida: {name}')
    return p


ABLATION_CONFIGS = ['FULL', 'sin_colonias', 'sin_inanicion', 'sin_quorum',
                    'sin_schedule', 'lineal', 'gompertz', 'sin_paso_relativo',
                    'sin_memoria', 'sin_greedy', 'sin_lag', 'ruido_normal',
                    'ruido_gauss', 'AGEO_REL']

# Fase A2: cadena de un cambio a la vez desde 'sin_colonias' (que ya esta corrida).
#   python ablacion_marea2.py run --cfgs A2_eje_student A2_eje_gauss A2_limit A2_lite A2_gompertz A2_todo \
#          --dims 10 50 --maxfes 5000 --runs 15 --workers 20
#   (y con --maxfes 50000 SIN A2_limit: a 50k su limit es el mismo que 'dim' y da
#    resultados identicos a sin_colonias)
#   python ablacion_marea2.py report --dims 10 50 --maxfes 50000     (muestra ABLATION + A2)
A2_CONFIGS = ['sin_colonias', 'A2_eje_student', 'A2_eje_gauss', 'A2_limit',
              'A2_lite', 'A2_gompertz', 'A2_todo']
