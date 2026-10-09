"""
Etapa B del plan de escritura: experimentos complementarios del paper OE1.

    python complementarias.py check
    python complementarias.py conv   --workers 20 [--runs 51]
    python complementarias.py random --workers 20
    python complementarias.py tiempo
    python complementarias.py sens   --workers 20 [--runs 15]
    python complementarias.py verificar

conv    B1 · Convergencia: mejor error hasta ahora en función de las EVALUACIONES, para
        los 10 algoritmos del paper, en 8 funciones de CEC2017 (2 por clase), D10 y D50,
        500k, con las semillas del protocolo (seed = run). El algoritmo corre por el mismo
        camino que el protocolo (experimento_config._import_alg) y solo se envuelve f con
        un contador, así que el último punto de cada curva tiene que coincidir con el JSON
        del protocolo: `verificar` lo comprueba.
random  B2 · Búsqueda aleatoria uniforme (línea base que pide Vermetten et al. 2024) en las
        71 funciones × 51 corridas × D10/D50. Una corrida de 500 080 evaluaciones da los tres
        presupuestos del protocolo como prefijos: (T+1)·N = 5 040 / 50 080 / 500 080, igual
        que los competidores. Escribe RANDOM_D{d}_FES{5000,50000,500000}.json con el formato
        del protocolo, para que analisis_paper.py lo lea como un algoritmo más.
tiempo  B3 · Tiempo de cómputo con el procedimiento de CEC2017 (T0, T1, T2 en F18 oficial,
        200 000 evaluaciones, 5 corridas), D10 y D50, secuencial. T1 se toma antes de cada
        algoritmo y se usa el mínimo. Hay que correrlo sin otras celdas pesadas en marcha.
sens    B4 · Sensibilidad de MAREA-2 (A2_hgt) a 6 parámetros, uno a la vez, en CEC2013
        (28 funciones, suite DISTINTA de las de evaluación), D10 y D50, 50k, 15 corridas.
        Incluye los dos números medidos de Cobetia (x_mid y nu).

Salida: resultados_experimento/complementarias/ (y los RANDOM_*.json en resultados_experimento/).
Todo es reanudable: lo ya guardado no se repite.
"""
import argparse, json, math, os, sys, time
import multiprocessing as mp
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import experimento_config as ec                      # noqa: E402

OUT = os.path.join(ec.RESULTS_DIR, "complementarias")
N = ec.POP_SIZE                                      # 80
ALGS = ["MAREA2_A2_hgt", "LSHADE", "JSO", "DE", "PSO", "ABC", "GWO", "WOA", "FA", "CLONALG"]

# B1: 2 funciones por clase de CEC2017 (numeración oficial; opfunu omite la F2 oficial,
# así que pid = 30 + k - 1 con k = 1 para F1 y k = oficial - 1 para las demás).
CONV_FUN = {30: "F1", 31: "F3", 33: "F5", 38: "F10", 41: "F13", 45: "F17", 49: "F21", 54: "F26"}
CONV_FES = 500_000
# puntos de registro: malla geométrica + los presupuestos del protocolo
CHECK = sorted(set(np.unique(np.round(np.geomspace(N, 500_080, 120))).astype(int).tolist()
                   + [4_960, 5_040, 50_000, 50_080, 500_000, 500_080]))

# B4: un parámetro a la vez alrededor de A2_hgt (valor por defecto entre paréntesis)
SENS = {
    "gamma_max":   [0.4, 0.6, 1.0, 1.2],          # (0.8)
    "alpha_max":   [0.3, 0.4, 0.7, 0.9],          # (0.5)
    "qs_gain":     [0.0, 1.0, 4.0, 8.0],          # (2.0); 0 apaga la regla 1/5
    "N0":          [40, 60, 120, 160],            # (80)
    "student_nu":  [1.5, 2.5, 5.0, 30.0],         # (3.05, Stehlik et al. 2025)
    "cob_xmid":    [2.5, 3.0, 3.9, 4.3],          # (3.429, Stehlik et al. 2025)
}
SENS_FES, SENS_DIMS = 50_000, (10, 50)


def _path(name):
    os.makedirs(OUT, exist_ok=True)
    return os.path.join(OUT, name)


def _load(path):
    return json.load(open(path)) if os.path.exists(path) else {}


def _save(path, data):
    tmp = path + ".tmp"
    json.dump(data, open(tmp, "w"))
    os.replace(tmp, path)


class Contador:
    """Envuelve f: cuenta evaluaciones y guarda el mejor valor en cada punto de CHECK."""
    def __init__(self, f, checks):
        self.f, self.checks, self.n, self.best, self.j = f, checks, 0, math.inf, 0
        self.curva = []
        self.bounds = getattr(f, "bounds", None)

    def __call__(self, x):
        v = self.f(x)
        self.n += 1
        if v < self.best:
            self.best = float(v)
        while self.j < len(self.checks) and self.n == self.checks[self.j]:
            self.curva.append(self.best); self.j += 1
        return v


def _funcion(pid, ndim):
    """La misma f que construye experimento_config._run_single."""
    entry = ec.CEC_REGISTRY[pid]
    eff = 20 if (ndim == 50 and 50 not in entry["dims"]) else ndim
    obj = entry["class"](ndim=eff)
    fg = obj.f_global
    def f(x): return obj.evaluate(x) - fg
    f.bounds = [(-100., 100.)] * eff
    return f


# ── B1 · convergencia ─────────────────────────────────────────────────────────
def _conv_uno(args):
    alg, d, pid, run = args
    fn = ec._import_alg(alg)
    g = Contador(_funcion(pid, d), CHECK)
    r = fn(g, g.bounds, N=N, T=CONV_FES // N, seed=run)
    curva = g.curva + [None] * (len(CHECK) - len(g.curva))   # puntos más allá del gasto real
    return alg, d, f"pid{pid}_run{run}", {"curva": curva, "fes": g.n, "best_f": float(r["best_f"]),
                                          "best_seen": g.best}


def conv(workers, runs):
    for d in ec.DIMS:
        for alg in ALGS:
            p = _path(f"CONV_{alg}_D{d}.json")
            data = _load(p); data["_checks"] = CHECK
            tareas = [(alg, d, pid, r) for pid in CONV_FUN for r in range(runs)
                      if f"pid{pid}_run{r}" not in data]
            print(f"conv {alg:14s} D{d}: {len(tareas)} corridas pendientes", flush=True)
            if not tareas:
                continue
            t0 = time.time()
            with mp.get_context("fork").Pool(workers) as pool:
                for i, (_, _, k, v) in enumerate(pool.imap_unordered(_conv_uno, tareas), 1):
                    data[k] = v
                    if i % 20 == 0 or i == len(tareas):
                        _save(p, data)
            print(f"   {time.time() - t0:7.0f} s", flush=True)


# ── B2 · búsqueda aleatoria ───────────────────────────────────────────────────
RS_FES = {5000: 5_040, 50000: 50_080, 500000: 500_080}


def _rs_uno(args):
    d, pid, run = args
    f = _funcion(pid, d)
    lo, hi = np.array(f.bounds).T
    rng = np.random.default_rng(run)
    best, n, out = math.inf, 0, {}
    for fes, tope in sorted(RS_FES.items(), key=lambda kv: kv[1]):
        while n < tope:
            X = lo + rng.random((N, lo.size)) * (hi - lo)
            for x in X[: tope - n]:
                v = f(x); n += 1
                if v < best: best = float(v)
        out[fes] = best
    return d, f"pid{pid}_run{run}", out


def random_search(workers):
    for d in ec.DIMS:
        paths = {fes: ec.result_path("RANDOM", d, fes) for fes in RS_FES}
        data = {fes: _load(p) for fes, p in paths.items()}
        tareas = [(d, pid, r) for pid in range(len(ec.CEC_REGISTRY)) for r in range(ec.N_RUNS)
                  if any(f"pid{pid}_run{r}" not in data[fes] for fes in RS_FES)]
        print(f"random D{d}: {len(tareas)} corridas pendientes", flush=True)
        if not tareas:
            continue
        t0 = time.time()
        with mp.get_context("fork").Pool(workers) as pool:
            for i, (_, k, v) in enumerate(pool.imap_unordered(_rs_uno, tareas), 1):
                for fes in RS_FES:
                    data[fes][k] = v[fes]
                if i % 200 == 0 or i == len(tareas):
                    for fes in RS_FES:
                        _save(paths[fes], data[fes])
                    print(f"   {i}/{len(tareas)}  {time.time() - t0:7.0f} s", flush=True)


# ── B3 · tiempo (procedimiento de CEC2017) ────────────────────────────────────
def _t0():
    t = time.perf_counter()
    x = 0.55
    for i in range(1, 1_000_001):
        x = 0.55 + i; x = x + x; x = x / 2; x = x * x; x = math.sqrt(x); x = math.log(x)
        x = math.exp(x); x = x / (x + 2)
    return time.perf_counter() - t


T_FES, T_RUNS = 200_000, 5


def _t1(f, X):
    t = time.perf_counter()
    for x in X: f(x)
    return time.perf_counter() - t


def tiempo():
    # T1 se mide una vez antes de cada algoritmo y T1_D{d} es el mínimo de esas 10 tomas:
    # con una sola toma T1 salió inflado en el mini y (T2-T1)/T0 dio negativo
    # para 4 algoritmos; en el MacBook T2 de PSO queda justo encima de T1, como debe.
    import opfunu.cec_based.cec2017 as c17
    p = _path("TIEMPO.json"); data = _load(p)
    data["T0"] = min(_t0() for _ in range(3))
    print(f"T0 = {data['T0']:.3f} s", flush=True)
    for d in ec.DIMS:
        obj = c17.F172017(ndim=d)                      # F18 oficial (opfunu omite la F2)
        def f(x): return obj.evaluate(x) - obj.f_global
        rng = np.random.default_rng(0)
        X = rng.uniform(-100, 100, (T_FES, d))
        _t1(f, X[: T_FES // 10])                       # calentamiento, se descarta
        for alg in ALGS:
            k, k1 = f"T2_{alg}_D{d}", f"T1_{alg}_D{d}"
            if k in data and k1 in data:
                continue
            data[k1] = _t1(f, X)
            fn = ec._import_alg(alg)
            ts = []
            for run in range(T_RUNS):
                t = time.perf_counter()
                fn(f, [(-100., 100.)] * d, N=N, T=T_FES // N, seed=run)
                ts.append(time.perf_counter() - t)
            data[k] = ts
            data[f"T1_D{d}"] = min(data[f"T1_{a}_D{d}"] for a in ALGS if f"T1_{a}_D{d}" in data)
            _save(p, data)
            print(f"D{d} {alg:14s} T1 = {data[k1]:6.2f} s   T2 = {np.mean(ts):8.2f} s   "
                  f"(T2-T1)/T0 = {(np.mean(ts) - data[f'T1_D{d}']) / data['T0']:8.1f}", flush=True)
    _save(p, data)


# ── B4 · sensibilidad en CEC2013 ──────────────────────────────────────────────
def _cec2013():
    import inspect
    import opfunu.cec_based.cec2013 as c13
    cls = sorted([(n, o) for n, o in inspect.getmembers(c13, inspect.isclass)
                  if n.startswith("F") and n[1].isdigit()], key=lambda t: int(t[0][1:-4]))
    return [o for _, o in cls]


def _sens_cfgs():
    yield "default", None, None
    for par, vals in SENS.items():
        for v in vals:
            yield f"{par}={v}", par, v


def _sens_uno(args):
    name, par, val, d, i, run = args
    from marea2 import marea2, config
    obj = _cec2013()[i](ndim=d)
    def f(x): return obj.evaluate(x) - obj.f_global
    p = config("A2_hgt"); p.N0 = N
    if par is not None:
        setattr(p, par, int(val) if par == "N0" else val)
    r = marea2(f, [(-100., 100.)] * d, p, maxfes=SENS_FES, seed=run)
    assert r["fes"] <= SENS_FES
    return name, d, f"fun{i}_run{run}", float(r["best_f"])


def sens(workers, runs):
    nfun = len(_cec2013())
    for d in SENS_DIMS:
        for name, par, val in _sens_cfgs():
            p = _path(f"SENS_{name}_D{d}.json"); data = _load(p)
            tareas = [(name, par, val, d, i, r) for i in range(nfun) for r in range(runs)
                      if f"fun{i}_run{r}" not in data]
            if not tareas:
                continue
            print(f"sens {name:18s} D{d}: {len(tareas)} corridas", flush=True)
            with mp.get_context("fork").Pool(workers) as pool:
                for _, _, k, v in pool.imap_unordered(_sens_uno, tareas):
                    data[k] = v
            _save(p, data)


# ── verificación ──────────────────────────────────────────────────────────────
def verificar(runs_conv, runs_sens):
    ok = True
    def mal(msg):
        nonlocal ok; ok = False; print("REVISAR", msg)
    # B1: completo, finito y el final coincide con el protocolo
    for d in ec.DIMS:
        for alg in ALGS:
            c = _load(_path(f"CONV_{alg}_D{d}.json"))
            prot = _load(ec.result_path(alg, d, CONV_FES))
            if not prot:
                print(f"        conv {alg:14s} D{d}: sin JSON del protocolo para comparar el final")
            n = dif = 0
            for pid in CONV_FUN:
                for r in range(runs_conv):
                    k = f"pid{pid}_run{r}"
                    if k not in c: continue
                    n += 1; v = c[k]
                    if not math.isfinite(v["best_seen"]): mal(f"conv {alg} D{d} {k} no finito")
                    if k in prot and abs(v["best_f"] - prot[k]) > 1e-9 * max(1.0, abs(prot[k])):
                        dif += 1
                    if abs(v["best_seen"] - v["best_f"]) > 1e-9 * max(1.0, abs(v["best_f"])):
                        mal(f"conv {alg} D{d} {k}: best_f del algoritmo != mejor evaluado")
            esperado = len(CONV_FUN) * runs_conv
            gastos = [c[k]["fes"] for k in c if k.startswith("pid")]
            if gastos:
                print(f"        conv {alg:14s} D{d}: evaluaciones gastadas {min(gastos)}-{max(gastos)}")
            if n != esperado: mal(f"conv {alg} D{d}: {n}/{esperado} corridas")
            elif dif: mal(f"conv {alg} D{d}: {dif}/{n} finales distintos del protocolo")
            else: print(f"OK      conv {alg:14s} D{d}: {n} corridas, final = protocolo")
    # B2
    for d in ec.DIMS:
        for fes in RS_FES:
            v = _load(ec.result_path("RANDOM", d, fes))
            m = sum(1 for x in v.values() if not math.isfinite(x) or x < 0)
            if len(v) != 71 * ec.N_RUNS or m: mal(f"RANDOM D{d} {fes}: {len(v)} corridas, {m} malos")
            else: print(f"OK      RANDOM_D{d}_FES{fes}: {len(v)} corridas")
    # B3
    t = _load(_path("TIEMPO.json"))
    falt = [f"T{i}_{a}_D{d}" for d in ec.DIMS for a in ALGS for i in (1, 2) if f"T{i}_{a}_D{d}" not in t]
    if falt or "T0" not in t: mal(f"tiempo: faltan {falt or 'T0'}")
    else:
        for d in ec.DIMS:
            t1s = [t[f"T1_{a}_D{d}"] for a in ALGS]
            print(f"        tiempo D{d}: T1 entre {min(t1s):.2f} y {max(t1s):.2f} s, se usa {t[f'T1_D{d}']:.2f}")
        neg = [f"{a} D{d}" for d in ec.DIMS for a in ALGS if np.mean(t[f"T2_{a}_D{d}"]) < t[f"T1_D{d}"]]
        if neg: mal(f"tiempo: T2 < T1 en {', '.join(neg)}")
        else: print("OK      TIEMPO.json completo, T2 > T1 en todos")
    # B4
    sens_ok = True
    for d in SENS_DIMS:
        for name, _, _ in _sens_cfgs():
            v = _load(_path(f"SENS_{name}_D{d}.json"))
            m = sum(1 for x in v.values() if not math.isfinite(x) or x < 0)
            if len(v) != 28 * runs_sens or m:
                sens_ok = False; mal(f"SENS {name} D{d}: {len(v)}/{28 * runs_sens}, {m} malos")
    if sens_ok: print("OK      SENS completo")
    print("\nTODO OK" if ok else "\nHAY QUE REVISAR: mandar esta salida a Claude")
    return ok


# ── chequeo rápido (segundos) ─────────────────────────────────────────────────
def check():
    import opfunu
    assert opfunu.__version__ == "1.0.4", "opfunu tiene que ser 1.0.4"
    assert len(_cec2013()) == 28
    # el contador no cambia el resultado y su mejor valor coincide con best_f
    for alg in ALGS:
        fn = ec._import_alg(alg)
        a = fn(_funcion(33, 10), [(-100., 100.)] * 10, N=N, T=10, seed=3)["best_f"]
        g = Contador(_funcion(33, 10), CHECK)
        b = fn(g, g.bounds, N=N, T=10, seed=3)["best_f"]
        assert a == b, (alg, a, b)
        assert abs(g.best - b) <= 1e-12 * max(1, abs(b)), (alg, g.best, b)
        print(f"  {alg:14s} contador OK  ({g.n} evaluaciones)")
    _, _, v = _rs_uno((10, 33, 0))
    assert v[5000] >= v[50000] >= v[500000] >= 0
    print("  random  OK", {k: round(x, 3) for k, x in v.items()})
    print("  sens    OK", _sens_uno(("default", None, None, 10, 0, 0))[3])
    print("\nCHEQUEO OK")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["check", "conv", "random", "tiempo", "sens", "verificar"])
    ap.add_argument("--workers", type=int, default=20)
    ap.add_argument("--runs", type=int, default=None)
    a = ap.parse_args()
    if a.cmd == "check": check()
    elif a.cmd == "conv": conv(a.workers, a.runs or ec.N_RUNS)
    elif a.cmd == "random": random_search(a.workers)
    elif a.cmd == "tiempo": tiempo()
    elif a.cmd == "sens": sens(a.workers, a.runs or 15)
    elif a.cmd == "verificar":
        sys.exit(0 if verificar(a.runs or ec.N_RUNS, 15) else 1)
