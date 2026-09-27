"""
Fase 1b · Validación de L-SHADE y jSO contra las tablas publicadas de CEC2017.

CEC2017 D50 con el presupuesto oficial (10 000·D = 500 000 FES exactos, no el
(T+1)·N del protocolo), 25 corridas por función. Criterio del plan: mismo orden
de magnitud por función y mismo patrón de funciones resueltas (error < 1e-8).

Numeración: opfunu omite la F2 oficial (Sum of Different Power, retirada por
inestable), así que su F_k es la F oficial k (k=1) o k+1 (k>=2).

Uso:
    python validar_shade.py run    [--algs LSHADE JSO] [--runs 25] [--workers 8]
    python validar_shade.py report [--algs LSHADE JSO]
"""

import argparse, json, os, sys
import multiprocessing as mp
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import opfunu.cec_based.cec2017 as c17          # noqa: E402
import competidores as comp                     # noqa: E402

D, MAXFES = 50, 500_000
OUT = os.path.join("resultados_experimento", "validacion")
OFICIAL = [1] + list(range(3, 31))              # F oficial de cada F_k de opfunu

# Media publicada del error a D50 (51 corridas), por F oficial.
# TRANSCRITAS DE MEMORIA: confirmar contra los PDF antes de citarlas.
#   JSO:    Brest et al. (2017), CEC 2017, tabla de resultados D=50.
#   LSHADE: tablas de CEC2017 D=50 en comparativas posteriores (p. ej. el
#           propio paper de jSO y LSHADE-cnEpSin, Awad et al. 2017).
PUBLICADO = {
    "JSO": {1: 0, 3: 0, 4: 5.62e1, 5: 1.64e1, 6: 0, 7: 6.65e1, 8: 1.70e1,
            9: 0, 10: 3.14e3, 11: 2.79e1, 12: 1.68e3, 13: 3.06e1, 14: 2.50e1,
            15: 2.39e1, 16: 4.51e2, 17: 2.83e2, 18: 2.43e1, 19: 1.41e1,
            20: 1.40e2, 21: 2.19e2, 22: 1.49e3, 23: 4.30e2, 24: 5.07e2,
            25: 4.81e2, 26: 1.13e3, 27: 5.11e2, 28: 4.60e2, 29: 3.63e2,
            30: 6.01e5},
    "LSHADE": {1: 0, 3: 0, 4: 5.9e1, 5: 1.2e1, 6: 0, 7: 6.3e1, 8: 1.2e1,
               9: 0, 10: 3.2e3, 11: 4.8e1, 12: 2.2e3, 13: 6.0e1, 14: 3.1e1,
               15: 4.2e1, 16: 3.4e2, 17: 2.4e2, 18: 4.6e1, 19: 2.5e1,
               20: 1.3e2, 21: 2.1e2, 22: 2.1e3, 23: 4.3e2, 24: 5.1e2,
               25: 4.8e2, 26: 1.1e3, 27: 5.2e2, 28: 4.6e2, 29: 3.6e2,
               30: 6.5e5},
}
RESUELTO = 1e-8


def _path(alg):
    return os.path.join(OUT, f"VALID_{alg}_CEC2017_D{D}_FES{MAXFES}.json")


def _uno(args):
    alg, k, run = args
    obj = getattr(c17, f"F{k}2017")(ndim=D)
    def f(x): return obj.evaluate(x) - obj.f_global
    r = getattr(comp, alg.lower())(f, [(-100., 100.)] * D, seed=run, maxfes=MAXFES)
    assert r["fes"] == MAXFES, (alg, k, run, r["fes"])
    return f"F{OFICIAL[k - 1]}_run{run}", float(r["best_f"])


def run(algs, runs, workers):
    os.makedirs(OUT, exist_ok=True)
    for alg in algs:
        data = json.load(open(_path(alg))) if os.path.exists(_path(alg)) else {}
        tareas = [(alg, k, r) for k in range(1, 30) for r in range(runs)
                  if f"F{OFICIAL[k - 1]}_run{r}" not in data]
        print(f"{alg}: {len(tareas)} corridas pendientes", flush=True)
        with mp.get_context("fork").Pool(workers) as pool:
            for i, (key, v) in enumerate(pool.imap_unordered(_uno, tareas), 1):
                data[key] = v
                if i % 25 == 0 or i == len(tareas):
                    json.dump(data, open(_path(alg), "w"))
                    print(f"  {alg} {i}/{len(tareas)}", flush=True)


def report(algs):
    for alg in algs:
        if not os.path.exists(_path(alg)):
            print(f"{alg}: sin resultados"); continue
        data = json.load(open(_path(alg)))
        pub = PUBLICADO[alg]
        print(f"\n{alg} — CEC2017 D{D}, {MAXFES:,} FES  (error f−f*; <1e-8 cuenta como 0)")
        print(f"  {'F':>4}{'n':>4}{'media':>11}{'mediana':>11}{'publicada':>11}{'log10 dif':>11}  aviso")
        malos = 0
        for Fo in OFICIAL:
            v = np.array([x for key, x in data.items() if key.startswith(f"F{Fo}_")])
            if not len(v): continue
            v = np.where(v < RESUELTO, 0.0, v)
            m, p = v.mean(), pub[Fo]
            if m == 0 and p == 0:
                d, aviso = 0.0, ""
            elif (m == 0) != (p == 0):
                # patrón de resueltas distinto, salvo que ambos estén por debajo de ~1e-3
                d = np.nan
                aviso = "" if max(m, p) < 1e-3 else "RESUELTA SOLO EN UNO"
            else:
                d = np.log10(m / p)
                aviso = "ORDEN DE MAGNITUD" if abs(d) >= 1 else ""
            malos += bool(aviso)
            print(f"  {Fo:>4}{len(v):>4}{m:>11.3e}{np.median(v):>11.3e}{p:>11.3e}{d:>11.2f}  {aviso}")
        print(f"  funciones con aviso: {malos}/29")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("modo", choices=["run", "report"])
    ap.add_argument("--algs", nargs="+", default=["LSHADE", "JSO"])
    ap.add_argument("--runs", type=int, default=25)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    run(a.algs, a.runs, a.workers) if a.modo == "run" else report(a.algs)
