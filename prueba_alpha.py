"""
Prueba exploratoria de alpha_max en CEC2013 (después de la sensibilidad B4).

    python prueba_alpha.py check
    python prueba_alpha.py correr --workers 20
    python prueba_alpha.py verificar

B4 (50k) mostró que alpha_max = 0.9 mejora mucho al valor por defecto 0.5. Pero el 0.5 se
eligió porque 0.99 (AGEO) colapsaba la colonia, y la convergencia prematura se nota a
presupuestos largos: hay que ver si la mejora se sostiene a 500k antes de recalibrar.

Mismo camino que B4 (A2_hgt, CEC2013, 28 funciones, D10 y D50, 15 corridas, seed = run).
  50k : default (tiene que reproducir SENS_default bit a bit), alpha_max 0.95 y 0.99, combo
  500k: default, alpha_max 0.7, 0.9 y 0.99, combo
combo = alpha_max 0.9 + qs_gain 4 + N0 40 (los tres que mejoraron a D50 en B4).

Salida: resultados_experimento/complementarias/ALPHA_{nombre}_D{d}_FES{fes}.json
Todo es reanudable. CEC2013 no se usa para evaluar, así que calibrar aquí no toca el protocolo.
"""
import argparse, math, os, sys, time
import multiprocessing as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import complementarias as cp                         # noqa: E402

RUNS, DIMS = 15, (10, 50)
COMBO = {"alpha_max": 0.9, "qs_gain": 4.0, "N0": 40}
CONFIGS = {
    50_000:  {"default": {}, "alpha_max=0.95": {"alpha_max": 0.95},
              "alpha_max=0.99": {"alpha_max": 0.99}, "combo": COMBO},
    500_000: {"default": {}, "alpha_max=0.7": {"alpha_max": 0.7},
              "alpha_max=0.9": {"alpha_max": 0.9}, "alpha_max=0.99": {"alpha_max": 0.99},
              "combo": COMBO},
}


def _path(name, d, fes):
    return cp._path(f"ALPHA_{name}_D{d}_FES{fes}.json")


def _uno(args):
    # igual que complementarias._sens_uno, con varios cambios y presupuesto variable
    name, cambios, d, fes, i, run = args
    from marea2 import marea2, config
    obj = cp._cec2013()[i](ndim=d)
    def f(x): return obj.evaluate(x) - obj.f_global
    p = config("A2_hgt"); p.N0 = cp.N
    for par, val in cambios.items():
        setattr(p, par, int(val) if par == "N0" else val)
    r = marea2(f, [(-100., 100.)] * d, p, maxfes=fes, seed=run)
    assert r["fes"] <= fes
    return f"fun{i}_run{run}", float(r["best_f"])


def correr(workers):
    nfun = len(cp._cec2013())
    for fes in sorted(CONFIGS):                      # 50k primero: da el chequeo de reproducción
        for d in DIMS:
            for name, cambios in CONFIGS[fes].items():
                p = _path(name, d, fes); data = cp._load(p)
                tareas = [(name, cambios, d, fes, i, r) for i in range(nfun) for r in range(RUNS)
                          if f"fun{i}_run{r}" not in data]
                if not tareas:
                    continue
                print(f"{fes:>7} D{d} {name:16s}: {len(tareas)} corridas", flush=True)
                t0 = time.time()
                with mp.get_context("fork").Pool(workers) as pool:
                    for j, (k, v) in enumerate(pool.imap_unordered(_uno, tareas), 1):
                        data[k] = v
                        if j % 50 == 0 or j == len(tareas):
                            cp._save(p, data)
                print(f"   {time.time() - t0:7.0f} s", flush=True)


def verificar():
    ok = True
    for fes in sorted(CONFIGS):
        for d in DIMS:
            for name in CONFIGS[fes]:
                v = cp._load(_path(name, d, fes))
                m = sum(1 for x in v.values() if not math.isfinite(x) or x < 0)
                if len(v) != 28 * RUNS or m:
                    ok = False; print(f"REVISAR {name} D{d} {fes}: {len(v)}/{28 * RUNS}, {m} malos")
    for d in DIMS:                                   # el default a 50k es la misma corrida que B4
        a, b = cp._load(_path("default", d, 50_000)), cp._load(cp._path(f"SENS_default_D{d}.json"))
        dif = sum(1 for k in b if a.get(k) != b[k])
        if not b or dif:
            ok = False; print(f"REVISAR default 50k D{d}: {dif}/{len(b)} distintos de SENS_default")
        else:
            print(f"OK      default 50k D{d} reproduce SENS_default ({len(b)} corridas)")
    print("\nTODO OK" if ok else "\nHAY QUE REVISAR: mandar esta salida a Claude")
    return ok


def check():
    import opfunu
    assert opfunu.__version__ == "1.0.4", "opfunu tiene que ser 1.0.4"
    assert len(cp._cec2013()) == 28
    sens = cp._load(cp._path("SENS_default_D10.json"))
    assert sens, "falta SENS_default_D10.json (resultado de B4)"
    k, v = _uno(("default", {}, 10, 50_000, 0, 0))
    assert v == sens[k], (k, v, sens[k])
    print(f"  default 50k D10 {k} = {v:.6g}, igual que SENS_default")
    _, v = _uno(("combo", COMBO, 10, 5_000, 0, 0))
    print(f"  combo 5k OK ({v:.6g})")
    print("\nCHEQUEO OK")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["check", "correr", "verificar"])
    ap.add_argument("--workers", type=int, default=20)
    a = ap.parse_args()
    if a.cmd == "check": check()
    elif a.cmd == "correr": correr(a.workers)
    elif a.cmd == "verificar": sys.exit(0 if verificar() else 1)
