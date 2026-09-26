"""
Verifica que MAREA2_sin_colonias sigue dando EXACTAMENTE los resultados del
protocolo guardado, despues de tocar marea2.py (Fase 0 del plan).

Recalcula D10 / 5k FES, las primeras --runs corridas de las 71 funciones, por
la misma ruta que el protocolo (experimento_config._run_single) y compara bit a
bit contra resultados_experimento/MAREA2_sin_colonias_D10_FES5000.json.

Solo tiene sentido en el Mac mini: en otra maquina las diferencias de BLAS
hacen fallar la comparacion aunque el codigo sea el mismo.

Uso:
    /ruta/al/python3 verificar_sin_colonias.py            # 5 corridas
    /ruta/al/python3 verificar_sin_colonias.py --runs 10 --workers 20
"""
import argparse, json, os, sys
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from experimento_config import (_run_single, get_problems_for_dim,  # noqa: E402
                                POP_SIZE, RESULTS_DIR)

ALG, DIM, MAXFES = "MAREA2_sin_colonias", 10, 5000


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", type=int, default=5)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 4))
    a = ap.parse_args()

    ref = json.load(open(os.path.join(RESULTS_DIR, f"{ALG}_D{DIM}_FES{MAXFES}.json")))
    T = max(1, MAXFES // POP_SIZE)
    pids = [p["pid"] for p in get_problems_for_dim(DIM)]
    jobs = [(ALG, pid, DIM, r, POP_SIZE, T) for pid in pids for r in range(a.runs)]

    with Pool(a.workers) as pool:
        vals = pool.map(_run_single, jobs)

    malas = [(i, r, v, ref[f"pid{i}_run{r}"])
             for (_, i, _, r, _, _), v in zip(jobs, vals) if v != ref[f"pid{i}_run{r}"]]
    print(f"{len(jobs) - len(malas)}/{len(jobs)} corridas identicas bit a bit")
    for i, r, v, w in malas[:10]:
        print(f"  pid{i}_run{r}: nuevo {v!r}  guardado {w!r}")
    sys.exit(1 if malas else 0)


if __name__ == "__main__":
    main()
