"""
Harness de ablación para MAREA-2
================================

Fase de screening (barata) antes del protocolo completo de 51 corridas:

    python ablacion_marea2.py run  --dims 10 50 --maxfes 50000 --runs 15 --workers 8
    python ablacion_marea2.py run  --dims 10 50 --maxfes 5000  --runs 15 --workers 8
    python ablacion_marea2.py report --dims 10 50 --maxfes 50000
    python ablacion_marea2.py ratio  --dims 10 50 --fes_a 5000 --fes_b 50000

Guarda JSON incremental en RESULTS_DIR/ablacion/ con el mismo formato que
experimento_config (claves 'fun{i}_run{r}'), asi que se puede interrumpir y
reanudar, y los JSON de MAREA_REL / ABC / etc. ya existentes se cargan tal cual
para comparar (se usan las primeras `runs` semillas).

Criterio de aceptacion de un componente (ver MAREA2_especificacion.md, sec. 5):
  el componente se queda si  (a) FULL gana a 'sin_<componente>' en B/E/W sobre
  la clase de funciones que ataca  y  (b) no pierde en unimodales.
"""
from __future__ import annotations
import argparse, json, os, sys, time
from multiprocessing import Pool
import numpy as np

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PROJECT_DIR)

from experimento_config import get_problems_for_dim, POP_SIZE, RESULTS_DIR   # noqa: E402
from marea2 import marea2, config, ABLATION_CONFIGS, A2_CONFIGS              # noqa: E402

TODAS = ABLATION_CONFIGS + [c for c in A2_CONFIGS if c not in ABLATION_CONFIGS]

try:
    from scipy.stats import wilcoxon, rankdata, friedmanchisquare
except ImportError:                                                          # pragma: no cover
    wilcoxon = rankdata = friedmanchisquare = None

ABL_DIR = os.path.join(RESULTS_DIR, 'ablacion')
COMPETIDORES = ['LSHADE', 'JSO', 'MAREA_REL', 'MAREABIO_REL', 'ABC', 'DE', 'FA', 'GWO', 'PSO', 'WOA']

# --------------------------------------------------------------------------
# Clasificacion de funciones (por nombre; fallback por indice en la lista)
# --------------------------------------------------------------------------
def clase(nombre: str, idx: int) -> str:
    import re
    m = re.search(r'F(\d+).*?(2014|2017|2022)', nombre) or re.search(r'(2014|2017|2022).*?F(\d+)', nombre)
    if m:
        g = m.groups()
        n, y = (int(g[0]), g[1]) if g[0].isdigit() and len(g[0]) <= 2 else (int(g[1]), g[0])
    else:                       # orden de get_problems_for_dim: 30 (2014) + 29 (2017) + 12 (2022)
        if idx < 30:   n, y = idx + 1, '2014'
        elif idx < 59: n, y = idx - 30 + 1 + (1 if idx - 30 >= 1 else 0), '2017'   # F2 retirada
        else:          n, y = idx - 59 + 1, '2022'
    if y == '2014':
        return 'unimodal' if n <= 3 else 'multimodal' if n <= 16 else 'hibrida' if n <= 22 else 'composicion'
    if y == '2017':
        return 'unimodal' if n <= 3 else 'multimodal' if n <= 10 else 'hibrida' if n <= 20 else 'composicion'
    return 'unimodal' if n == 1 else 'multimodal' if n <= 5 else 'hibrida' if n <= 8 else 'composicion'

CLASES = ['unimodal', 'multimodal', 'hibrida', 'composicion']

# --------------------------------------------------------------------------
# Ejecucion
# --------------------------------------------------------------------------
def _path(cfg, dim, maxfes):
    os.makedirs(ABL_DIR, exist_ok=True)
    return os.path.join(ABL_DIR, f'MAREA2_{cfg}_D{dim}_FES{maxfes}.json')


_PROB_CACHE = {}

def _worker(args):
    cfg, dim, maxfes, i, r = args
    if dim not in _PROB_CACHE:                    # una vez por proceso worker
        _PROB_CACHE[dim] = get_problems_for_dim(dim)
    prob = _PROB_CACHE[dim][i]
    p = config(cfg); p.N0 = POP_SIZE
    out = marea2(prob['f'], prob['bounds'], p, maxfes=maxfes, seed=r)
    assert out['fes'] <= maxfes, f'{cfg} se paso del presupuesto: {out["fes"]} > {maxfes}'
    return f'fun{i}_run{r}', out['best_f']


def run(cfgs, dims, maxfes, runs, workers):
    for dim in dims:
        n_prob = len(get_problems_for_dim(dim))
        for cfg in cfgs:
            path = _path(cfg, dim, maxfes)
            done = json.load(open(path)) if os.path.exists(path) else {}
            todo = [(cfg, dim, maxfes, i, r) for i in range(n_prob) for r in range(runs)
                    if f'fun{i}_run{r}' not in done]
            print(f'[{cfg:18s} D{dim:2d} FES={maxfes:>7,}] pendientes {len(todo):4d} / {n_prob*runs}',
                  flush=True)
            if not todo:
                continue
            t0 = time.time()
            with Pool(workers) as pool:
                for k, (key, val) in enumerate(pool.imap_unordered(_worker, todo, chunksize=4), 1):
                    done[key] = val
                    if k % 50 == 0 or k == len(todo):
                        json.dump(done, open(path, 'w'))
                        print(f'    {k}/{len(todo)}  {time.time()-t0:6.0f}s', flush=True)
            json.dump(done, open(path, 'w'))


# --------------------------------------------------------------------------
# Carga
# --------------------------------------------------------------------------
def load(alg, dim, maxfes, runs, n_prob):
    """Matriz (n_prob x runs) o None. Busca primero en ablacion/, luego en RESULTS_DIR."""
    cands = [_path(alg, dim, maxfes),
             os.path.join(RESULTS_DIR, f'{alg}_D{dim}_FES{maxfes}.json'),
             os.path.join(PROJECT_DIR, f'{alg}_D{dim}_FES{maxfes}.json')]
    for c in cands:
        if os.path.exists(c):
            M = np.full((n_prob, runs), np.nan)
            for k, v in json.load(open(c)).items():
                a, b = k.split('_'); i, r = int(a[3:]), int(b[3:])
                if i < n_prob and r < runs:
                    M[i, r] = v
            return M if np.isfinite(M).mean() > 0.95 else None
    return None


def _ranks_tie(col, tol=1e-8):
    r = rankdata(col); o = np.argsort(col); out = r.copy(); i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and abs(col[o[j+1]] - col[o[i]]) < tol:
            j += 1
        if j > i:
            a = np.mean([r[o[k]] for k in range(i, j+1)])
            for k in range(i, j+1):
                out[o[k]] = a
        i = j + 1
    return out


def _bew(A, B, sel):
    Bc = Ec = Wc = 0
    for q in sel:
        a, b = A[q], B[q]
        m = np.isfinite(a) & np.isfinite(b)
        if np.all(a[m] == b[m]):       # identicos: el componente no se activo
            Ec += 1; continue
        try:
            pv = wilcoxon(a[m], b[m])[1]
        except Exception:
            pv = 1.0
        if not np.isfinite(pv) or pv >= 0.05: Ec += 1
        elif np.nanmean(a) < np.nanmean(b):   Bc += 1
        else:                                 Wc += 1
    return f'{Bc}/{Ec}/{Wc}'


# --------------------------------------------------------------------------
# Reporte
# --------------------------------------------------------------------------
def report(dims, maxfes, runs, cfgs=None):
    cfgs = cfgs or TODAS
    for dim in dims:
        probs = get_problems_for_dim(dim); n = len(probs)
        cls = [clase(pr['name'], i) for i, pr in enumerate(probs)]
        pool = {}
        for cfg in cfgs:
            M = load(cfg, dim, maxfes, runs, n)
            if M is not None: pool[cfg] = M
        for a in COMPETIDORES:
            M = load(a, dim, maxfes, runs, n)
            if M is not None: pool[a] = M
        if 'FULL' not in pool:
            print(f'D{dim}: falta FULL'); continue
        print(f'\n{"="*78}\nD={dim}  MaxFES={maxfes:,}  runs={runs}  algoritmos={len(pool)}')

        # ranking global
        names = list(pool)
        means = np.array([np.nanmean(pool[a], 1) for a in names])
        R = np.column_stack([_ranks_tie(means[:, q]) for q in range(n)])
        mr = R.mean(1); o = np.argsort(mr)
        print('\nRank medio (menor=mejor):')
        for i in o:
            tag = '  <- variante MAREA-2' if names[i] in cfgs else ''
            print(f'  {names[i]:18s} {mr[i]:6.2f}{tag}')
        if friedmanchisquare is not None and len(names) >= 3:
            st, pv = friedmanchisquare(*means)
            print(f'  Friedman: chi2={st:.1f}  p={pv:.2e}')

        # ablacion: referencia vs cada variante, por clase. Las A2_* son una cadena
        # de un cambio desde 'sin_colonias': se comparan contra ella, no contra FULL.
        print(f'\nAblacion — referencia vs variante  (B/E/W de la referencia, Wilcoxon signed-rank, alpha=0.05)')
        print(f'  {"variante":18s}{"ref":>14s}' + ''.join(f'{c:>14s}' for c in CLASES) + f'{"total":>12s}')
        for cfg in cfgs:
            if cfg == 'FULL' or cfg not in pool: continue
            ref = 'sin_colonias' if cfg.startswith('A2_') else 'FULL'
            if ref not in pool: continue
            row = [_bew(pool[ref], pool[cfg], [q for q in range(n) if cls[q] == c]) for c in CLASES]
            tot = _bew(pool[ref], pool[cfg], range(n))
            print(f'  {cfg:18s}{ref:>14s}' + ''.join(f'{r:>14s}' for r in row) + f'{tot:>12s}')

        # FULL vs competidores, por clase
        comp = [a for a in COMPETIDORES if a in pool]
        if comp:
            print(f'\nFULL vs competidores (B/E/W de FULL)')
            print(f'  {"competidor":18s}' + ''.join(f'{c:>14s}' for c in CLASES) + f'{"total":>12s}')
            for a in comp:
                row = [_bew(pool['FULL'], pool[a], [q for q in range(n) if cls[q] == c]) for c in CLASES]
                print(f'  {a:18s}' + ''.join(f'{r:>14s}' for r in row) + f'{_bew(pool["FULL"], pool[a], range(n)):>12s}')


def ratio(dims, fes_a, fes_b, runs, cfgs=None):
    cfgs = cfgs or TODAS
    """Mediana de err(fes_b)/err(fes_a) por clase: cuanto aprovecha cada algoritmo
    el presupuesto extra. ~1 => estancado (sin mecanismo de escape)."""
    for dim in dims:
        probs = get_problems_for_dim(dim); n = len(probs)
        cls = [clase(pr['name'], i) for i, pr in enumerate(probs)]
        print(f'\nD={dim}  mediana err({fes_b:,}) / err({fes_a:,})   (menor = mejor aprovecha presupuesto)')
        print(f'  {"algoritmo":18s}' + ''.join(f'{c:>14s}' for c in CLASES))
        for a in cfgs + COMPETIDORES:
            A = load(a, dim, fes_a, runs, n); B = load(a, dim, fes_b, runs, n)
            if A is None or B is None: continue
            ma, mb = np.nanmean(A, 1), np.nanmean(B, 1)
            row = []
            for c in CLASES:
                q = [i for i in range(n) if cls[i] == c and ma[i] > 0]
                row.append(np.median(mb[q] / ma[q]) if q else np.nan)
            print(f'  {a:18s}' + ''.join(f'{r:14.3f}' for r in row))


# --------------------------------------------------------------------------
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['run', 'report', 'ratio'])
    ap.add_argument('--cfgs', nargs='*', default=None)
    ap.add_argument('--dims', nargs='*', type=int, default=[10, 50])
    ap.add_argument('--maxfes', type=int, default=50000)
    ap.add_argument('--fes_a', type=int, default=5000)
    ap.add_argument('--fes_b', type=int, default=50000)
    ap.add_argument('--runs', type=int, default=15)
    ap.add_argument('--workers', type=int, default=8)
    a = ap.parse_args()
    if a.cmd == 'run':
        run(a.cfgs or ABLATION_CONFIGS, a.dims, a.maxfes, a.runs, a.workers)
    elif a.cmd == 'report':
        report(a.dims, a.maxfes, a.runs, a.cfgs)
    else:
        ratio(a.dims, a.fes_a, a.fes_b, a.runs, a.cfgs)
