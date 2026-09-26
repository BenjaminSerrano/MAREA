"""
Huella digital de las 71 funciones CEC.

Evalúa cada función en puntos FIJOS y deterministas (sin RNG) e imprime el
valor con precisión completa. Correr en dos máquinas y comparar las salidas:
cualquier función cuyo valor difiera tiene una DEFINICIÓN distinta entre
ambos entornos, y por tanto los resultados no son comparables.

Uso:
    python3 diagnostico_cec.py > huella_<maquina>.txt
"""
import sys, hashlib
import warnings
warnings.filterwarnings("ignore")
import numpy as np

sys.path.insert(0, ".")
from experimento_config import CEC_REGISTRY

print("=" * 78)
print("HUELLA DIGITAL CEC")
print("=" * 78)
print(f"python  : {sys.version.split()[0]}")
print(f"numpy   : {np.__version__}")
try:
    import opfunu
    print(f"opfunu  : {opfunu.__version__}")
except Exception as e:
    print(f"opfunu  : ??? ({e})")
try:
    cfg = np.show_config(mode="dicts")
    blas = cfg.get("Build Dependencies", {}).get("blas", {})
    print(f"blas    : {blas.get('name','?')} {blas.get('version','?')}")
except Exception:
    print("blas    : (no disponible)")
print("=" * 78)
print(f"{'func':<12}{'dim':>4}  {'f(centro)':>26}  {'f(rampa)':>26}")
print("-" * 78)

h = hashlib.sha256()
problemas = []

for entry in CEC_REGISTRY:
    for ndim in (10, 20 if 50 not in entry["dims"] else 50):
        try:
            obj = entry["class"](ndim=ndim)
            g = obj.f_global
            # Puntos deterministas, sin RNG: el centro del dominio y una rampa.
            centro = np.zeros(ndim)
            rampa = np.linspace(-50.0, 50.0, ndim)
            v1 = float(obj.evaluate(centro)) - g
            v2 = float(obj.evaluate(rampa)) - g
        except Exception as e:
            print(f"{entry['name']:<12}{ndim:>4}  ERROR: {type(e).__name__}: {e}")
            problemas.append((entry["name"], ndim, f"excepcion {type(e).__name__}"))
            continue

        malo = ""
        if not np.isfinite(v1) or not np.isfinite(v2):
            malo = "   <-- NO FINITO"
            problemas.append((entry["name"], ndim, "NaN/inf"))

        print(f"{entry['name']:<12}{ndim:>4}  {v1:>26.17e}  {v2:>26.17e}{malo}")
        h.update(f"{entry['name']}|{ndim}|{v1!r}|{v2!r}".encode())

print("-" * 78)
print(f"HUELLA GLOBAL sha256: {h.hexdigest()}")
print()

# ─────────────────────────────────────────────────────────────────────────────
# BARRIDO DE NaN — la parte decisiva.
# default_rng tiene stream estable entre versiones de numpy, asi que los 2000
# puntos son LOS MISMOS en cualquier maquina. Si una funcion devuelve NaN aqui
# y no en la otra maquina, su implementacion cambio.
# ─────────────────────────────────────────────────────────────────────────────
print("=" * 78)
print("BARRIDO DE NaN — 2000 puntos deterministas por funcion")
print("=" * 78)
rotos = []
for entry in CEC_REGISTRY:
    for ndim in (10, 20 if 50 not in entry["dims"] else 50):
        rng = np.random.default_rng(12345)
        P = rng.uniform(-100, 100, size=(2000, ndim))
        try:
            obj = entry["class"](ndim=ndim)
            g = obj.f_global
            v = np.array([obj.evaluate(x) - g for x in P])
        except Exception as e:
            print(f"  {entry['name']:<12} dim={ndim:<3} EXCEPCION {type(e).__name__}: {e}")
            rotos.append((entry["name"], ndim, "excepcion"))
            continue
        nf = int(np.sum(~np.isfinite(v)))
        if nf:
            print(f"  {entry['name']:<12} dim={ndim:<3} {nf:>5}/2000 NO FINITOS  <-- ROTA")
            rotos.append((entry["name"], ndim, f"{nf}/2000 no finitos"))

if rotos:
    print(f"\nRESUMEN: {len(rotos)} funcion(es) rota(s) en esta maquina:")
    for n, d, q in rotos:
        print(f"  {n} dim={d}: {q}")
    print("\nUna funcion rota aqui invalida TODO resultado calculado sobre ella")
    print("en esta maquina, para cualquier algoritmo.")
else:
    print("  Ninguna. Las 71 funciones son finitas en los 2000 puntos.")
print()
if problemas:
    print(f"FUNCIONES CON PROBLEMA ({len(problemas)}):")
    for n, d, q in problemas:
        print(f"  {n} dim={d}: {q}")
else:
    print("Todas las funciones evaluaron a valores finitos.")
