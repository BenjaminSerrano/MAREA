"""
¿Es ESTE intérprete el mismo que produjo los resultados guardados?

Re-ejecuta una muestra de corridas de los COMPETIDORES (mismas semillas,
mismo presupuesto) y las compara bit a bit contra lo que hay en los JSON.

  - 100% idénticas  -> entorno correcto, se puede seguir usando.
  - algo distinto   -> entorno distinto: todo lo que se calcule aquí NO es
                       comparable con lo ya guardado.

Uso:
    python3 verificar_entorno.py
    /ruta/a/otro/python3 verificar_entorno.py     # para probar candidatos
"""
import sys, os, json
import warnings
warnings.filterwarnings("ignore")
import numpy as np

sys.path.insert(0, ".")
from experimento_config import CEC_REGISTRY, POP_SIZE, result_path

MAXFES = 5000
NDIM = 10
T = MAXFES // POP_SIZE
PIDS = [0, 40]
RUNS = [0, 1]

print("=" * 74)
print("VERIFICACION DE ENTORNO")
print("=" * 74)
print(f"python : {sys.executable}")
print(f"version: {sys.version.split()[0]}")
print(f"numpy  : {np.__version__}")
try:
    import opfunu
    print(f"opfunu : {opfunu.__version__}")
except Exception as e:
    print(f"opfunu : ??? ({e})")
print("=" * 74)

import competidores as comp

total = ok = 0
por_alg = {}

for alg in ["PSO", "DE", "CLONALG", "GWO", "WOA", "FA", "ABC"]:
    path = result_path(alg, NDIM, MAXFES)
    if not os.path.exists(path):
        print(f"  {alg:<9} (sin archivo guardado — se omite)")
        continue
    guardado = json.load(open(path))
    fn = getattr(comp, alg.lower())
    a_ok = a_tot = 0
    peor = 0.0
    for pid in PIDS:
        entry = CEC_REGISTRY[pid]
        obj = entry["class"](ndim=NDIM)
        g = obj.f_global
        def f(x, obj=obj, g=g):
            return obj.evaluate(x) - g
        bounds = [(-100.0, 100.0)] * NDIM
        for run in RUNS:
            clave = f"pid{pid}_run{run}"
            if clave not in guardado:
                continue
            ref = guardado[clave]
            got = float(fn(f, bounds, N=POP_SIZE, T=T, seed=run)["best_f"])
            a_tot += 1
            if got == ref:
                a_ok += 1
            else:
                rel = abs(got - ref) / abs(ref) if ref else float("inf")
                peor = max(peor, rel)
    if a_tot:
        estado = "IDENTICO" if a_ok == a_tot else f"DIFIERE (peor rel={peor:.2e})"
        print(f"  {alg:<9} {a_ok:>3}/{a_tot} identicas   {estado}")
        por_alg[alg] = (a_ok, a_tot)
        ok += a_ok
        total += a_tot

print("-" * 74)
if total == 0:
    print("No habia resultados guardados para comparar.")
elif ok == total:
    print(f"VEREDICTO: {ok}/{total} identicas -> ESTE ES EL ENTORNO CORRECTO.")
    print("Se puede seguir calculando aqui; los resultados son comparables.")
else:
    print(f"VEREDICTO: solo {ok}/{total} identicas -> ENTORNO DISTINTO.")
    print("Lo que se calcule aqui NO es comparable con lo ya guardado.")
    print("Buscar el interprete que usaba el servidor Jupyter y usar ese.")
