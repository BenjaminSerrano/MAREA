"""
Ejecuta MAREA con el interprete que lo invoca — sin kernels de por medio.

ANTES de calcular nada, verifica que el entorno no sea uno de los rotos
(opfunu 1.0.1 devuelve no-finito en F17-2017 y F29-2017 a dim=10, al 100%).
Si lo esta, aborta sin escribir un solo archivo.

Uso:
    /ruta/al/python3 correr_marea.py --check
    /ruta/al/python3 correr_marea.py --maxfes 5000
    /ruta/al/python3 correr_marea.py --maxfes 50000 --workers 20
"""
import argparse, sys, os
import warnings
warnings.filterwarnings("ignore")


def revisar_entorno():
    """Devuelve True si el entorno puede evaluar las 71 funciones."""
    import numpy as np
    print("=" * 70)
    print(f"interprete : {sys.executable}")
    print(f"python     : {sys.version.split()[0]}")
    print(f"numpy      : {np.__version__}")
    try:
        import opfunu
        print(f"opfunu     : {opfunu.__version__}")
    except Exception as e:
        print(f"opfunu     : ??? ({e})")
    print(f"directorio : {os.getcwd()}")
    print("=" * 70)

    from experimento_config import CEC_REGISTRY, result_path
    print(f"escribe en : {os.path.abspath(os.path.dirname(result_path('MAREA', 10, 5000)))}")

    rng = np.random.default_rng(12345)
    rotas = []
    for entry in CEC_REGISTRY:
        for ndim in (10, 20 if 50 not in entry["dims"] else 50):
            P = rng.uniform(-100, 100, size=(50, ndim))
            try:
                obj = entry["class"](ndim=ndim)
                g = obj.f_global
                v = np.array([obj.evaluate(x) - g for x in P])
                if not np.all(np.isfinite(v)):
                    rotas.append((entry["name"], ndim))
            except Exception as exc:
                rotas.append((entry["name"], f"{ndim} ({type(exc).__name__})"))

    if rotas:
        print()
        print("ENTORNO ROTO — funciones que no evaluan a valores finitos:")
        for n, d in rotas:
            print(f"  {n} dim={d}")
        print()
        print("No se calcula nada. Los resultados serian NaN y no comparables")
        print("con lo ya guardado. Usar el interprete que reproduce los datos")
        print("existentes (ver verificar_entorno.py).")
        return False

    print("entorno OK : las 71 funciones evaluan finitas")
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--maxfes", type=int, help="presupuesto (5000, 50000, 500000...)")
    ap.add_argument("--dims", type=int, nargs="+", default=[10, 50])
    ap.add_argument("--workers", type=int,
                    default=max(1, (os.cpu_count() or 4) - 4))
    ap.add_argument("--check", action="store_true",
                    help="solo verificar el entorno y salir")
    ap.add_argument("--algo", default="MAREA",
                    choices=["MAREA", "MAREA_REL", "MAREABIO", "MAREABIO_REL",
                             "PSO", "DE", "ABC", "GWO", "WOA", "FA", "CLONALG"],
                    help="algoritmo a correr (default: MAREA)")
    ap.add_argument("--outdir", default=None,
                    help="directorio de resultados (default: resultados_experimento)")
    args = ap.parse_args()

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

    if args.outdir:
        import experimento_config as _cfg
        os.makedirs(args.outdir, exist_ok=True)
        _cfg.RESULTS_DIR = args.outdir

    if not revisar_entorno():
        sys.exit(1)
    if args.check:
        sys.exit(0)
    if not args.maxfes:
        ap.error("hace falta --maxfes (o usar --check)")

    from experimento_config import run_experiment, _import_alg

    print(f"algoritmo  : {args.algo}")
    if args.algo.startswith("MAREA"):
        if args.algo.startswith("MAREABIO"):
            from marea_bio_benchmark import MAREABioParams as _P
        else:
            from marea_benchmark import MAREAParams as _P
        p = _P()
        if p.early_stop:
            print(f"\nERROR: {_P.__name__}.early_stop es True. El benchmark exige")
            print("agotar MaxFES. Falta copiar el modulo actualizado.")
            sys.exit(1)
        print(f"early_stop : {p.early_stop} (agota MaxFES)")
    else:
        # FA y CLONALG cortan solos al agotar (T+1)*N evaluaciones; el resto
        # gasta N por iteracion por construccion. Ver competidores._presupuesto.
        import competidores as _c
        if args.algo in ("FA", "CLONALG") and not hasattr(_c, "_presupuesto"):
            print(f"\nERROR: competidores.py no tiene el control de presupuesto.")
            print(f"{args.algo} gastaria mucho mas que MaxFES. Falta copiar el modulo.")
            sys.exit(1)
    print(f"workers    : {args.workers}\n")

    # El worker reimporta por nombre; esta alg_fn es solo para la firma.
    alg_fn = _import_alg(args.algo)

    for ndim in args.dims:
        run_experiment(args.algo, alg_fn, ndim=ndim,
                       maxfes=args.maxfes, n_workers=args.workers)


if __name__ == "__main__":
    main()
