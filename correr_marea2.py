"""
Ejecuta MAREA-2 con el interprete que lo invoca — sin kernels de por medio.

Mismo contrato que correr_marea.py: ANTES de calcular nada verifica que el
entorno no sea uno de los rotos (opfunu 1.0.1 devuelve no-finito en F17-2017 y
F29-2017 a dim=10, al 100%). Si lo esta, aborta sin escribir un solo archivo.

Añade una comprobacion propia de MAREA-2: el algoritmo controla el presupuesto
en EVALUACIONES, asi que se hace una corrida de prueba y se exige fes <= MaxFES
antes de gastar horas. Es el equivalente al chequeo de early_stop de MAREA.

Uso:
    /ruta/al/python3 correr_marea2.py --check
    /ruta/al/python3 correr_marea2.py --ablacion --maxfes 5000  --workers 20
    /ruta/al/python3 correr_marea2.py --ablacion --maxfes 50000 --workers 20
    /ruta/al/python3 correr_marea2.py --maxfes 5000  --workers 20
    /ruta/al/python3 correr_marea2.py --maxfes 50000 --config sin_memoria

El protocolo (sin --ablacion) escribe MAREA2_D{dim}_FES{maxfes}.json con la
config FULL, o MAREA2_<variante>_D... si se pasa --config. La ablacion escribe
en resultados_experimento/ablacion/ y se puede interrumpir y reanudar.
"""
import argparse, sys, os
import warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from correr_marea import revisar_entorno    # noqa: E402  (mismo chequeo de opfunu)


def revisar_presupuesto(cfg):
    """Una corrida barata: MAREA-2 no debe pasarse de MaxFES nunca."""
    from experimento_config import get_problems_for_dim, POP_SIZE
    from marea2 import marea2, config

    T = 62
    maxfes = POP_SIZE * T
    prob = get_problems_for_dim(10)[0]
    p = config(cfg); p.N0 = POP_SIZE
    r = marea2(prob["f"], prob["bounds"], p, maxfes=maxfes, seed=0)

    if r["fes"] > maxfes:
        print(f"\nERROR: {cfg} gasto {r['fes']} FES con MaxFES={maxfes}.")
        print("El benchmark exige no pasarse. Falta copiar marea2.py actualizado.")
        return False
    print(f"presupuesto: {r['fes']}/{maxfes} FES en {r['iters']} iteraciones (OK)")
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--maxfes", type=int, help="presupuesto (5000, 50000, 500000...)")
    ap.add_argument("--dims", type=int, nargs="+", default=[10, 50])
    ap.add_argument("--workers", type=int,
                    default=max(1, (os.cpu_count() or 4) - 4))
    ap.add_argument("--check", action="store_true",
                    help="solo verificar el entorno y salir")
    ap.add_argument("--config", default="FULL",
                    help="configuracion de marea2.config() (default: FULL)")
    ap.add_argument("--ablacion", action="store_true",
                    help="Fase A: corre las 14 variantes en vez del protocolo")
    ap.add_argument("--cfgs", nargs="*", default=None,
                    help="solo con --ablacion: subconjunto de variantes")
    ap.add_argument("--runs", type=int, default=15,
                    help="solo con --ablacion: corridas por funcion (default: 15)")
    ap.add_argument("--outdir", default=None,
                    help="directorio de resultados (default: resultados_experimento)")
    args = ap.parse_args()

    if args.outdir:
        import experimento_config as _cfg
        os.makedirs(args.outdir, exist_ok=True)
        _cfg.RESULTS_DIR = args.outdir

    from marea2 import ABLATION_CONFIGS, config
    cfgs = args.cfgs or ABLATION_CONFIGS
    for c in ([args.config] if not args.ablacion else cfgs):
        try:
            config(c)
        except ValueError as e:
            ap.error(str(e))

    if not revisar_entorno():
        sys.exit(1)
    if not revisar_presupuesto(args.config if not args.ablacion else "FULL"):
        sys.exit(1)
    if args.check:
        sys.exit(0)
    if not args.maxfes:
        ap.error("hace falta --maxfes (o usar --check)")

    print(f"workers    : {args.workers}")

    if args.ablacion:
        import ablacion_marea2 as abl
        print(f"modo       : ablacion, {len(cfgs)} variantes, {args.runs} corridas")
        print(f"escribe en : {os.path.abspath(abl.ABL_DIR)}\n")
        abl.run(cfgs, dims=args.dims, maxfes=args.maxfes,
                runs=args.runs, workers=args.workers)
        return

    from experimento_config import run_experiment, _import_alg
    alg = "MAREA2" if args.config == "FULL" else f"MAREA2_{args.config}"
    print(f"modo       : protocolo completo, config={args.config} -> {alg}\n")

    # El worker reimporta por nombre; esta alg_fn es solo para la firma.
    alg_fn = _import_alg(alg)

    for ndim in args.dims:
        run_experiment(alg, alg_fn, ndim=ndim,
                       maxfes=args.maxfes, n_workers=args.workers)


if __name__ == "__main__":
    main()
