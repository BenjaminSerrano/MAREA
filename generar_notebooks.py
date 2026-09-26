"""
Genera un notebook Jupyter (.ipynb) por algoritmo.
Ejecutar: python3 generar_notebooks.py
"""

import json, os

# ══════════════════════════════════════════════════════════════════════════════
# DEFINICIÓN DE ALGORITMOS
# ══════════════════════════════════════════════════════════════════════════════

MAXFES = [5_000, 50_000, 500_000, 5_000_000]
DIMS   = [10, 50]

ALGORITMOS = [
    {
        "name":     "MAREA",
        "import":   "from marea_benchmark import marea, MAREAParams\n"
                    "p = MAREAParams()\n"
                    "def alg_fn(f, bounds, N, T, seed):\n"
                    "    p.N = N; p.T = T\n"
                    "    return marea(f, bounds, p, seed=seed)",
        "desc":     "MAREA-base — Metaheuristic Algorithm for Regulated Exploration and Adaptation\n"
                    "Parámetros calibrados: N=80, α_max=0.99, β_max=0.3, γ_max=0.2",
        "color":    "#1f77b4",
    },
    {
        "name":     "PSO",
        "import":   "from competidores import pso\n"
                    "def alg_fn(f, bounds, N, T, seed):\n"
                    "    return pso(f, bounds, N=N, T=T, seed=seed)",
        "desc":     "PSO — Particle Swarm Optimization (Kennedy & Eberhart, 1995)\n"
                    "Parámetros: w=0.729, c1=1.494, c2=1.494 (Clerc & Kennedy, 2002)",
        "color":    "#ff7f0e",
    },
    {
        "name":     "DE",
        "import":   "from competidores import de\n"
                    "def alg_fn(f, bounds, N, T, seed):\n"
                    "    return de(f, bounds, N=N, T=T, seed=seed)",
        "desc":     "DE/rand/1/bin — Differential Evolution (Storn & Price, 1997)\n"
                    "Parámetros: F=0.8, CR=0.9",
        "color":    "#2ca02c",
    },
    {
        "name":     "CLONALG",
        "import":   "from competidores import clonalg\n"
                    "def alg_fn(f, bounds, N, T, seed):\n"
                    "    return clonalg(f, bounds, N=N, T=T, seed=seed)",
        "desc":     "CLONALG — Artificial Immune System (de Castro & Von Zuben, 2002)\n"
                    "Parámetros: n_clones=10, beta=100, d=0.1",
        "color":    "#d62728",
    },
    {
        "name":     "GWO",
        "import":   "from competidores import gwo\n"
                    "def alg_fn(f, bounds, N, T, seed):\n"
                    "    return gwo(f, bounds, N=N, T=T, seed=seed)",
        "desc":     "GWO — Grey Wolf Optimizer (Mirjalili et al., 2014)\n"
                    "Sin parámetros libres (a decrece de 2 a 0)",
        "color":    "#9467bd",
    },
    {
        "name":     "WOA",
        "import":   "from competidores import woa\n"
                    "def alg_fn(f, bounds, N, T, seed):\n"
                    "    return woa(f, bounds, N=N, T=T, seed=seed)",
        "desc":     "WOA — Whale Optimization Algorithm (Mirjalili & Lewis, 2016)\n"
                    "Parámetros: b=1.0",
        "color":    "#8c564b",
    },
    {
        "name":     "ABC",
        "import":   "from competidores import abc\n"
                    "def alg_fn(f, bounds, N, T, seed):\n"
                    "    return abc(f, bounds, N=N, T=T, seed=seed)",
        "desc":     "ABC — Artificial Bee Colony (Karaboga & Basturk, 2007)\n"
                    "Parámetros: limit=N*D/2 (recomendado por Karaboga)",
        "color":    "#e377c2",
    },
    {
        "name":     "FA",
        "import":   "from competidores import fa\n"
                    "def alg_fn(f, bounds, N, T, seed):\n"
                    "    return fa(f, bounds, N=N, T=T, seed=seed)",
        "desc":     "FA — Firefly Algorithm (Yang, 2010)\n"
                    "Parámetros: alpha=0.5, beta0=1.0, gamma=1.0",
        "color":    "#7f7f7f",
    },
]


# ══════════════════════════════════════════════════════════════════════════════
# PLANTILLA DE CELDAS
# ══════════════════════════════════════════════════════════════════════════════

def make_notebook(alg):
    name  = alg["name"]
    color = alg["color"]

    cells = []

    # ── Celda 1: Título y descripción ─────────────────────────────────────────
    cells.append({
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            f"# Experimento MAREA — Competidor: **{name}**\n",
            "\n",
            f"**{alg['desc']}**\n",
            "\n",
            "### Protocolo experimental\n",
            "- **Benchmarks:** 71 funciones CEC 2014 (30) + CEC 2017 (29) + CEC 2022 (12)\n",
            "- **Dimensiones:** 10D y 50D (CEC2022: 10D y 20D)\n",
            "- **MaxFES:** 5,000 / 50,000 / 500,000 / 5,000,000\n",
            "- **Corridas:** 51 independientes por configuración\n",
            "- **Guardado:** incremental en JSON — se puede interrumpir y reanudar\n",
            "\n",
            "> Ejecutar las celdas en orden. Cada celda de experimento es independiente\n",
            "> y guarda sus resultados automáticamente al terminar.\n",
        ]
    })

    # ── Celda 2: Imports y setup ───────────────────────────────────────────────
    cells.append({
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "import sys, os\n",
            "import numpy as np\n",
            "import matplotlib.pyplot as plt\n",
            "import warnings\n",
            "warnings.filterwarnings('ignore')\n",
            "\n",
            "# Asegurar que el directorio del proyecto esté en el path\n",
            "PROJECT_DIR = os.path.dirname(os.path.abspath('__file__'))\n",
            "if PROJECT_DIR not in sys.path:\n",
            "    sys.path.insert(0, PROJECT_DIR)\n",
            "\n",
            "from experimento_config import (\n",
            "    run_experiment, load_partial, CEC_REGISTRY,\n",
            "    MAXFES, DIMS, N_RUNS, POP_SIZE, RESULTS_DIR\n",
            ")\n",
            "\n",
            f"ALG_NAME = '{name}'\n",
            "\n",
            "# Cargar el algoritmo\n",
            f"{alg['import']}\n",
            "\n",
            f"print(f'✓ {{ALG_NAME}} cargado correctamente')\n",
            f"print(f'  Resultados en: {{RESULTS_DIR}}/')\n",
            f"print(f'  Funciones CEC disponibles: {{len(CEC_REGISTRY)}}')\n",
        ]
    })

    # ── Celda 3: Verificación rápida (smoke test) ─────────────────────────────
    cells.append({
        "cell_type": "markdown",
        "metadata": {},
        "source": ["## Verificación rápida (smoke test)\n",
                   "Corre 1 función × 3 semillas para confirmar que todo funciona antes del experimento completo."]
    })
    cells.append({
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "from experimento_config import get_problems_for_dim\n",
            "import time\n",
            "\n",
            "problems_10d = get_problems_for_dim(10)\n",
            "test_prob    = problems_10d[0]   # F1 de CEC2014\n",
            "print(f'Probando con: {test_prob[\"name\"]} | dim={test_prob[\"dim\"]}')\n",
            "\n",
            "t0 = time.time()\n",
            "for seed in range(3):\n",
            "    r = alg_fn(test_prob['f'], test_prob['bounds'], N=POP_SIZE, T=62, seed=seed)\n",
            "    print(f'  seed={seed}  best_f={r[\"best_f\"]:.4e}')\n",
            "print(f'Tiempo: {time.time()-t0:.2f}s  ✓ Listo para experimento completo')\n",
        ]
    })

    # ── Celdas 4–7: Una por MaxFES ────────────────────────────────────────────
    maxfes_labels = {
        5_000:       ("MaxFES = 5,000",      "Caso de uso: optimización costosa / tiempo-real"),
        50_000:      ("MaxFES = 50,000",     "Caso de uso: aplicaciones prácticas estándar"),
        500_000:     ("MaxFES = 500,000",    "Caso de uso: convención CEC estándar"),
        5_000_000:   ("MaxFES = 5,000,000",  "Caso de uso: búsqueda exploratoria profunda"),
    }

    for mf in MAXFES:
        title, uso = maxfes_labels[mf]
        cells.append({
            "cell_type": "markdown",
            "metadata": {},
            "source": [f"## Experimento — {title}\n", f"*{uso}*\n"]
        })
        cells.append({
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                f"# ── {title} ──────────────────────────────────────────────\n",
                f"# Se corre sobre 10D y 50D. Cada corrida guarda su resultado.\n",
                f"# Si se interrumpe, la siguiente ejecución retoma desde donde quedó.\n",
                "\n",
                f"for ndim in DIMS:\n",
                f"    run_experiment(ALG_NAME, alg_fn, ndim=ndim, maxfes={mf})\n",
            ]
        })

    # ── Celda final: Resumen de resultados ────────────────────────────────────
    cells.append({
        "cell_type": "markdown",
        "metadata": {},
        "source": ["## Resumen de resultados\n",
                   "Muestra el ranking promedio de este algoritmo en cada configuración."]
    })
    cells.append({
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "import json\n",
            "\n",
            "print(f'\\n=== RESUMEN {ALG_NAME} ===')\n",
            "print(f'{\"Config\":<30} {\"Completadas\":>12} {\"Mean\":>14} {\"Best\":>14}')\n",
            "print('-'*72)\n",
            "\n",
            "for mf in MAXFES:\n",
            "    for nd in DIMS:\n",
            "        data = load_partial(ALG_NAME, nd, mf)\n",
            "        vals = list(data.values())\n",
            "        if vals:\n",
            "            label = f'D={nd:2d} MaxFES={mf:>10,}'\n",
            "            print(f'{label:<30} {len(vals):>12} {np.mean(vals):>14.4e} {np.min(vals):>14.4e}')\n",
            "        else:\n",
            "            print(f'D={nd:2d} MaxFES={mf:>10,}  (sin resultados aún)')\n",
        ]
    })

    # ── Estructura .ipynb ─────────────────────────────────────────────────────
    return {
        "nbformat": 4,
        "nbformat_minor": 5,
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3"
            },
            "language_info": {
                "name": "python",
                "version": "3.11.0"
            }
        },
        "cells": cells
    }


# ══════════════════════════════════════════════════════════════════════════════
# GENERAR TODOS LOS NOTEBOOKS
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    output_dir = os.path.dirname(os.path.abspath(__file__))
    print("Generando notebooks...")

    for alg in ALGORITMOS:
        nb   = make_notebook(alg)
        path = os.path.join(output_dir, f"experimento_{alg['name']}.ipynb")
        with open(path, "w") as fh:
            json.dump(nb, fh, indent=1)
        print(f"  ✓ {path}")

    print(f"\n{len(ALGORITMOS)} notebooks generados.")
    print("Para abrir en Jupyter: jupyter notebook")
