"""
Pipeline de análisis del paper OE1 (Fase 2 del plan, conforme al setup §5–§9).

    python3 analisis_paper.py              # todas las tablas y figuras en paper/
    python3 analisis_paper.py --regresion  # control: reproduce la tabla 2.1 de la memoria

Lee la lista de algoritmos de `analisis_config.json` y los JSON de
`resultados_experimento/`. Métrica principal: media de las 51 corridas;
mediana como robustez. Errores < 1e-8 cuentan como 0 (regla de CEC 2017).

Salida:
  paper/tablas/        CSV + LaTeX (booktabs) de cada análisis
  paper/suplementario/ media, mediana, SD, mejor y peor por función
  paper/figuras/       rangos por MaxFES, diagramas de diferencia crítica (Nemenyi), W/T/L
"""
import argparse, csv, json, os, sys
from itertools import combinations
from math import comb
import numpy as np
from scipy.stats import (rankdata, friedmanchisquare, ranksums, wilcoxon,
                         norm, f as fdist, studentized_range)

D_RES = "resultados_experimento"
OUT = "paper"
CFG = json.load(open("analisis_config.json"))
ALPHA = CFG["alpha"]
CELDAS = [(d, f) for d in CFG["dims"] for f in CFG["maxfes"]]
BLAB = {5000: "5k", 50000: "50k", 500000: "500k"}
DLAB = {10: "D10", 50: "D50"}            # D50 lleva CEC2022 en D20 (setup §1)

# 71 problemas en el orden de experimento_config.CEC_REGISTRY: 30 CEC2014 + 29 CEC2017 + 12 CEC2022.
# opfunu omite la F2 oficial de CEC2017: su F_k es la F oficial k (k=1) o k+1 (k>=2).
SUITES = {"CEC2014": range(0, 30), "CEC2017": range(30, 59), "CEC2022": range(59, 71)}

def funcion(pid):
    if pid < 30: return "CEC2014", pid + 1
    if pid < 59: k = pid - 29; return "CEC2017", k if k == 1 else k + 1
    return "CEC2022", pid - 58

def clase(pid):
    s, n = funcion(pid)          # numeración oficial
    lim = {"CEC2014": (3, 16, 22), "CEC2017": (3, 10, 20), "CEC2022": (1, 5, 8)}[s]
    return ("unimodal" if n <= lim[0] else "multimodal" if n <= lim[1]
            else "hibrida" if n <= lim[2] else "composicion")

CLASES = ["unimodal", "multimodal", "hibrida", "composicion"]

def et(a): return CFG["etiquetas"].get(a, a)
def celda(d, f): return f"{DLAB[d]} {BLAB[f]}"

# ── carga ───────────────────────────────────────────────────────────────────
_cache = {}
def load(a, d, f, cero=True, abl=False):
    """abl=True: brazo de la ablación (resultados_experimento/ablacion/MAREA2_<a>, 15 corridas)."""
    key = (a, d, f, cero, abl)
    if key not in _cache:
        ruta = f"{D_RES}/ablacion/MAREA2_{a}" if abl else f"{D_RES}/{a}"
        M = np.full((71, CFG["corridas_ablacion"] if abl else CFG["corridas"]), np.nan)
        for k, v in json.load(open(f"{ruta}_D{d}_FES{f}.json")).items():   # claves pid/fun{i}_run{r}
            M[int(k.split("_")[0][3:]), int(k.split("run")[1])] = v
        if np.isnan(M).any() or not np.isfinite(M).all():
            sys.exit(f"{a} D{d} {f}: faltan corridas o hay no finitos")
        if cero: M[M < CFG["umbral_cero"]] = 0.0
        _cache[key] = M
    return _cache[key]

# ── Friedman, Shaffer y Nemenyi ────────────────────────────────────────────
def _S(k, memo={}):
    """Números posibles de hipótesis verdaderas entre k algoritmos (Shaffer 1986)."""
    if k <= 1: return {0}
    if k not in memo:
        memo[k] = {comb(j, 2) + x for j in range(1, k + 1) for x in _S(k - j)}
    return memo[k]

def shaffer(p):
    """p ajustados por el procedimiento estático de Shaffer (García y Herrera 2008)."""
    p = np.asarray(p); m = len(p); k = int(round((1 + np.sqrt(1 + 8 * m)) / 2))
    S = sorted(_S(k)); orden = np.argsort(p); adj = np.empty(m); run = 0.0
    for i, j in enumerate(orden):
        t = max(s for s in S if s <= m - i)
        run = max(run, min(1.0, t * p[j])); adj[j] = run
    return adj

def ranking(algs, d, f, idx=range(71), stat=np.mean, cero=True):
    """Rangos por problema (media de 51 corridas) y todas las pruebas sobre ellos."""
    idx = list(idx)
    X = np.array([stat(load(a, d, f, cero), axis=1)[idx] for a in algs])   # k × N
    R = np.apply_along_axis(rankdata, 0, X)
    k, N = R.shape; rm = R.mean(1)
    chi2, p_f = friedmanchisquare(*X)
    F_id = (N - 1) * chi2 / (N * (k - 1) - chi2)
    p_id = fdist.sf(F_id, k - 1, (k - 1) * (N - 1))
    se = np.sqrt(k * (k + 1) / (6 * N))
    cd = studentized_range.ppf(1 - ALPHA, k, np.inf) / np.sqrt(2) * se
    pares = list(combinations(range(k), 2))
    diff = np.array([abs(rm[i] - rm[j]) for i, j in pares])
    p_z = 2 * norm.sf(diff / se)
    p_nem = studentized_range.sf(diff / se * np.sqrt(2), k, np.inf)
    p_sh = shaffer(p_z)
    pos = np.empty(k, int); pos[np.argsort(rm, kind="stable")] = np.arange(1, k + 1)
    return dict(algs=algs, X=X, R=R, rm=rm, sd=R.std(1, ddof=1), pos=pos, N=N, k=k,
                chi2=chi2, p_f=p_f, F_id=F_id, p_id=p_id, cd=cd, pares=pares,
                diff=diff, p_z=p_z, p_sh=p_sh, p_nem=np.minimum(p_nem, 1.0))

def p_par(r, i, j, cual):
    for n, (a, b) in enumerate(r["pares"]):
        if {a, b} == {i, j}: return r[cual][n]

def pct_distintos(r, cual):
    k = r["k"]
    return np.array([100 * sum(p_par(r, i, j, cual) < ALPHA for j in range(k) if j != i) / (k - 1)
                     for i in range(k)])

# ── W/T/L y signed-rank ────────────────────────────────────────────────────
def wtl(A, B, idx=range(71)):
    """Rank-sum por función (51 vs 51), α=0,05; el sentido lo decide la media."""
    w = t = l = 0
    for i in idx:
        if np.array_equal(A[i], B[i]) or ranksums(A[i], B[i]).pvalue >= ALPHA: t += 1
        elif A[i].mean() < B[i].mean(): w += 1
        else: l += 1
    return w, t, l

def signed_rank(a, b):
    """Wilcoxon signed-rank pareado por problema sobre las medias (Derrac et al. 2011).
    R+ = suma de rangos donde gana la propuesta (a < b); los empates se reparten."""
    dlt = b - a
    r = rankdata(np.abs(dlt))
    rp = r[dlt > 0].sum() + r[dlt == 0].sum() / 2
    rn = r[dlt < 0].sum() + r[dlt == 0].sum() / 2
    p = 1.0 if np.all(dlt == 0) else wilcoxon(a, b, zero_method="zsplit").pvalue
    return rp, rn, p

def holm(p):
    p = np.asarray(p); o = np.argsort(p); m = len(p); adj = np.empty(m); run = 0.0
    for i, j in enumerate(o):
        run = max(run, min(1.0, (m - i) * p[j])); adj[j] = run
    return adj

# ── salida CSV + LaTeX ─────────────────────────────────────────────────────
def escribir(nombre, cab, filas, titulo, sub="tablas", tex=True, nota=None):
    os.makedirs(f"{OUT}/{sub}", exist_ok=True)
    with open(f"{OUT}/{sub}/{nombre}.csv", "w", newline="") as fh:
        csv.writer(fh, lineterminator="\n").writerows([cab] + filas)
    if not tex: return
    esc = lambda s: str(s).replace("_", r"\_").replace("%", r"\%").replace("&", r"\&")
    env = "longtable" if len(filas) > 40 else "tabular"
    L = [f"% {titulo}" + (f"\n% {nota}" if nota else "")]
    if env == "tabular": L += [r"\begin{table}[ht]\centering\small", f"\\caption{{{esc(titulo)}}}"]
    L += [f"\\begin{{{env}}}{{l{'r' * (len(cab) - 1)}}}"]
    if env == "longtable": L += [f"\\caption{{{esc(titulo)}}}\\\\"]
    L += [r"\toprule", " & ".join(esc(c) for c in cab) + r" \\", r"\midrule"]
    if env == "longtable": L += [r"\endhead"]
    L += [" & ".join(esc(x) for x in fl) + r" \\" for fl in filas]
    L += [r"\bottomrule", f"\\end{{{env}}}"]
    if nota: L += [f"\\par\\footnotesize {esc(nota)}"]
    if env == "tabular": L += [r"\end{table}"]
    open(f"{OUT}/{sub}/{nombre}.tex", "w").write("\n".join(L) + "\n")

def fp(p): return "<1e-15" if p < 1e-15 else f"{p:.2e}" if p < 1e-3 else f"{p:.3f}"

# ── control de regresión ───────────────────────────────────────────────────
def regresion():
    rg = CFG["regresion"]; algs = rg["algoritmos"]; mal = 0
    for c, (d, f) in enumerate(CELDAS):
        r = ranking(algs, d, f, stat=np.median, cero=False)
        for i, a in enumerate(algs):
            esp = rg["esperado"][a][c]
            if abs(round(r["rm"][i], 2) - esp) > 0.005:
                mal += 1; print(f"  DISTINTO {a} {celda(d, f)}: {r['rm'][i]:.3f} vs {esp}")
    print(f"Regresión (tabla 2.1, {len(algs)} algoritmos × 6 celdas): "
          + ("OK, 72/72 iguales" if not mal else f"{mal} diferencias"))
    return mal == 0

# ── análisis principal ─────────────────────────────────────────────────────
def main():
    algs = CFG["algoritmos"]; P = CFG["propuesta"]; ip = algs.index(P)
    comp = [a for a in algs if a != P]
    AGG = {(d, f): ranking(algs, d, f) for d, f in CELDAS}
    SUB = {(s, d, f): ranking(algs, d, f, ix) for s, ix in SUITES.items() for d, f in CELDAS}
    MED = {(d, f): ranking(algs, d, f, stat=np.median) for d, f in CELDAS}
    cab6 = ["Algoritmo"] + [celda(d, f) for d, f in CELDAS]
    nota_d50 = "D50: CEC2014 y CEC2017 en D50, CEC2022 en D20."

    # T1 · rangos agregados (CEC_71) con puesto
    orden = np.argsort(np.mean([AGG[c]["rm"] for c in CELDAS], axis=0))
    filas = [[et(algs[i])] + [f"{AGG[c]['rm'][i]:.2f} ({AGG[c]['pos'][i]})" for c in CELDAS]
             for i in orden]
    filas.append(["Friedman p"] + [fp(AGG[c]["p_f"]) for c in CELDAS])
    filas.append(["CD Nemenyi"] + [f"{AGG[c]['cd']:.2f}" for c in CELDAS])
    escribir("t1_rangos_agregados", cab6, filas,
             "Rango medio de Friedman (media de 51 corridas), 71 funciones", nota=nota_d50)

    # T2 · rangos por benchmark (18 rankings)
    cab = ["Benchmark", "Algoritmo"] + [celda(d, f) for d, f in CELDAS]
    filas = []
    for s in SUITES:
        for i in orden:
            filas.append([s, et(algs[i])] + [f"{SUB[(s,) + c]['rm'][i]:.2f} ({SUB[(s,) + c]['pos'][i]})"
                                             for c in CELDAS])
        filas.append([s, "Friedman p"] + [fp(SUB[(s,) + c]["p_f"]) for c in CELDAS])
    escribir("t2_rangos_por_benchmark", cab, filas,
             "Rango medio de Friedman por benchmark", nota=nota_d50 + " CEC2022 tiene 12 funciones.")

    # T3 · resumen de posiciones y robustez (§9): 6 rankings agregados y 18 por benchmark
    cab = ["Algoritmo"] + [celda(d, f) for d, f in CELDAS] + [
        "Pos. media (6)", "SD pos. (6)", "Pos. media (18)", "SD pos. (18)", "SD rangos"]
    filas = []
    for i in orden:
        p6 = [AGG[c]["pos"][i] for c in CELDAS]
        p18 = [SUB[k]["pos"][i] for k in SUB]
        sdr = np.mean([AGG[c]["sd"][i] for c in CELDAS])
        filas.append([et(algs[i])] + p6 + [f"{np.mean(p6):.2f}", f"{np.std(p6, ddof=1):.2f}",
                                           f"{np.mean(p18):.2f}", f"{np.std(p18, ddof=1):.2f}",
                                           f"{sdr:.2f}"])
    escribir("t3_resumen_posiciones", cab, filas, "Puestos por celda y robustez de la posición",
             nota="SD rangos: desviación estándar de los rangos por función, promediada en las 6 celdas.")

    # T4 · Friedman / Iman-Davenport / CD por ranking (24)
    cab = ["Ranking", "N", "k", "chi2", "p Friedman", "F Iman-Davenport", "p ID", "CD Nemenyi"]
    filas = [["CEC_71 " + celda(*c), AGG[c]["N"], AGG[c]["k"], f"{AGG[c]['chi2']:.1f}",
              fp(AGG[c]["p_f"]), f"{AGG[c]['F_id']:.1f}", fp(AGG[c]["p_id"]), f"{AGG[c]['cd']:.2f}"]
             for c in CELDAS]
    filas += [[f"{k[0]} {celda(*k[1:])}", r["N"], r["k"], f"{r['chi2']:.1f}", fp(r["p_f"]),
               f"{r['F_id']:.1f}", fp(r["p_id"]), f"{r['cd']:.2f}"] for k, r in SUB.items()]
    escribir("t4_friedman", cab, filas, "Pruebas de Friedman e Iman-Davenport, alpha = 0,05")

    # T5 · post-hoc: propuesta vs cada competidor (Shaffer y Nemenyi)
    cab = ["Competidor"] + [f"{celda(*c)} {t}" for c in CELDAS for t in ("Shaffer", "Nemenyi")]
    filas = []
    for a in comp:
        j = algs.index(a); fl = [et(a)]
        for c in CELDAS:
            r = AGG[c]; sgn = "+" if r["rm"][ip] < r["rm"][j] else "-"
            for cual in ("p_sh", "p_nem"):
                p = p_par(r, ip, j, cual)
                fl.append(fp(p) + (sgn if p < ALPHA else ""))
        filas.append(fl)
    escribir("t5_posthoc_propuesta", cab, filas,
             f"Post-hoc de Friedman: {et(P)} contra cada competidor (p ajustados)",
             nota="+ / -: diferencia significativa a favor / en contra de la propuesta (alpha = 0,05).")

    # T5b · todos los pares, los 24 rankings (solo CSV)
    filas = []
    for nom, r in [("CEC_71 " + celda(*c), AGG[c]) for c in CELDAS] + \
                  [(f"{k[0]} {celda(*k[1:])}", r) for k, r in SUB.items()]:
        for n, (i, j) in enumerate(r["pares"]):
            filas.append([nom, algs[i], algs[j], f"{r['rm'][i]:.3f}", f"{r['rm'][j]:.3f}",
                          f"{r['p_z'][n]:.3e}", f"{r['p_sh'][n]:.3e}", f"{r['p_nem'][n]:.3e}"])
    escribir("t5b_posthoc_todos_los_pares",
             ["ranking", "alg_i", "alg_j", "rango_i", "rango_j", "p_z", "p_shaffer", "p_nemenyi"],
             filas, "", tex=False)

    # T6 · % de competidores de los que difiere cada algoritmo
    cab = ["Algoritmo"] + [f"{celda(*c)} {t}" for c in CELDAS for t in ("Sh", "Nem")]
    filas = [[et(algs[i])] + [f"{pct_distintos(AGG[c], cual)[i]:.0f}" for c in CELDAS
                              for cual in ("p_sh", "p_nem")] for i in orden]
    escribir("t6_pct_distintos", cab, filas,
             "% de los demás algoritmos de los que difiere significativamente (Shaffer / Nemenyi)")

    # T7 · W/T/L de la propuesta contra cada competidor
    filas = []
    for a in comp:
        filas.append([et(a)] + ["/".join(map(str, wtl(load(P, *c), load(a, *c)))) for c in CELDAS])
    escribir("t7_wtl", ["Competidor"] + cab6[1:], filas,
             f"Victorias/empates/derrotas de {et(P)} (Wilcoxon rank-sum por función, alpha = 0,05)")

    # T8 · W/T/L por clase de función
    cab = ["Competidor", "Clase", "n"] + [celda(*c) for c in CELDAS]
    filas = []
    for a in comp:
        for cl in CLASES:
            ix = [p for p in range(71) if clase(p) == cl]
            filas.append([et(a), cl, len(ix)] +
                         ["/".join(map(str, wtl(load(P, *c), load(a, *c), ix))) for c in CELDAS])
    escribir("t8_wtl_por_clase", cab, filas, f"W/T/L de {et(P)} por clase de función",
             nota="Clases con la numeración oficial (opfunu omite la F2 de CEC2017).")

    # T9 · signed-rank pareado por problema (p global por par y celda), Holm sobre los competidores
    res = {c: [signed_rank(AGG[c]["X"][ip], AGG[c]["X"][algs.index(a)]) for a in comp] for c in CELDAS}
    adj = {c: holm([x[2] for x in res[c]]) for c in CELDAS}
    cab = ["Competidor"] + [f"{celda(*c)} {t}" for c in CELDAS for t in ("R+/R-", "p Holm")]
    filas = []
    for n, a in enumerate(comp):
        fl = [et(a)]
        for c in CELDAS:
            rp, rn, _ = res[c][n]
            fl += [f"{rp:.0f}/{rn:.0f}", fp(adj[c][n])]
        filas.append(fl)
    escribir("t9_signed_rank", cab, filas,
             f"Wilcoxon signed-rank sobre las 71 medias: {et(P)} contra cada competidor",
             nota="R+: rangos a favor de la propuesta. p ajustado por Holm dentro de cada celda.")

    # T10 · robustez: rangos con mediana
    filas = [[et(algs[i])] + [f"{MED[c]['rm'][i]:.2f} ({MED[c]['pos'][i]})" for c in CELDAS]
             for i in orden]
    escribir("t10_rangos_mediana", cab6, filas, "Robustez: rango medio de Friedman con la mediana")

    # Suplementario: estadísticos por función
    cab = ["algoritmo", "D", "maxfes", "benchmark", "funcion", "clase",
           "media", "mediana", "sd", "mejor", "peor"]
    filas = []
    for d, f in CELDAS:
        for a in algs:
            M = load(a, d, f)
            for pid in range(71):
                s, n = funcion(pid)
                filas.append([a, d, f, s, f"F{n}", clase(pid)] +
                             [f"{v:.4e}" for v in (M[pid].mean(), np.median(M[pid]), M[pid].std(ddof=1),
                                                   M[pid].min(), M[pid].max())])
    escribir("estadisticos_por_funcion", cab, filas, "", sub="suplementario", tex=False)
    for d, f in CELDAS:
        filas = []
        for pid in range(71):
            s, n = funcion(pid)
            mus = [load(a, d, f)[pid].mean() for a in algs]
            filas.append([f"{s} F{n}"] + [f"{m:.2e} ({load(a, d, f)[pid].std(ddof=1):.1e})"
                                          for a, m in zip(algs, mus)])
        escribir(f"medias_{DLAB[d]}_{BLAB[f]}", ["Función"] + [et(a) for a in algs], filas,
                 f"Media (SD) del error de 51 corridas, {celda(d, f)}", sub="suplementario")

    ablacion()
    parametros()
    figuras(algs, P, AGG)
    resumen(algs, P, AGG, comp)

# ── tabla de parámetros (valores leídos del código, no transcritos) ───────
def parametros():
    import inspect, competidores as C, marea2
    from experimento_config import POP_SIZE as N
    dfl = lambda fn: {k: v.default for k, v in inspect.signature(fn).parameters.items()
                      if v.default is not inspect.Parameter.empty}
    pre = lambda g: " / ".join(f"{g(f):,}".replace(",", " ") for f in CFG["maxfes"])
    clas = pre(lambda f: (f // N + 1) * N)                     # (T+1)·N
    m = marea2.config(CFG["propuesta"].replace("MAREA2_", ""))
    pso, de, cl, woa, fa = dfl(C.pso), dfl(C.de), dfl(C.clonalg), dfl(C.woa), dfl(C.fa)
    filas = [
        ["MAREA-2", f"{m.N0} -> {m.N_min} (logística de Cobetia)", pre(lambda f: f // N * N),
         f"logística xmid={m.cob_xmid}, s={m.cob_scal}, ventana {m.cob_tmax:g} h; lag={m.lag}; "
         f"alpha {m.alpha_min}-{m.alpha_max}, beta_max={m.beta_max}, gamma_max={m.gamma_max}; "
         f"dif. de vecinas + {m.diff_axis_mix} por eje Student-t nu={m.student_nu}; "
         f"HGT adaptativa (media inicial {m.hgt_cr}, sd {m.hgt_sd}, EMA {m.hgt_rate}); "
         f"quorum: objetivo {m.qs_target}, EMA {m.qs_rate}, ganancia {m.qs_gain}, "
         f"factor en [{m.qs_clip[0]}, {m.qs_clip[1]}]; inanición tras {m.limit_base}+2D iteraciones, "
         f"dispersa {m.p_disperse} con Lévy beta={m.levy_beta}, escala {m.disp_scale}(U-L); "
         f"memoria de {m.memory_size} cuencas, radio {m.depletion_radius}(U-L); greedy; límites por reflexión",
         "Stehlik et al. (2025), esta tesis", "Propuesta; sin colonias (ablación)"],
        ["L-SHADE", "18D -> 4 (LPSR)", clas,
         "H=6, M_F=M_CR=0.5, archivo 2.6N, p=0.11, current-to-pbest/1",
         "Tanabe y Fukunaga (2014)", "Población propia"],
        ["jSO", "round(25 ln(D) sqrt(D)) -> 4 (LPSR)", clas,
         "H=5, M_F=0.3, M_CR=0.8, archivo 1.0N, p 0.125->0.25, current-to-pbest-w/1",
         "Brest et al. (2017)", "Población propia"],
        ["DE", N, clas, f"DE/rand/1/bin, F={de['F_scale']}, CR={de['CR']}", "Storn y Price (1997)", "N=80"],
        ["PSO", N, clas, f"constricción, w={pso['w']}, c1={pso['c1']}, c2={pso['c2']}",
         "Clerc y Kennedy (2002)", "N=80"],
        ["ABC", f"{N} ({N // 2} fuentes)", clas, "limit = fuentes·D/2", "Karaboga y Basturk (2007)", "N=80"],
        ["GWO", N, clas, "a: 2 -> 0 lineal", "Mirjalili et al. (2014)", "N=80"],
        ["WOA", N, clas, f"b={woa['b']}, a: 2 -> 0", "Mirjalili y Lewis (2016)", "N=80"],
        ["FA", N, clas, f"alpha={fa['alpha']}, beta0={fa['beta0']}, gamma={fa['gamma']}/(D(U-L)^2)",
         "Yang (2010)", "N=80; gamma reescalado (underflow)"],
        ["CLONALG", N, clas, f"beta={cl['beta']}, máx. {cl['n_clones']} clones, reemplazo d={cl['d']}",
         "de Castro y Von Zuben (2002)", "N=80"],
    ]
    escribir("t12_parametros", ["Algoritmo", "Población", "Evaluaciones 5k / 50k / 500k",
                                "Parámetros", "Fuente", "Desviación del setup §4"], filas,
             "Parámetros de los algoritmos",
             nota="MAREA-2 recibe N·T evaluaciones y los demás (T+1)·N, con T = MaxFES // N: 80 menos para "
                  "la propuesta (-1,6 % a 5k; sesgo conservador, declarado). 'N=80': tamaño uniforme en los 7 clásicos, decisión de diseño.")

# ── tabla de ablación (15 corridas, semillas 0-14 en todos los brazos) ────
def ablacion():
    A = CFG["ablacion"]; filas = []
    for comp, fase, con, sin, decision in A["filas"]:
        celdas = []
        for d, f in CELDAS:
            try: celdas.append("/".join(map(str, wtl(load(con, d, f, abl=True), load(sin, d, f, abl=True)))))
            except FileNotFoundError: celdas.append("—")
        filas.append([comp, fase, con, sin] + celdas + [decision])
    escribir("t11_ablacion", ["Componente", "Fase", "Con", "Sin"] + [celda(d, f) for d, f in CELDAS]
             + ["Decisión"], filas, "Ablación de MAREA-2 (W/T/L, rank-sum, 15 corridas por brazo)",
             nota=A["nota"] + " 'ruido_normal' y 'eje_student' son perturbación por eje Student-t, no gaussiana.")

# ── figuras ────────────────────────────────────────────────────────────────
def figuras(algs, P, AGG):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
    INK, INK2, INK3, SURF, GRIS = "#0b0b0b", "#52514e", "#8a8985", "#fcfcfb", "#c3c2bd"
    plt.rcParams.update({
        "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
        "font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": INK3,
        "axes.linewidth": 0.6, "axes.labelcolor": INK2, "text.color": INK,
        "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
        "axes.spines.right": False})
    COL = {P: BLUE, "LSHADE": ORANGE, "JSO": AQUA, "DE": VIOLET}
    os.makedirs(f"{OUT}/figuras", exist_ok=True)
    def guardar(fig, n):
        for ext in ("png", "pdf"): fig.savefig(f"{OUT}/figuras/{n}.{ext}", dpi=220, bbox_inches="tight",
                                         metadata={"CreationDate": None} if ext == "pdf" else None)
        plt.close(fig)

    # Fig. A · rango medio por MaxFES, D10 y D50; etiquetas directas
    fig, axs = plt.subplots(1, 2, figsize=(9, 4.2), sharey=True)
    xs = np.arange(len(CFG["maxfes"]))
    for ax, d in zip(axs, CFG["dims"]):
        Y = {a: [AGG[(d, f)]["rm"][i] for f in CFG["maxfes"]] for i, a in enumerate(algs)}
        for a in sorted(algs, key=lambda a: a in COL):
            c = COL.get(a, GRIS); lw = 2.4 if a == P else 2 if a in COL else 1.2
            ax.plot(xs, Y[a], color=c, lw=lw, marker="o", ms=4.5 if a in COL else 3,
                    zorder=3 if a in COL else 2)
        ult = sorted(algs, key=lambda a: Y[a][-1]); yl = []
        for a in ult:                                   # separar etiquetas del borde derecho
            y = Y[a][-1]
            if yl and y - yl[-1] < 0.42: y = yl[-1] + 0.42
            yl.append(y)
            ax.text(xs[-1] + 0.1, y, et(a), va="center", fontsize=8,
                    color=INK if a in COL else INK2, fontweight="bold" if a == P else "normal")
        ax.set_xticks(xs, [BLAB[f] for f in CFG["maxfes"]]); ax.set_xlim(-0.15, xs[-1] + 0.75)
        ax.set_title("D10" if d == 10 else "D50 (CEC2022 en D20)", fontsize=10, loc="left")
        ax.set_xlabel("MaxFES"); ax.grid(axis="y", color="#e6e5e1", lw=0.6)
    axs[0].set_ylabel("Rango medio de Friedman (menor = mejor)"); axs[0].invert_yaxis()
    axs[0].set_ylim(len(algs) + 0.4, 0.6)
    guardar(fig, "fig_rangos")

    # Fig. B · diagramas de diferencia crítica (Nemenyi), 6 celdas
    fig, axs = plt.subplots(3, 2, figsize=(10, 8.4))
    for ax, (d, f) in zip(axs.T.flat, CELDAS):
        r = AGG[(d, f)]; k = r["k"]; rm = r["rm"]; o = np.argsort(rm); cd = r["cd"]
        ax.set_xlim(0.5, k + 0.5); ax.set_ylim(-1.2 - (k // 2) * 0.5, 1.3); ax.axis("off")
        ax.plot([1, k], [0, 0], color=INK3, lw=0.8)
        for t in range(1, k + 1):
            ax.plot([t, t], [0, 0.12], color=INK3, lw=0.8)
            ax.text(t, 0.25, str(t), ha="center", fontsize=7.5, color=INK2)
        ax.plot([1, 1 + cd], [0.95, 0.95], color=INK, lw=1.2)
        for x in (1, 1 + cd): ax.plot([x, x], [0.88, 1.02], color=INK, lw=1.2)
        ax.text(1 + cd / 2, 1.08, f"CD = {cd:.2f}", ha="center", fontsize=7.5, color=INK2)
        mitad = (k + 1) // 2
        for n, i in enumerate(o):
            izq = n < mitad; fila = n if izq else k - 1 - n
            y = -0.45 - fila * 0.5; xt = 0.55 if izq else k + 0.45
            a = algs[i]; c = COL.get(a, INK2)
            ax.plot([rm[i], rm[i], xt], [0, y, y], color=c, lw=1.6 if a in COL else 0.8)
            ax.text(xt + (-0.05 if izq else 0.05), y, f"{et(a)} {rm[i]:.2f}",
                    ha="right" if izq else "left", va="center", fontsize=7.5,
                    color=INK if a in COL else INK2, fontweight="bold" if a == P else "normal")
        # grupos sin diferencia significativa: tramos maximales con rango_max - rango_min < CD
        grupos = []
        for s in range(k):
            e = s
            while e + 1 < k and rm[o[e + 1]] - rm[o[s]] < cd: e += 1
            if e > s and not any(gs <= s and e <= ge for gs, ge in grupos): grupos.append((s, e))
        for g, (s, e) in enumerate(grupos):
            yg = -0.18 - g * 0.11
            ax.plot([rm[o[s]] - 0.05, rm[o[e]] + 0.05], [yg, yg], color=INK, lw=2.2,
                    solid_capstyle="round")
        ax.set_title(f"{celda(d, f)}  (Friedman p = {fp(r['p_f'])})", fontsize=9.5, loc="left")
    fig.text(0.01, 0.005, "Barras negras: algoritmos sin diferencia significativa (Nemenyi, alpha = 0,05). "
             "D50: CEC2022 en D20.", fontsize=7.5, color=INK2)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    guardar(fig, "fig_cd_nemenyi")

    # Fig. C · W/T/L de la propuesta contra cada competidor (barras apiladas por celda)
    comp = [a for a in algs if a != P]
    fig, axs = plt.subplots(2, 3, figsize=(10.5, 5.6), sharex=True, sharey=True)
    for ax, (d, f) in zip(axs.flat, CELDAS):
        W = np.array([wtl(load(P, d, f), load(a, d, f)) for a in comp])
        y = np.arange(len(comp))[::-1]
        izq = np.zeros(len(comp))
        for v, c, lab in zip(W.T, (BLUE, "#d8d7d3", "#e34948"), ("Gana", "Empata", "Pierde")):
            ax.barh(y, v, left=izq, color=c, height=0.72, edgecolor=SURF, linewidth=1, label=lab)
            izq += v
        for yy, (w, t, l) in zip(y, W):
            ax.text(72, yy, f"{w}/{t}/{l}", va="center", fontsize=6.8, color=INK2)
        ax.set_yticks(y, [et(a) for a in comp], fontsize=8); ax.set_xlim(0, 85)
        ax.set_title(celda(d, f), fontsize=9.5, loc="left"); ax.tick_params(axis="y", length=0)
    axs[0, 0].legend(ncol=3, fontsize=7.5, frameon=False, loc="lower left", bbox_to_anchor=(0, 1.12))
    for ax in axs[1]: ax.set_xlabel("Funciones (de 71)")
    fig.tight_layout()
    guardar(fig, "fig_wtl")

# ── resumen en consola ─────────────────────────────────────────────────────
def resumen(algs, P, AGG, comp):
    ip = algs.index(P)
    print(f"Propuesta: {et(P)}  ·  {len(algs)} algoritmos  ·  media de 51 corridas")
    for c in CELDAS:
        r = AGG[c]; o = np.argsort(r["rm"])
        peores = [et(algs[j]) for j in range(r["k"])
                  if j != ip and p_par(r, ip, j, "p_nem") >= ALPHA]
        print(f"  {celda(*c):9s} {r['pos'][ip]}º ({r['rm'][ip]:.2f})  1º {et(algs[o[0]])} "
              f"{r['rm'][o[0]]:.2f}  CD={r['cd']:.2f}  sin dif. Nemenyi con: {', '.join(peores) or '-'}")
    print(f"Salida en {OUT}/tablas, {OUT}/suplementario, {OUT}/figuras")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--regresion", action="store_true", help="solo el control de la tabla 2.1")
    a = ap.parse_args()
    ok = regresion()
    if a.regresion: sys.exit(0 if ok else 1)
    if not ok: sys.exit("El control de regresión falla: no se generan tablas")
    main()
