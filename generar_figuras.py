import json, os, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import ranksums

D = "resultados_experimento"
OUT = "figuras"

# ── Paleta validada (dataviz, light mode) ────────────────────────────────────
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
RED, GRAY_MID = "#e34948", "#d8d7d3"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8985"
SURF = "#fcfcfb"
GRAY_LINE = "#c3c2bd"

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.edgecolor": INK3, "axes.linewidth": 0.6,
    "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
    "axes.spines.top": False, "axes.spines.right": False,
})

def load(a, d, f):
    p = f"{D}/{a}_D{d}_FES{f}.json"
    if not os.path.exists(p): return None
    raw = json.load(open(p)); M = np.full((71, 51), np.nan)
    for k, v in raw.items():
        M[int(k.split("_")[0][3:]), int(k.split("run")[1])] = v
    return M

def wtl(A, B, alpha=0.05):
    w = t = l = 0
    for i in range(71):
        p = ranksums(A[i], B[i]).pvalue
        if p >= alpha: t += 1
        elif np.median(A[i]) < np.median(B[i]): w += 1
        else: l += 1
    return w, t, l

BUDGETS = [5000, 50000, 500000]
BLAB = {5000: "5k", 50000: "50k", 500000: "500k"}

# ══════════════════════════════════════════════════════════════════════════
# FIG 1 — la inversión bio → base
# ══════════════════════════════════════════════════════════════════════════
def fig_inversion(pares, titulo, subt, fname, etq_izq, etq_der):
    filas = []
    for (nomA, nomB, lab) in pares:
        for d in (10, 50):
            for fes in BUDGETS:
                A, B = load(nomA, d, fes), load(nomB, d, fes)
                if A is None or B is None: continue
                w, t, l = wtl(A, B)
                filas.append((f"{lab}  D{d} · {BLAB[fes]}", w, t, l))
    n = len(filas)
    fig, ax = plt.subplots(figsize=(7.4, 0.42 * n + 1.5))
    y = np.arange(n)[::-1]
    G = 0.55  # hueco de 2px aprox entre segmentos
    for yi, (_, w, t, l) in zip(y, filas):
        ax.barh(yi, w, color=BLUE, height=0.62, zorder=3)
        ax.barh(yi, t, left=w + G, color=GRAY_MID, height=0.62, zorder=3)
        ax.barh(yi, l, left=w + t + 2 * G, color=RED, height=0.62, zorder=3)
        if w >= 5: ax.text(w / 2, yi, str(w), ha="center", va="center",
                           color="white", fontsize=8, fontweight="bold", zorder=4)
        if l >= 5: ax.text(w + t + 2 * G + l / 2, yi, str(l), ha="center", va="center",
                           color="white", fontsize=8, fontweight="bold", zorder=4)
    ax.axvline(35.5, color=INK3, lw=0.7, ls=(0, (3, 3)), zorder=2)
    ax.set_yticks(y); ax.set_yticklabels([f[0] for f in filas], fontsize=8.5)
    ax.set_xlim(0, 73); ax.set_xlabel("funciones (de 71)", fontsize=8.5)
    ax.grid(axis="x", color=GRAY_LINE, lw=0.5, alpha=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.legend(handles=[Patch(color=BLUE, label=etq_izq),
                       Patch(color=GRAY_MID, label="sin diferencia (p≥0,05)"),
                       Patch(color=RED, label=etq_der)],
              loc="lower center", bbox_to_anchor=(0.5, -0.055 - 0.30/ (0.42*n+1.5) * 3),
              ncol=3, frameon=False, fontsize=8.5)
    ax.set_title(titulo, fontsize=11.5, fontweight="bold", color=INK, loc="left", pad=24)
    ax.text(0, 1.028, subt, transform=ax.transAxes, fontsize=8.5, color=INK2, va="bottom")
    fig.tight_layout()
    fig.savefig(f"{OUT}/{fname}", dpi=220, bbox_inches="tight")
    plt.close(fig)
    return filas

f1 = fig_inversion(
    [("MAREABIO_REL", "MAREA_REL", "")],
    "El ciclo biológico ayuda sólo con presupuesto corto",
    "MAREABIO_REL vs MAREA_REL (mismo escalado) · Wilcoxon rank-sum por función, α=0,05",
    "fig1_inversion.png", "gana MAREA-bio", "gana MAREA-base")

f4 = fig_inversion(
    [("MAREA_REL", "MAREA", "base"), ("MAREABIO_REL", "MAREABIO", "bio ")],
    "El reescalado del ruido: decisivo en la base, se agota en bio",
    "variante reescalada vs unidades absolutas · Wilcoxon rank-sum por función, α=0,05",
    "fig4_escalado.png", "gana la versión reescalada", "gana la versión original")

# ══════════════════════════════════════════════════════════════════════════
# FIG 2 — rango medio por presupuesto
# ══════════════════════════════════════════════════════════════════════════
ALGS = ["MAREA", "MAREA_REL", "MAREABIO", "MAREABIO_REL",
        "PSO", "DE", "ABC", "GWO", "WOA", "FA", "CLONALG"]
COL = {"MAREA_REL": BLUE, "MAREA": ORANGE, "MAREABIO_REL": AQUA, "MAREABIO": VIOLET}

def rangos(d, fes):
    med, names = [], []
    for a in ALGS:
        M = load(a, d, fes)
        if M is None: continue
        med.append(np.median(M, axis=1)); names.append(a)
    med = np.array(med)
    R = np.empty_like(med)
    for j in range(med.shape[1]):
        col = med[:, j]
        srt = np.argsort(col, kind="stable")
        r = np.empty(len(col)); r[srt] = np.arange(1, len(col) + 1)
        u, inv, cnt = np.unique(col, return_inverse=True, return_counts=True)
        for ui in np.where(cnt > 1)[0]:
            m = inv == ui; r[m] = r[m].mean()
        R[:, j] = r
    return dict(zip(names, R.mean(axis=1)))

def separar(vals, minsep):
    """Empuja etiquetas para que no se solapen, preservando el orden."""
    orden = sorted(range(len(vals)), key=lambda i: vals[i])
    y = [vals[i] for i in orden]
    for k in range(1, len(y)):
        if y[k] - y[k-1] < minsep: y[k] = y[k-1] + minsep
    exceso = y[-1] - max(vals)
    if exceso > 0:
        for k in range(len(y)): y[k] -= exceso * 0.5
    out = [0]*len(vals)
    for pos, i in enumerate(orden): out[i] = y[pos]
    return out

fig, axes = plt.subplots(1, 2, figsize=(10.6, 5.2), sharey=True)
xs = np.arange(3)
for ax, d in zip(axes, (10, 50)):
    serie = {a: [] for a in ALGS}
    for fes in BUDGETS:
        r = rangos(d, fes)
        for a in ALGS: serie[a].append(r.get(a, np.nan))
    finales = [serie[a][2] for a in ALGS]
    yl = separar(finales, 0.42)
    for a in ALGS:
        if a in COL: continue
        ax.plot(xs, serie[a], color=GRAY_LINE, lw=1.4, marker="o", ms=4,
                zorder=2, solid_capstyle="round")
    for a in ALGS:
        if a not in COL: continue
        ax.plot(xs, serie[a], color=COL[a], lw=2.2, marker="o", ms=7,
                zorder=4, solid_capstyle="round",
                markeredgecolor=SURF, markeredgewidth=1.6)
    for a, yy in zip(ALGS, yl):
        es = a in COL
        ax.annotate(a, xy=(2, serie[a][2]), xytext=(2.14, yy),
                    fontsize=8.5 if es else 7.5,
                    color=COL[a] if es else INK3,
                    fontweight="bold" if es else "normal",
                    va="center", ha="left", zorder=5,
                    arrowprops=dict(arrowstyle="-", lw=0.6,
                                    color=COL[a] if es else GRAY_LINE,
                                    shrinkA=0, shrinkB=2))
    ax.set_xticks(xs); ax.set_xticklabels([BLAB[b] for b in BUDGETS])
    ax.set_xlim(-0.15, 3.35); ax.set_xlabel("presupuesto (MaxFES)", fontsize=8.5)
    ax.grid(axis="y", color=GRAY_LINE, lw=0.5, alpha=0.55, zorder=0)
    ax.set_axisbelow(True)
    ax.set_title(f"D = {d}", fontsize=10, color=INK2, loc="left")
axes[0].invert_yaxis()
axes[0].set_ylabel("rango medio sobre 71 funciones  (1 = mejor)", fontsize=8.5)
fig.suptitle("MAREA_REL es 2\u00ba de 11 en alta dimensi\u00f3n; en D10 no destaca",
             fontsize=11.5, fontweight="bold", color=INK, x=0.008, ha="left", y=1.005)
fig.text(0.008, 0.955, "rango medio de Friedman \u00b7 71 funciones CEC 2014/2017/2022 \u00b7 51 corridas",
         fontsize=8.5, color=INK2, ha="left")
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(f"{OUT}/fig2_rangos.png", dpi=220, bbox_inches="tight")
plt.close(fig)

# ══════════════════════════════════════════════════════════════════════════
# FIG 3 — colapso de B~
# ══════════════════════════════════════════════════════════════════════════
import sys; sys.path.insert(0, os.path.abspath("."))
from marea_bio_benchmark import MAREABioParams, indicadores
import marea_bio_benchmark as MBB

def traza_B(n, seed=1, T=40):
    p = MAREABioParams(); p.N = 80; p.T = T
    rng = np.random.default_rng(seed)
    L, U = np.full(n, -100.0), np.full(n, 100.0)
    f = lambda x: float(np.sum(x**2))
    X = rng.uniform(L, U, size=(p.N, n)); F = np.array([f(x) for x in X])
    isp = np.mean(np.std(X, axis=0)) + 1e-12; isf = np.std(F) + 1e-12
    fg = F.min(); Xg = X[F.argmin()].copy(); no_imp = 0; Bs = []
    for t in range(1, T + 1):
        mu_t, B_t, V_t = indicadores(X, F, no_imp, isp, isf, p)
        Bs.append(B_t)
        mu_pop = X.mean(axis=0); Xn = np.empty_like(X)
        for i in range(p.N):
            xi = rng.uniform(-1, 1, n)
            if B_t > p.B_threshold:
                Xn[i] = np.clip(X[i] + p.a2*B_t*(Xg-X[i]) + p.c2*V_t*(mu_pop-X[i]) + p.sigma2*xi, L, U)
            else:
                Xn[i] = np.clip(X[i] + p.a1*mu_t*(Xg-X[i]) + p.sigma1*(1-B_t)*xi, L, U)
        X = Xn; F = np.array([f(x) for x in X])
        if F.min() < fg: fg = F.min(); Xg = X[F.argmin()].copy(); no_imp = 0
        else: no_imp += 1
    return np.array(Bs)

fig, ax = plt.subplots(figsize=(7.4, 4.3))
NT = 15
trazas = [(10, BLUE, "D = 10", "-"), (50, ORANGE, "D = 50", (0, (5, 2)))]
for n, col, lab, ls in trazas:
    B = traza_B(n, T=NT)
    ax.plot(np.arange(1, len(B)+1), B, color=col, lw=2.2, ls=ls, marker="o", ms=6,
            markeredgecolor=SURF, markeredgewidth=1.4, zorder=4, label=lab,
            solid_capstyle="round")
ax.axhline(0.45, color=RED, lw=1.4, ls=(0, (4, 3)), zorder=3)
ax.text(NT + 0.3, 0.45, "umbral 0,45", fontsize=8, color=RED, ha="right", va="bottom")
ax.axvspan(0.55, 1.45, color=GRAY_MID, alpha=0.7, zorder=1)
ax.annotate("Fase I dura\n1 iteraci\u00f3n", xy=(1.45, 0.30), xytext=(2.6, 0.22),
            fontsize=8.5, color=INK2, va="center", ha="left",
            arrowprops=dict(arrowstyle="->", color=INK3, lw=0.9))
ax.annotate("B~ = 0,70 en la iteraci\u00f3n 2", xy=(2, 0.70), xytext=(4.2, 0.60),
            fontsize=8.5, color=INK2, va="center", ha="left",
            arrowprops=dict(arrowstyle="->", color=INK3, lw=0.9))
ax.set_xlim(0.55, NT + 0.45); ax.set_ylim(-0.04, 1.09)
ax.set_xticks(range(1, NT + 1, 2))
ax.set_xlabel("iteraci\u00f3n", fontsize=8.5)
ax.set_ylabel("B~  (indicador de biofilm)", fontsize=8.5)
ax.grid(color=GRAY_LINE, lw=0.5, alpha=0.55, zorder=0); ax.set_axisbelow(True)
ax.legend(loc="lower right", frameon=False, fontsize=8.5, ncol=2,
          bbox_to_anchor=(1.0, 0.03))
ax.set_title("B~ cruza el umbral en la iteraci\u00f3n 2 y satura: la Fase I no existe",
             fontsize=11.5, fontweight="bold", color=INK, loc="left", pad=22)
ax.text(0, 1.035,
        "a1 = 0,7 contrae la poblaci\u00f3n al 30 % del spread inicial en un paso  \u2192  B~ = 1 \u2212 0,30 = 0,70",
        transform=ax.transAxes, fontsize=8.5, color=INK2, va="bottom")
fig.tight_layout()
fig.savefig(f"{OUT}/fig3_colapso_fases.png", dpi=220, bbox_inches="tight")
plt.close(fig)
print("listo")
