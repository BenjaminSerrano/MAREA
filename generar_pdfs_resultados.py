#!/usr/bin/env python3
"""
Consolida resultados_experimento/*.json en un PDF por algoritmo.

Un PDF por algoritmo (11 en total) en pdfs_resultados/. Cada uno lleva:
  - portada con el protocolo y el resumen de rangos de Friedman por celda
  - una tabla por dimension: 71 funciones x (3 presupuestos x best/worst/median/mean/std)

Los valores son ERROR al optimo, f(x) - f*, tal como los guarda experimento_config.
Idempotente: se vuelve a correr entero cuando cambian los datos.
"""
import json, os, re, sys
from collections import defaultdict
import numpy as np
from scipy.stats import rankdata
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.units import mm
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import (SimpleDocTemplate, Table, TableStyle, Paragraph,
                                Spacer, PageBreak)

R        = "resultados_experimento"
OUT      = "pdfs_resultados"
ALGOS    = ["MAREA", "MAREA_REL", "MAREABIO", "MAREABIO_REL",
            "PSO", "DE", "ABC", "GWO", "WOA", "FA", "CLONALG"]
DIMS     = [10, 50]
BUDGETS  = [5000, 50000, 500000]
N_RUNS   = 51

INK, INK2 = colors.HexColor("#1a1a1a"), colors.HexColor("#5a5a5a")
AZUL      = colors.HexColor("#2a78d6")
FONDO     = colors.HexColor("#f0f4f9")

NOTAS = {
    "MAREA":        "Base, tres fuerzas con horario fijo en t/T. gamma_max=0.2 en unidades "
                    "ABSOLUTAS: el ruido queda inerte en el dominio CEC [-100,100].",
    "MAREA_REL":    "Base con el ruido reescalado por (U-L)/10.24. El reescalado es el efecto "
                    "dominante del estudio: gana 68-70 de 71 funciones contra MAREA.",
    "MAREABIO":     "Tres fases biologicas (planctonica / biofilm / desprendimiento), sigma sin "
                    "reescalar. La Fase I dura una sola iteracion (a1=0.7).",
    "MAREABIO_REL": "Tres fases biologicas con sigma_rel=True. Gana a MAREA_REL a 5k FES y "
                    "pierde a 50k: el ciclo biologico es arranque rapido, no convergencia.",
    "FA":           "Firefly Algorithm. CORREGIDO 2026-09-04: gastaba 39.6x el MaxFES y su "
                    "atraccion hacia underflow a 0.0 (gamma sin reescalar). Estos resultados "
                    "son los posteriores al arreglo.",
    "CLONALG":      "Algoritmo clonal. CORREGIDO 2026-09-04: gastaba 3.5x el MaxFES. Estos "
                    "resultados usan el presupuesto honesto, por eso son PEORES que los "
                    "previos: el arreglo le quito presupuesto, no le arreglo un mecanismo.",
}

def cargar(alg, dim, bud):
    with open(f"{R}/{alg}_D{dim}_FES{bud}.json") as fh:
        d = json.load(fh)
    g = defaultdict(list)
    for k, v in d.items():
        g[int(re.match(r"pid(\d+)_run", k).group(1))].append(float(v))
    return {p: np.array(v) for p, v in g.items()}

def registro_funciones():
    """pid -> (nombre, suite, dim efectiva por dimension nominal). Sin instanciar."""
    import inspect
    import opfunu.cec_based.cec2014 as c14
    import opfunu.cec_based.cec2017 as c17
    import opfunu.cec_based.cec2022 as c22
    out = []
    for mod, pref, dims in ((c14, "CEC2014", [10, 50]),
                            (c17, "CEC2017", [10, 50]),
                            (c22, "CEC2022", [10, 20])):
        cls = sorted([(n, o) for n, o in inspect.getmembers(mod, inspect.isclass)
                      if n.startswith("F") and n[1].isdigit()],
                     key=lambda x: int("".join(filter(str.isdigit,
                                                      x[0].replace(pref, "")))))
        out += [{"name": n, "suite": pref, "dims": dims} for n, _ in cls]
    return out

def num(x):
    """Notacion cientifica de ancho fijo — el formato estandar de las tablas CEC."""
    if not np.isfinite(x): return "n/f"
    if x == 0:             return "0.00e+00"
    return f"{x:.2e}"

def rangos_friedman():
    """(dim,bud) -> {alg: (rango_medio, posicion)}"""
    out = {}
    for dim in DIMS:
        for bud in BUDGETS:
            M = {a: cargar(a, dim, bud) for a in ALGOS}
            pids = sorted(set.intersection(*[set(M[a]) for a in ALGOS]))
            mat = np.array([[np.median(M[a][p]) for a in ALGOS] for p in pids])
            mr = np.apply_along_axis(rankdata, 1, mat).mean(0)
            orden = np.argsort(mr)
            pos = {ALGOS[j]: i + 1 for i, j in enumerate(orden)}
            out[(dim, bud)] = {a: (mr[i], pos[a]) for i, a in enumerate(ALGOS)}
    return out

def gana_de_11(alg, dim, bud):
    """En cuantas funciones este algoritmo tiene la mejor mediana de los 11."""
    M = {a: cargar(a, dim, bud) for a in ALGOS}
    pids = sorted(set.intersection(*[set(M[a]) for a in ALGOS]))
    n = 0
    for p in pids:
        meds = {a: np.median(M[a][p]) for a in ALGOS}
        if meds[alg] <= min(meds.values()):
            n += 1
    return n, len(pids)

def construir(alg, funcs, rk, ganas):
    ss = getSampleStyleSheet()
    h1 = ParagraphStyle("h1", parent=ss["Title"], fontSize=17, textColor=INK,
                        alignment=0, spaceAfter=2)
    sub = ParagraphStyle("sub", parent=ss["Normal"], fontSize=9.5, textColor=INK2,
                         spaceAfter=10)
    cuerpo = ParagraphStyle("b", parent=ss["Normal"], fontSize=8.5, textColor=INK,
                            leading=12, spaceAfter=6)
    h2 = ParagraphStyle("h2", parent=ss["Normal"], fontSize=11, textColor=INK,
                        fontName="Helvetica-Bold", spaceBefore=8, spaceAfter=5)

    doc = SimpleDocTemplate(f"{OUT}/{alg}.pdf", pagesize=landscape(A4),
                            leftMargin=12*mm, rightMargin=12*mm,
                            topMargin=12*mm, bottomMargin=12*mm,
                            title=f"Resultados CEC — {alg}", author="Benjamin Serrano")
    el = []
    el.append(Paragraph(f"{alg} — resultados CEC", h1))
    el.append(Paragraph("Benchmark MAREA · paper OE1 · 71 funciones CEC2014/2017/2022 · "
                        f"{N_RUNS} corridas por funcion", sub))

    el.append(Paragraph("Protocolo", h2))
    el.append(Paragraph(
        "Valores = <b>error al optimo</b>, f(x) - f*. Menor es mejor; 0 es el optimo global. "
        "Dominio [-100,100] en todas las funciones. Poblacion N=80, iteraciones T=MaxFES/N. "
        "Presupuestos 5.000 / 50.000 / 500.000 evaluaciones. Las 30 funciones de CEC2014 y las "
        "29 de CEC2017 corren en D=10 y D=50; las 12 de CEC2022 solo llegan a D=20, asi que en "
        "la columna D50 corren en <b>D=20</b> (marcado en la tabla). Todos los resultados vienen "
        "del mismo entorno del Mac mini (python 3.11.9, numpy 2.3.3, opfunu 1.0.4).", cuerpo))

    if alg in NOTAS:
        el.append(Paragraph("Nota sobre este algoritmo", h2))
        el.append(Paragraph(NOTAS[alg], cuerpo))

    el.append(Paragraph("Resumen — posicion entre los 11 algoritmos", h2))
    fil = [["dimension", "presupuesto", "rango medio de Friedman", "posicion de 11",
            "funciones donde es el mejor"]]
    for dim in DIMS:
        for bud in BUDGETS:
            mr, pos = rk[(dim, bud)][alg]
            g, tot = ganas[(dim, bud)]
            fil.append([f"D{dim}", f"{bud//1000}k FES", f"{mr:.2f}",
                        f"{pos}º", f"{g} de {tot}"])
    t = Table(fil, hAlign="LEFT", colWidths=[24*mm, 28*mm, 46*mm, 26*mm, 42*mm])
    t.setStyle(TableStyle([
        ("FONT", (0,0), (-1,0), "Helvetica-Bold", 8),
        ("FONT", (0,1), (-1,-1), "Helvetica", 8),
        ("TEXTCOLOR", (0,0), (-1,0), colors.white),
        ("BACKGROUND", (0,0), (-1,0), AZUL),
        ("ROWBACKGROUNDS", (0,1), (-1,-1), [colors.white, FONDO]),
        ("ALIGN", (2,1), (-1,-1), "RIGHT"),
        ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
        ("GRID", (0,0), (-1,-1), 0.25, colors.HexColor("#c8d0d8")),
        ("TOPPADDING", (0,0), (-1,-1), 3), ("BOTTOMPADDING", (0,0), (-1,-1), 3),
    ]))
    el.append(t)

    # ---- tablas de estadisticas, una por dimension ----
    for dim in DIMS:
        el.append(PageBreak())
        el.append(Paragraph(f"{alg} — estadisticas por funcion, D = {dim}", h1))
        el.append(Paragraph("error al optimo sobre 51 corridas · b=best  w=worst  "
                            "med=median  avg=mean  sd=std", sub))
        datos = {b: cargar(alg, dim, b) for b in BUDGETS}
        enc1 = ["func", "suite", "dim"]
        for b in BUDGETS:
            enc1 += [f"{b//1000}k", "", "", "", ""]
        enc2 = ["", "", ""] + ["b", "w", "med", "avg", "sd"] * 3
        fil = [enc1, enc2]
        for pid, fmeta in enumerate(funcs):
            dim_ef = 20 if (dim == 50 and 50 not in fmeta["dims"]) else dim
            fila = [fmeta["name"], fmeta["suite"].replace("CEC", ""), str(dim_ef)]
            for b in BUDGETS:
                v = datos[b].get(pid)
                if v is None:
                    fila += ["-"] * 5
                else:
                    fila += [num(v.min()), num(v.max()), num(np.median(v)),
                             num(v.mean()), num(v.std(ddof=1))]
            fil.append(fila)
        anchos = [17*mm, 13*mm, 9*mm] + [15.5*mm] * 15
        t = Table(fil, repeatRows=2, hAlign="LEFT", colWidths=anchos)
        est = [
            ("FONT", (0,0), (-1,1), "Helvetica-Bold", 6.2),
            ("FONT", (0,2), (-1,-1), "Helvetica", 5.8),
            ("TEXTCOLOR", (0,0), (-1,1), colors.white),
            ("BACKGROUND", (0,0), (-1,1), AZUL),
            ("ROWBACKGROUNDS", (0,2), (-1,-1), [colors.white, FONDO]),
            ("ALIGN", (2,0), (-1,-1), "RIGHT"),
            ("ALIGN", (0,0), (1,-1), "LEFT"),
            ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
            ("GRID", (0,0), (-1,-1), 0.2, colors.HexColor("#c8d0d8")),
            ("TOPPADDING", (0,0), (-1,-1), 1.1),
            ("BOTTOMPADDING", (0,0), (-1,-1), 1.1),
            ("LEFTPADDING", (0,0), (-1,-1), 2), ("RIGHTPADDING", (0,0), (-1,-1), 2),
        ]
        for i, b in enumerate(BUDGETS):          # une la cabecera de cada presupuesto
            c0 = 3 + i*5
            est += [("SPAN", (c0,0), (c0+4,0)), ("ALIGN", (c0,0), (c0+4,0), "CENTER")]
            if i % 2 == 1:
                est.append(("BACKGROUND", (c0,2), (c0+4,-1), colors.HexColor("#eef1f4")))
        t.setStyle(TableStyle(est))
        el.append(t)

    doc.build(el)

def main():
    os.makedirs(OUT, exist_ok=True)
    funcs = registro_funciones()
    print(f"{len(funcs)} funciones en el registro")
    print("calculando rangos de Friedman...")
    rk = rangos_friedman()
    for alg in ALGOS:
        ganas = {(d, b): gana_de_11(alg, d, b) for d in DIMS for b in BUDGETS}
        construir(alg, funcs, rk, ganas)
        print(f"  {OUT}/{alg}.pdf")
    print("listo")

if __name__ == "__main__":
    main()
