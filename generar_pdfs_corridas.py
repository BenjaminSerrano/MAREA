#!/usr/bin/env python3
"""
Un PDF por algoritmo con las 51 CORRIDAS INDIVIDUALES de cada funcion.

Complementa generar_pdfs_resultados.py (que solo lleva estadisticos resumidos).
Salida en pdfs_corridas/. Por algoritmo: 6 celdas (2 dims x 3 presupuestos),
cada una con las 71 funciones x 51 corridas, partidas en dos bloques de
columnas para que quepan a lo ancho.

Los valores son ERROR al optimo, f(x) - f*, redondeados a %.3e (4 cifras
significativas). Los JSON de resultados_experimento/ siguen siendo la fuente
de verdad en precision completa.

Idempotente: se vuelve a correr entero cuando cambian los datos.
"""
import json, os, re
from collections import defaultdict
import numpy as np
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.units import mm
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import (SimpleDocTemplate, Table, TableStyle, Paragraph,
                                PageBreak)

R       = "resultados_experimento"
OUT     = "pdfs_corridas"
ALGOS   = ["MAREA", "MAREA_REL", "MAREABIO", "MAREABIO_REL",
           "PSO", "DE", "ABC", "GWO", "WOA", "FA", "CLONALG"]
DIMS    = [10, 50]
BUDGETS = [5000, 50000, 500000]
N_RUNS  = 51
BLOQUES = [(0, 26), (26, 51)]          # las 51 corridas partidas a lo ancho

INK, INK2 = colors.HexColor("#1a1a1a"), colors.HexColor("#5a5a5a")
AZUL      = colors.HexColor("#2a78d6")
FONDO     = colors.HexColor("#f0f4f9")


def cargar(alg, dim, bud):
    with open(f"{R}/{alg}_D{dim}_FES{bud}.json") as fh:
        d = json.load(fh)
    g = defaultdict(dict)
    for k, v in d.items():
        m = re.match(r"pid(\d+)_run(\d+)$", k)
        g[int(m.group(1))][int(m.group(2))] = float(v)
    return g


def registro_funciones():
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
    if x is None:          return "-"
    if not np.isfinite(x): return "n/f"
    if x == 0:             return "0.000e+00"
    return f"{x:.3e}"


def construir(alg, funcs):
    ss = getSampleStyleSheet()
    h1 = ParagraphStyle("h1", parent=ss["Title"], fontSize=15, textColor=INK,
                        alignment=0, spaceAfter=2)
    sub = ParagraphStyle("sub", parent=ss["Normal"], fontSize=8.5, textColor=INK2,
                         spaceAfter=8)
    cuerpo = ParagraphStyle("b", parent=ss["Normal"], fontSize=8.5, textColor=INK,
                            leading=12, spaceAfter=6)

    doc = SimpleDocTemplate(f"{OUT}/{alg}_corridas.pdf", pagesize=landscape(A4),
                            leftMargin=9*mm, rightMargin=9*mm,
                            topMargin=10*mm, bottomMargin=10*mm,
                            title=f"Corridas individuales — {alg}",
                            author="Benjamin Serrano")
    el = []
    el.append(Paragraph(f"{alg} — corridas individuales", h1))
    el.append(Paragraph("Benchmark MAREA · paper OE1 · 71 funciones CEC2014/2017/2022 · "
                        "51 corridas · 6 celdas (D10/D50 × 5k/50k/500k)", sub))
    el.append(Paragraph(
        "Cada celda de la tabla es el <b>error al optimo</b> f(x) - f* de UNA corrida. "
        "Menor es mejor. Las columnas r0..r50 son los indices de corrida (la semilla de cada "
        "una es su indice), partidos en dos bloques por ancho de pagina: r0-r25 y luego r26-r50. "
        "Dominio [-100,100], poblacion N=80, T=MaxFES/80. Las 12 funciones de CEC2022 solo "
        "llegan a D=20, asi que en las celdas D50 corren en <b>D=20</b>. "
        "Valores redondeados a 4 cifras significativas: para reanalisis exacto usar los JSON "
        "de resultados_experimento/, que van en precision completa.", cuerpo))

    for dim in DIMS:
        for bud in BUDGETS:
            datos = cargar(alg, dim, bud)
            for bi, (r0, r1) in enumerate(BLOQUES):
                el.append(PageBreak())
                el.append(Paragraph(
                    f"{alg} — D{dim}, {bud//1000}k FES — corridas r{r0}–r{r1-1}", h1))
                el.append(Paragraph(
                    f"error al optimo por corrida · bloque {bi+1} de {len(BLOQUES)} "
                    f"· 71 funciones", sub))
                enc = ["func"] + [f"r{r}" for r in range(r0, r1)]
                fil = [enc]
                for pid, fmeta in enumerate(funcs):
                    fila = [fmeta["name"]]
                    runs = datos.get(pid, {})
                    fila += [num(runs.get(r)) for r in range(r0, r1)]
                    fil.append(fila)
                ncol = r1 - r0
                anchos = [15*mm] + [(263.0 - 15) / ncol * mm] * ncol
                t = Table(fil, repeatRows=1, hAlign="LEFT", colWidths=anchos)
                t.setStyle(TableStyle([
                    ("FONT", (0,0), (-1,0), "Helvetica-Bold", 5.2),
                    ("FONT", (0,1), (-1,-1), "Helvetica", 4.9),
                    ("TEXTCOLOR", (0,0), (-1,0), colors.white),
                    ("BACKGROUND", (0,0), (-1,0), AZUL),
                    ("ROWBACKGROUNDS", (0,1), (-1,-1), [colors.white, FONDO]),
                    ("ALIGN", (1,0), (-1,-1), "RIGHT"),
                    ("ALIGN", (0,0), (0,-1), "LEFT"),
                    ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
                    ("GRID", (0,0), (-1,-1), 0.15, colors.HexColor("#ccd4dc")),
                    ("TOPPADDING", (0,0), (-1,-1), 0.9),
                    ("BOTTOMPADDING", (0,0), (-1,-1), 0.9),
                    ("LEFTPADDING", (0,0), (-1,-1), 1.2),
                    ("RIGHTPADDING", (0,0), (-1,-1), 1.2),
                ]))
                el.append(t)
    doc.build(el)


def main():
    os.makedirs(OUT, exist_ok=True)
    funcs = registro_funciones()
    for alg in ALGOS:
        construir(alg, funcs)
        p = f"{OUT}/{alg}_corridas.pdf"
        print(f"  {p}  {os.path.getsize(p)/1024:.0f} KB")
    print("listo")


if __name__ == "__main__":
    main()
