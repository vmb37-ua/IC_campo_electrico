#!/usr/bin/env python3
"""
Practica 3 - Ingenieria de los Computadores (2026-27)
Lee el CSV de scripts/medir.py, calcula medias, speed-up y eficiencia, y
genera las graficas para la memoria.

  python3 scripts/graficas.py                       # resultados/medidas.csv
  python3 scripts/graficas.py a.csv b.csv --salida resultados/graficas

Requiere matplotlib (pip install matplotlib).

Speed-up   S = T_secuencial / T_paralelo   (T = media de t_total del caso)
Eficiencia E = S / np

Sirve para otros programas (p. ej. el kNN) si su CSV tiene al menos las
columnas: programa, caso, np, nodos, t_calculo, t_comunicacion, t_total.
Las columnas ancho, alto y cargas son opcionales (ordenan los casos por talla).
"""

import argparse
import csv
import os
import statistics
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- Estilo -------------------------------------------------------------------
SUPERFICIE = "#fcfcfb"
TINTA = "#0b0b0b"
TINTA_2 = "#52514e"
REJILLA = "#e4e3df"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
          "#e87ba4", "#008300", "#4a3aa7", "#e34948"]   # orden fijo
CALCULO, COMUNICACION = SERIES[0], SERIES[1]

NOMBRES = {"campo": "Secuencial", "campo_mpi": "A · datos",
           "campo_mpi_func": "B · funcional"}

plt.rcParams.update({
    "figure.facecolor": SUPERFICIE, "axes.facecolor": SUPERFICIE,
    "savefig.facecolor": SUPERFICIE, "figure.dpi": 110, "savefig.dpi": 200,
    "font.size": 10, "text.color": TINTA, "axes.labelcolor": TINTA_2,
    "axes.edgecolor": TINTA_2, "axes.linewidth": 0.8,
    "xtick.color": TINTA_2, "ytick.color": TINTA_2,
    "axes.grid": True, "axes.axisbelow": True, "grid.color": REJILLA, "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "legend.frameon": False, "lines.linewidth": 2, "lines.markersize": 7,
})


def nombre(p):
    return NOMBRES.get(p, p)


def leer(ficheros):
    filas = []
    for fich in ficheros:
        with open(fich, newline="") as f:
            filas += list(csv.DictReader(f))
    return filas


def flotante(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def agregar(filas):
    """Agrupa por (programa, caso, np, nodos) y calcula medias y desviaciones."""
    grupos = defaultdict(list)
    for r in filas:
        clave = (r["programa"], r["caso"], int(r["np"]), int(r.get("nodos") or 1))
        grupos[clave].append(r)
    res = {}
    for clave, rs in grupos.items():
        tot = [flotante(r["t_total"]) for r in rs]
        cal = [flotante(r.get("t_calculo")) or 0.0 for r in rs]
        com = [flotante(r.get("t_comunicacion")) or 0.0 for r in rs]
        r0 = rs[0]
        talla = None
        if r0.get("ancho") and r0.get("alto") and r0.get("cargas"):
            talla = int(r0["ancho"]) * int(r0["alto"]) * int(r0["cargas"])
        res[clave] = {
            "reps": len(tot), "media": statistics.mean(tot),
            "desv": statistics.stdev(tot) if len(tot) > 1 else 0.0,
            "calculo": statistics.mean(cal), "comunicacion": statistics.mean(com),
            "talla": talla,
            "dims": f"{r0.get('ancho')}x{r0.get('alto')}, {r0.get('cargas')} cargas"
                    if talla else "",
        }
    return res


def series_paralelas(agg, secuencial):
    """Claves (programa, nodos) de las versiones paralelas, en orden fijo."""
    claves = sorted({(p, nd) for (p, _, _, nd) in agg if p != secuencial})
    varios_nodos = len({nd for _, nd in claves}) > 1
    etiquetas = {k: nombre(k[0]) + (f" · {k[1]} nodos" if varios_nodos else "")
                 for k in claves}
    colores = {k: SERIES[i % len(SERIES)] for i, k in enumerate(claves)}
    return claves, etiquetas, colores


def t_secuencial(agg, secuencial):
    return {caso: v["media"] for (p, caso, _, _), v in agg.items() if p == secuencial}


def guardar(fig, carpeta, nombre_fich):
    fig.tight_layout()
    ruta = os.path.join(carpeta, nombre_fich)
    fig.savefig(ruta)
    plt.close(fig)
    print("  ", ruta)


def base_caso(caso):
    return os.path.splitext(os.path.basename(caso))[0]


# --- Graficas ---------------------------------------------------------------------
def grafica_speedup_eficiencia(agg, tseq, caso, claves, etiq, col, carpeta):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    nps_todos = set()
    for k in claves:
        pts = sorted((np_, v["media"]) for (p, c, np_, nd), v in agg.items()
                     if (p, nd) == k and c == caso)
        if not pts:
            continue
        nps = [n for n, _ in pts]
        nps_todos.update(nps)
        s = [tseq[caso] / t for _, t in pts]
        e = [si / n for si, n in zip(s, nps)]
        ax1.plot(nps, s, marker="o", color=col[k], label=etiq[k])
        ax2.plot(nps, e, marker="o", color=col[k], label=etiq[k])
        ax1.annotate(f"{s[-1]:.2f}", (nps[-1], s[-1]), xytext=(6, 0),
                     textcoords="offset points", va="center", color=TINTA, fontsize=9)
    if not nps_todos:
        plt.close(fig)
        return
    lim = [1, max(nps_todos)]
    ax1.plot(lim, lim, linestyle="--", linewidth=1.2, color=TINTA_2, label="Ideal (S = np)")
    ax2.axhline(1.0, linestyle="--", linewidth=1.2, color=TINTA_2, label="Ideal (E = 1)")
    for ax, titulo, ylab in [(ax1, "Speed-up", "S = T_sec / T_par"),
                             (ax2, "Eficiencia", "E = S / np")]:
        ax.set_title(titulo)
        ax.set_xlabel("Número de procesos (np)")
        ax.set_ylabel(ylab)
        ax.set_xticks(sorted(nps_todos | {1}))
        ax.set_ylim(bottom=0)
        ax.legend(loc="best")
    fig.suptitle(f"{base_caso(caso)}   (T secuencial = {tseq[caso]:.2f} s)",
                 x=0.01, ha="left", fontsize=10, color=TINTA_2)
    guardar(fig, carpeta, f"speedup_{base_caso(caso)}.png")


def grafica_comunicacion(agg, caso, claves, etiq, carpeta):
    """Barras apiladas calculo + comunicacion por np (preguntas 3 y 9)."""
    claves = [k for k in claves if any((p, nd) == k and c == caso for (p, c, _, nd) in agg)]
    if not claves:
        return
    fig, axes = plt.subplots(1, len(claves), figsize=(5 * len(claves), 3.8),
                             sharey=True, squeeze=False)
    for ax, k in zip(axes[0], claves):
        pts = sorted((np_, v) for (p, c, np_, nd), v in agg.items()
                     if (p, nd) == k and c == caso)
        etiquetas_x = [str(n) for n, _ in pts]
        cal = [v["calculo"] for _, v in pts]
        com = [max(v["media"] - v["calculo"], 0.0) for _, v in pts]
        x = range(len(pts))
        ax.bar(x, cal, width=0.6, color=CALCULO, label="Cálculo",
               edgecolor=SUPERFICIE, linewidth=1)
        ax.bar(x, com, width=0.6, bottom=cal, color=COMUNICACION,
               label="Comunicación y espera", edgecolor=SUPERFICIE, linewidth=1)
        for xi, c1, c2 in zip(x, cal, com):
            tot = c1 + c2
            ax.annotate(f"{100 * c2 / tot:.1f} %", (xi, tot), xytext=(0, 3),
                        textcoords="offset points", ha="center", fontsize=8, color=TINTA)
        ax.set_xticks(list(x), etiquetas_x)
        ax.set_xlabel("Número de procesos (np)")
        ax.set_title(etiq[k])
        ax.grid(axis="x", visible=False)
    axes[0][0].set_ylabel("Tiempo (s)")
    axes[0][0].legend(loc="upper right")
    fig.suptitle(f"{base_caso(caso)}: reparto del tiempo (el % es la parte de comunicación)",
                 x=0.01, ha="left", fontsize=10, color=TINTA_2)
    guardar(fig, carpeta, f"comunicacion_{base_caso(caso)}.png")


def grafica_tallas(agg, tseq, serie, casos, claves, etiq, col, carpeta):
    """Speed-up y eficiencia frente a la talla, un panel por np (small multiples)."""
    casos = [c for c in casos if c in tseq]
    if len(casos) < 2:
        return
    tallas = {c: next(v["talla"] for (p, cc, _, _), v in agg.items() if cc == c) for c in casos}
    casos.sort(key=lambda c: (tallas[c] or 0, c))
    nps = sorted({np_ for (p, c, np_, nd) in agg
                  if c in casos and (p, nd) in claves and np_ > 1})
    if not nps:
        return
    fig, axes = plt.subplots(2, len(nps), figsize=(3.4 * len(nps) + 1, 6.4),
                             sharex=True, sharey="row", squeeze=False)
    xs = list(range(len(casos)))
    etiquetas_x = [base_caso(c).split("_", 1)[-1] for c in casos]
    for j, n in enumerate(nps):
        for k in claves:
            pts = [(i, agg.get((k[0], c, n, k[1]))) for i, c in enumerate(casos)]
            pts = [(i, v) for i, v in pts if v]
            if not pts:
                continue
            s = [tseq[casos[i]] / v["media"] for i, v in pts]
            axes[0][j].plot([i for i, _ in pts], s, marker="o", color=col[k], label=etiq[k])
            axes[1][j].plot([i for i, _ in pts], [si / n for si in s], marker="o",
                            color=col[k], label=etiq[k])
        axes[0][j].axhline(n, linestyle="--", linewidth=1.2, color=TINTA_2)
        axes[1][j].axhline(1.0, linestyle="--", linewidth=1.2, color=TINTA_2)
        axes[0][j].set_title(f"np = {n}")
        axes[1][j].set_xticks(xs, etiquetas_x, rotation=45 if len(casos) > 4 else 0)
        axes[1][j].set_xlabel(f"Talla ({serie})")
    axes[0][0].set_ylabel("Speed-up")
    axes[1][0].set_ylabel("Eficiencia")
    axes[0][0].set_ylim(bottom=0)
    axes[1][0].set_ylim(0, 1.1)
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper right", ncol=len(l))
    fig.suptitle(f"Escalado con el tamaño del problema ({serie}); línea discontinua = ideal",
                 x=0.01, ha="left", fontsize=10, color=TINTA_2)
    guardar(fig, carpeta, f"tallas_{serie}.png")


def escribir_resumen(agg, tseq, secuencial, ruta):
    with open(ruta, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["programa", "caso", "np", "nodos", "reps", "t_total_media", "t_total_desv",
                    "t_calculo", "t_comunicacion", "pct_comunicacion", "speedup", "eficiencia"])
        for (p, c, n, nd), v in sorted(agg.items(), key=lambda kv: (kv[0][1], kv[0][0], kv[0][2], kv[0][3])):
            s = tseq[c] / v["media"] if c in tseq else None
            com = max(v["media"] - v["calculo"], 0.0) if p != secuencial else 0.0
            w.writerow([p, c, n, nd, v["reps"], f"{v['media']:.4f}", f"{v['desv']:.4f}",
                        f"{v['calculo']:.4f}", f"{com:.4f}", f"{100 * com / v['media']:.2f}",
                        f"{s:.3f}" if s else "", f"{s / n:.3f}" if s else ""])
    print("  ", ruta)


def main():
    ap = argparse.ArgumentParser(description="Graficas de speed-up y eficiencia.")
    ap.add_argument("csv", nargs="*", default=["resultados/medidas.csv"])
    ap.add_argument("--salida", default="resultados/graficas", help="carpeta de las imagenes")
    ap.add_argument("--secuencial", default="campo", help="programa de referencia (T_sec)")
    args = ap.parse_args()

    agg = agregar(leer(args.csv))
    tseq = t_secuencial(agg, args.secuencial)
    if not tseq:
        raise SystemExit(f"No hay medidas del programa secuencial '{args.secuencial}'.")
    claves, etiq, col = series_paralelas(agg, args.secuencial)
    os.makedirs(args.salida, exist_ok=True)

    print("Generando:")
    escribir_resumen(agg, tseq, args.secuencial, os.path.join(args.salida, "resumen.csv"))
    casos = sorted(tseq)
    for caso in casos:
        grafica_speedup_eficiencia(agg, tseq, caso, claves, etiq, col, args.salida)
        grafica_comunicacion(agg, caso, claves, etiq, args.salida)

    # Series de tallas: casos agrupados por el prefijo del nombre (espacio_*, cargas_*)
    series = defaultdict(list)
    for c in casos:
        b = base_caso(c)
        if "_" in b:
            series[b.split("_", 1)[0]].append(c)
    for serie, cs in series.items():
        grafica_tallas(agg, tseq, serie, cs, claves, etiq, col, args.salida)


if __name__ == "__main__":
    main()
