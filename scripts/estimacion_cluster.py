#!/usr/bin/env python3
"""
Practica 3 - Ingenieria de los Computadores (2026-27)
Estima el tiempo que tardarian campo_mpi y campo_mpi_func en un cluster de
varios equipos a partir de las medidas hechas en UNA sola maquina.

En una sola maquina los mensajes MPI se copian por memoria compartida, asi
que la comunicacion es casi gratis. En un cluster viajan por la red. Se
suma a cada medida local el coste de red de los mensajes que envia el
programa, con el modelo clasico  t = latencia + bytes / ancho_de_banda.

Supuestos (se pueden cambiar con los argumentos):
  - Un proceso MPI por equipo (caso peor: todos los mensajes van por la red).
  - Ethernet a 1 Gbit/s: ~117 MB/s utiles y ~100 us de latencia por mensaje.
  - El raiz recibe por un unico enlace, asi que los trozos que le llegan se
    transmiten uno detras de otro (no se solapan).

Estrategia A (por proceso raiz):
  2 x Bcast pequenos (dims y cargas) + 2 x Gatherv de la malla (V y E):
  T_red = 2 * ceil(log2 np) * lat  +  2 * [ (np-1) * lat + (np-1)/np * M / BW ]
Estrategia B:
  Bcast igual; Gatherv de V (grupo de nV) y de E (grupo de nE) a la vez en
  raices distintas, y despues Send de E completo al rank 0:
  T_red = Bcast + max(Gather_V, Gather_E) + (lat + M / BW)
donde M = ancho * alto * 8 bytes (una malla de doubles).

Uso:
  python3 scripts/estimacion_cluster.py                 # 1 GbE
  python3 scripts/estimacion_cluster.py --bw 1170e6     # 10 GbE
"""

import argparse
import csv
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from graficas import (agregar, leer, t_secuencial, base_caso, guardar, plt,  # noqa: E402
                      SERIES, TINTA_2, nombre)


def t_red(programa, np_, ancho, alto, lat, bw):
    malla = ancho * alto * 8.0
    bcast = 2 * math.ceil(math.log2(np_)) * lat if np_ > 1 else 0.0

    def gather(p):
        return (p - 1) * lat + (p - 1) / p * malla / bw if p > 1 else 0.0

    if programa == "campo_mpi":
        return bcast + 2 * gather(np_)
    nV = np_ // 2
    nE = np_ - nV
    return bcast + max(gather(nV), gather(nE)) + lat + malla / bw


def main():
    ap = argparse.ArgumentParser(description="Estimacion del tiempo en un cluster.")
    ap.add_argument("csv", nargs="*", default=["resultados/medidas.csv"])
    ap.add_argument("--bw", type=float, default=117e6, help="ancho de banda util (bytes/s)")
    ap.add_argument("--lat", type=float, default=100e-6, help="latencia por mensaje (s)")
    ap.add_argument("--salida", default="resultados/graficas")
    ap.add_argument("--secuencial", default="campo")
    args = ap.parse_args()

    filas = leer(args.csv)
    dims = {r["caso"]: (int(r["ancho"]), int(r["alto"])) for r in filas}
    agg = agregar(filas)
    tseq = t_secuencial(agg, args.secuencial)
    os.makedirs(args.salida, exist_ok=True)

    progs = ["campo_mpi", "campo_mpi_func"]
    est = {}  # (programa, caso, np) -> dict
    for (p, c, n, _), v in agg.items():
        if p not in progs:
            continue
        red = t_red(p, n, *dims[c], args.lat, args.bw)
        tot = v["media"] + red
        est[(p, c, n)] = {"local": v["media"], "calculo": v["calculo"], "red": red,
                          "total": tot, "s_local": tseq[c] / v["media"],
                          "s_est": tseq[c] / tot}

    ruta = os.path.join(args.salida, "estimacion_cluster.csv")
    with open(ruta, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["programa", "caso", "np", "t_local", "t_red_estimado", "t_cluster_estimado",
                    "pct_red", "speedup_local", "speedup_cluster", "eficiencia_cluster"])
        for (p, c, n), e in sorted(est.items(), key=lambda kv: (kv[0][1], kv[0][0], kv[0][2])):
            w.writerow([p, c, n, f"{e['local']:.4f}", f"{e['red']:.4f}", f"{e['total']:.4f}",
                        f"{100 * e['red'] / e['total']:.1f}", f"{e['s_local']:.3f}",
                        f"{e['s_est']:.3f}", f"{e['s_est'] / n:.3f}"])
    print("Generando:\n  ", ruta)

    gbit = args.bw * 8 / 1e9 / 0.936  # 117 MB/s utiles ~ 1 Gbit/s de linea
    for caso in sorted(tseq):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
        nps = []
        for k, p in enumerate(progs):
            pts = sorted((n, e) for (pp, c, n), e in est.items() if pp == p and c == caso)
            if not pts:
                continue
            x = [n for n, _ in pts]
            nps = sorted(set(nps) | set(x))
            ax1.plot(x, [e["s_local"] for _, e in pts], marker="o", color=SERIES[k],
                     label=f"{nombre(p)} · medido (1 máquina)")
            ax1.plot(x, [e["s_est"] for _, e in pts], marker="s", linestyle="--",
                     color=SERIES[k], label=f"{nombre(p)} · estimado (clúster)")
            ax2.plot(x, [100 * e["red"] / e["total"] for _, e in pts], marker="s",
                     linestyle="--", color=SERIES[k], label=nombre(p))
        ax1.plot(nps, nps, linestyle=":", linewidth=1.2, color=TINTA_2, label="ideal")
        ax1.set_xlabel("Procesos MPI (np)")
        ax1.set_ylabel("Speed-up")
        ax1.set_ylim(bottom=0)
        ax1.set_title("Speed-up medido y estimado")
        ax1.legend(fontsize=8)
        ax2.set_xlabel("Procesos MPI (np)")
        ax2.set_ylabel("% del tiempo total en la red")
        ax2.set_ylim(0, 100)
        ax2.set_title("Peso estimado de la comunicación")
        ax2.legend(fontsize=8)
        for ax in (ax1, ax2):
            ax.set_xticks(nps)
        fig.suptitle(f"{base_caso(caso)}: estimación con red de ~{gbit:.0f} Gbit/s "
                     f"({args.bw / 1e6:.0f} MB/s, {args.lat * 1e6:.0f} µs), un proceso por equipo",
                     x=0.01, ha="left", fontsize=10, color=TINTA_2)
        guardar(fig, args.salida, f"cluster_{base_caso(caso)}.png")


if __name__ == "__main__":
    main()
