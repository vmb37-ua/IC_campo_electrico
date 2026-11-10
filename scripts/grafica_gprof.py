#!/usr/bin/env python3
"""Grafica del perfil plano de gprof (porcentaje de tiempo propio por funcion).

Uso: python3 scripts/grafica_gprof.py resultados/gprof/perfil_O0.txt resultados/gprof/perfil_O2_noinline.txt
"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from graficas import plt, SERIES, TINTA_2, guardar  # noqa: E402


def perfil_plano(ruta, n=5):
    filas = []
    with open(ruta) as f:
        texto = f.read().split("Call graph")[0]
    for linea in texto.splitlines():
        m = re.match(r"\s*([\d.]+)\s+[\d.]+\s+([\d.]+)\s+(?:\d+\s+[\d.]+\s+[\d.]+\s+)?(.+)$", linea)
        if m:
            nombre = re.sub(r"\(.*\)", "()", m.group(3)).replace("std::vector<Carga, std::allocator<Carga> >", "vector<Carga>")
            filas.append((float(m.group(1)), nombre.strip()))
    return filas[:n]


def main():
    titulos = {"perfil_O0.txt": "-O0 -g -pg", "perfil_O2_noinline.txt": "-O2 -fno-inline -g -pg",
               "perfil_O2.txt": "-O2 -g -pg"}
    rutas = sys.argv[1:]
    fig, axes = plt.subplots(len(rutas), 1, figsize=(9, 2.3 * len(rutas)))
    for ax, ruta in zip(axes, rutas):
        datos = perfil_plano(ruta)[::-1]
        ax.barh([d[1] for d in datos], [d[0] for d in datos], color=SERIES[0], height=0.6)
        for i, (p, _) in enumerate(datos):
            ax.text(p + 1, i, f"{p:.1f} %", va="center", fontsize=9, color=TINTA_2)
        ax.set_xlim(0, 110)
        ax.set_title(titulos.get(os.path.basename(ruta), ruta))
        ax.grid(axis="y", visible=False)
    axes[-1].set_xlabel("% del tiempo propio (self) muestreado por gprof")
    guardar(fig, "resultados/graficas", "gprof_campo.png")


if __name__ == "__main__":
    main()
