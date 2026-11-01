#!/usr/bin/env python3
"""
Practica 3 - Ingenieria de los Computadores (2026-27)
Lanza las versiones secuencial y MPI de campo con distintos casos y numeros
de procesos, repite cada medida y guarda los tiempos en un CSV.

Ejemplos (desde la raiz del repositorio):

  # En local, caso de referencia, 1-8 procesos, 3 repeticiones
  python3 scripts/medir.py --np 2 4 8 --reps 3

  # Serie A completa en el cluster
  python3 scripts/medir.py --casos serieA --np 2 4 6 8 12 --hostfile hosts

  # Probar con mas procesos que nucleos (los tiempos no son validos)
  python3 scripts/medir.py --casos ejemplo.txt --np 16 --oversubscribe --reps 1

El CSV se abre en modo anadir: se pueden juntar medidas de varias sesiones o
de distintos hostfiles en el mismo fichero (la columna 'nodos' los distingue).
"""

import argparse
import csv
import datetime
import glob
import os
import re
import socket
import subprocess
import sys

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

COLUMNAS = ["fecha", "etiqueta", "maquina", "programa", "caso", "ancho", "alto",
            "cargas", "np", "nodos", "rep", "t_calculo", "t_calculo_v",
            "t_calculo_e", "t_comunicacion", "t_total"]

# Atajos para --casos
SERIES = {
    "serieA": sorted(glob.glob(os.path.join(RAIZ, "casos", "espacio_*.txt")),
                     key=lambda f: int(re.search(r"(\d+)", os.path.basename(f)).group(1))),
    "serieB": sorted(glob.glob(os.path.join(RAIZ, "casos", "cargas_*.txt")),
                     key=lambda f: int(re.search(r"(\d+)", os.path.basename(f)).group(1))),
}


def numero(patron, texto):
    m = re.search(patron, texto)
    return float(m.group(1)) if m else None


def interpretar(programa, salida):
    """Extrae dimensiones y tiempos de la salida de campo / campo_mpi / campo_mpi_func."""
    m = re.search(r"Espacio:\s*(\d+)\s*x\s*(\d+)\s*um\s+Cargas:\s*(\d+)", salida)
    if not m:
        return None
    fila = {"ancho": int(m.group(1)), "alto": int(m.group(2)), "cargas": int(m.group(3))}

    if programa == "campo":
        # Secuencial: no hay comunicacion
        tv = numero(r"Tiempo de calculo del potencial:\s*([\d.]+)", salida)
        te = numero(r"Tiempo de calculo del campo:\s*([\d.]+)", salida)
        tt = numero(r"Tiempo de calculo total:\s*([\d.]+)", salida)
        fila.update(t_calculo=tt, t_calculo_v=tv, t_calculo_e=te,
                    t_comunicacion=0.0, t_total=tt)
    elif programa == "campo_mpi_func":
        tv = numero(r"Tiempo de calculo grupo V \(max\):\s*([\d.]+)", salida)
        te = numero(r"Tiempo de calculo grupo E \(max\):\s*([\d.]+)", salida)
        fila.update(t_calculo=max(tv, te), t_calculo_v=tv, t_calculo_e=te,
                    t_comunicacion=numero(r"Tiempo de comunicacion:\s*([\d.]+)", salida),
                    t_total=numero(r"Tiempo total:\s*([\d.]+)", salida))
    else:  # campo_mpi
        fila.update(t_calculo=numero(r"Tiempo de calculo \(max\):\s*([\d.]+)", salida),
                    t_calculo_v="", t_calculo_e="",
                    t_comunicacion=numero(r"Tiempo de comunicacion:\s*([\d.]+)", salida),
                    t_total=numero(r"Tiempo total:\s*([\d.]+)", salida))
    return fila


def contar_nodos(hostfile):
    if not hostfile:
        return 1
    with open(hostfile) as f:
        hosts = {l.split()[0] for l in f if l.strip() and not l.strip().startswith("#")}
    return len(hosts)


def comando(programa, caso, np_, args):
    binario = os.path.join(".", programa)
    if programa == "campo":
        return [binario, caso]
    cmd = ["mpirun"]
    if args.hostfile:
        cmd += ["-mca", "plm_rsh_no_tree_spawn", "1", "-hostfile", args.hostfile]
    if args.oversubscribe:
        cmd += ["--oversubscribe"]
    return cmd + ["-np", str(np_), binario, caso]


def main():
    ap = argparse.ArgumentParser(description="Mide campo, campo_mpi y campo_mpi_func.")
    ap.add_argument("--programas", nargs="+", default=["campo", "campo_mpi", "campo_mpi_func"],
                    help="binarios a medir (por defecto los tres)")
    ap.add_argument("--casos", nargs="+", default=["espacio.txt"],
                    help="ficheros de entrada, o los atajos serieA / serieB")
    ap.add_argument("--np", nargs="+", type=int, default=[2, 4, 8],
                    help="numeros de procesos MPI (el secuencial siempre usa 1)")
    ap.add_argument("--reps", type=int, default=3, help="repeticiones de cada medida")
    ap.add_argument("--hostfile", help="fichero de hosts para ejecutar en el cluster")
    ap.add_argument("--oversubscribe", action="store_true",
                    help="permite mas procesos que nucleos (solo para pruebas)")
    ap.add_argument("--etiqueta", default="", help="texto libre para identificar la tanda")
    ap.add_argument("--salida", default="resultados/medidas.csv", help="CSV de salida")
    args = ap.parse_args()

    os.chdir(RAIZ)
    casos = []
    for c in args.casos:
        casos += [os.path.relpath(f, RAIZ) for f in SERIES[c]] if c in SERIES else [c]

    # Compila lo necesario
    if subprocess.call(["make", "-s"] + args.programas) != 0:
        sys.exit("Error al compilar")

    nodos = contar_nodos(args.hostfile)
    maquina = socket.gethostname()
    os.makedirs(os.path.dirname(args.salida) or ".", exist_ok=True)
    nuevo = not os.path.exists(args.salida)

    # Lista de ejecuciones: el secuencial una vez por caso; B necesita np >= 2
    trabajos = []
    for caso in casos:
        for prog in args.programas:
            for np_ in ([1] if prog == "campo" else args.np):
                if prog == "campo_mpi_func" and np_ < 2:
                    continue
                for rep in range(1, args.reps + 1):
                    trabajos.append((prog, caso, np_, rep))

    with open(args.salida, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNAS)
        if nuevo:
            w.writeheader()
        for i, (prog, caso, np_, rep) in enumerate(trabajos, 1):
            cmd = comando(prog, caso, np_, args)
            print(f"[{i}/{len(trabajos)}] {prog:15s} {caso:25s} np={np_:<3d} rep={rep} ... ",
                  end="", flush=True)
            r = subprocess.run(cmd, capture_output=True, text=True)
            fila = interpretar(prog, r.stdout) if r.returncode == 0 else None
            if fila is None or fila.get("t_total") is None:
                print("ERROR")
                print(r.stderr.strip()[-500:], file=sys.stderr)
                continue
            fila.update(fecha=datetime.datetime.now().isoformat(timespec="seconds"),
                        etiqueta=args.etiqueta, maquina=maquina, programa=prog, caso=caso,
                        np=np_, nodos=nodos, rep=rep)
            w.writerow(fila)
            f.flush()  # no perder medidas si se interrumpe
            print(f"{fila['t_total']:.3f} s")

    print(f"\nMedidas guardadas en {args.salida}")


if __name__ == "__main__":
    main()
