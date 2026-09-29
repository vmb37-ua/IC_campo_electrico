#!/usr/bin/env python3
"""Genera los ficheros de entrada de los casos de prueba.

Serie A: se varia el tamano del espacio con 100 cargas fijas.
Serie B: se varia el numero de cargas con un espacio de 4000x4000.
Las posiciones son deterministas (semilla fija) para que las medidas
sean reproducibles.
"""
import random

def generar(nombre, lado, n_cargas):
    random.seed(1234)
    lineas = ["%dx%d" % (lado, lado)]
    for _ in range(n_cargas):
        q = random.uniform(1.0, 10.0) * random.choice((1, -1))
        x = random.uniform(0, lado - 1)
        y = random.uniform(0, lado - 1)
        lineas.append("(%.2f, %.1f, %.1f)" % (q, x, y))
    with open(nombre, "w") as f:
        f.write("\n".join(lineas) + "\n")
    print("%-24s %5dx%-5d %4d cargas  %11d evaluaciones"
          % (nombre, lado, lado, n_cargas, lado * lado * n_cargas))

# Serie A - escalado del tamano del espacio (100 cargas)
for lado in (2000, 3000, 5000, 7000, 9000, 11000):
    generar("casos/espacio_%d.txt" % lado, lado, 100)

# Serie B - escalado del numero de cargas (4000x4000)
for n in (50, 100, 200, 400):
    generar("casos/cargas_%d.txt" % n, 4000, n)
