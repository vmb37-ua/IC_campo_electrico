# IC_campo_electrico

Cálculo secuencial del **potencial eléctrico** y del **campo eléctrico** en una región
plana del espacio. Dadas *n* cargas puntuales colocadas en esa región, el programa
evalúa ambas magnitudes en cada punto de una malla y las representa en el terminal,
una debajo de otra, con colores ANSI de 24 bits:

- **Potencial V:** rojo = potencial alto, azul = potencial bajo, verde oscuro ≈ potencial nulo.
- **Módulo del campo |E|:** amarillo = campo intenso, azul oscuro = campo débil.

Al terminar imprime el tiempo de cálculo de cada magnitud y el total, en segundos.

Es la aplicación secuencial de referencia de la Práctica 2 de Ingeniería de los
Computadores: **no usa hebras, procesos ni ningún tipo de paralelismo explícito**.

## Modelo físico

El potencial electrostático creado por *n* cargas puntuales en un punto **r** es la
suma escalar de las contribuciones individuales (principio de superposición):

```
             n        q_i
V(x, y) = k · Σ  ─────────────        r_i = √( (x − x_i)² + (y − y_i)² )
            i=1       r_i
```

- `k = 8.9875517923 · 10⁹ N·m²/C²` — constante de Coulomb (vacío).
- El potencial es una magnitud **escalar**, por lo que las contribuciones se suman
  directamente sin descomponer en componentes.

El campo eléctrico es una magnitud **vectorial**: cada carga aporta un vector que
apunta en la dirección radial, así que se acumulan las componentes por separado
y al final se calcula el módulo:

```
              n    q_i · (x − x_i)                 n    q_i · (y − y_i)
E_x(x, y) = k · Σ  ───────────────     E_y(x, y) = k · Σ  ───────────────
             i=1       r_i³                       i=1       r_i³

|E| = √( E_x² + E_y² )
```

El campo se calcula **directamente a partir de las cargas**, no como −∇V, de modo
que el cálculo del campo y el del potencial son independientes entre sí: los dos
solo dependen de la lectura de las cargas.

- Se trabaja en 2D: la malla es un corte plano de la región, y las cargas se
  consideran contenidas en ese mismo plano.

### Unidades

| Magnitud | Unidad en el fichero | Conversión interna |
|---|---|---|
| Tamaño del espacio | micrómetros (µm) | 1 celda de la malla = 1 µm |
| Posición de las cargas `(x, y)` | micrómetros (µm) | × 10⁻⁶ → metros |
| Valor de la carga `q` | nanoculombios (nC) | × 10⁻⁹ → culombios |
| Potencial calculado `V` | **voltios (V)** | — |
| Campo calculado `|E|` | **voltios por metro (V/m)** | — |

Todo el cálculo se hace en `double`.

### Singularidad en las cargas

`V → ∞` y `|E| → ∞` cuando `r → 0`. Para evitar la división por cero en las celdas que caen justo
sobre una carga se aplica un radio mínimo de `R_MIN = 0.5 µm` (la mitad del paso de la
malla): si `r < R_MIN`, se usa `R_MIN`. Es el equivalente al *softening* habitual en
simulaciones de N cuerpos.

### Escala de color

Potencial y campo crecen sin límite junto a las cargas, así que una escala lineal
saturaría la imagen. Se normaliza con

```
t = tanh( D / D_ref )        D_ref = media de |D| sobre toda la malla   (D = V o |E|)
```

- **Potencial** (escala divergente): `t ∈ (−1, 1)`; `t = +1 → rojo`, `t = −1 → azul`,
  `t = 0 → verde oscuro`.
- **Campo** (escala secuencial, `|E| ≥ 0`): `t ∈ [0, 1)`; `t = 0 → azul oscuro`,
  `t → 1 → amarillo`.

## Formato del fichero de entrada

Primera línea: dimensiones de la región en µm, `AnchoxAlto`.
Resto de líneas: una carga por línea, `(carga_en_nC, x_en_µm, y_en_µm)`.

```
100x100
(5, 34, 5)
(-4, 20, 60)
```

La coordenada `(0, 0)` es la **esquina superior izquierda**; `x` crece hacia la derecha
e `y` hacia abajo. Las líneas que no encajen en ninguno de los dos formatos se ignoran,
así que se pueden usar como comentarios.

## Compilación y ejecución

```sh
make                 # compila (g++ -O2 -std=c++11)
make run             # ejecuta el caso de prueba de referencia (espacio.txt)
make clean           # borra el binario
```

Para ejecutar otro caso:

```sh
./campo casos/espacio_9000.txt
make run CASO=casos/cargas_400.txt
```

Para el estudio de opciones de compilación basta con redefinir `CXXFLAGS`:

```sh
make clean && make CXXFLAGS="-O3 -march=native -ffast-math"
```

## Coste computacional

Hay dos núcleos de cálculo, `calcular_potencial()` y `calcular_campo()`. Cada uno es
un triple bucle de coste

```
O(ancho · alto · n_cargas)
```

El del potencial hace una raíz cuadrada y una división por evaluación. El del campo
hace además el cubo de `r` y acumula dos componentes, por lo que es algo más caro.
Los dos núcleos solo leen las cargas y escriben cada uno en su propia malla, así que
son **independientes entre sí**. La memoria ocupada por las dos mallas (`V` y `|E|`)
es `2 · ancho · alto · 8 bytes`.

### Casos de prueba incluidos

Se regeneran con `python3 generar_casos.py` (posiciones deterministas, semilla fija).
Tiempos medidos con `g++ -O2` en un único núcleo, en una sola ejecución:

**Serie A — escalado del tamaño del espacio (100 cargas fijas)**

| Fichero | Espacio (µm) | Evaluaciones por magnitud | Potencial | Campo | Total | Memoria de las mallas |
|---|---|---|---|---|---|---|
| `casos/espacio_2000.txt`  | 2000 × 2000   | 4,0 · 10⁸   | 0,50 s  | 0,57 s  | 1,07 s  | 64 MB |
| `casos/espacio_3000.txt`  | 3000 × 3000   | 9,0 · 10⁸   | 1,04 s  | 1,28 s  | 2,32 s  | 144 MB |
| `casos/espacio_5000.txt`  | 5000 × 5000   | 2,5 · 10⁹   | 2,90 s  | 3,56 s  | 6,46 s  | 400 MB |
| `casos/espacio_7000.txt`  | 7000 × 7000   | 4,9 · 10⁹   | 5,67 s  | 6,99 s  | 12,66 s | 784 MB |
| `casos/espacio_9000.txt`  | 9000 × 9000   | 8,1 · 10⁹   | 9,39 s  | 11,58 s | 20,96 s | 1,30 GB |
| `casos/espacio_11000.txt` | 11000 × 11000 | 1,21 · 10¹⁰ | 14,03 s | 17,24 s | 31,27 s | 1,94 GB |

**Serie B — escalado del número de cargas (espacio fijo de 4000 × 4000 µm)**

| Fichero | Cargas | Evaluaciones por magnitud | Potencial | Campo | Total |
|---|---|---|---|---|---|
| `casos/cargas_50.txt`  | 50  | 8,0 · 10⁸ | 0,93 s | 1,15 s | 2,08 s |
| `casos/cargas_100.txt` | 100 | 1,6 · 10⁹ | 1,86 s | 2,28 s | 4,14 s |
| `casos/cargas_200.txt` | 200 | 3,2 · 10⁹ | 3,71 s | 4,53 s | 8,24 s |
| `casos/cargas_400.txt` | 400 | 6,4 · 10⁹ | 7,41 s | 9,05 s | 16,46 s |

El **caso de referencia** es `espacio.txt` (5000 × 5000 µm, 100 cargas, ≈ 6,5 s), que es
el que ejecuta `make run`. `ejemplo.txt` es el caso pequeño del formato de arriba,
útil para ver el dibujo con detalle.

En ambas series el tiempo escala linealmente con el número de evaluaciones
(≈ 1,2 ns por evaluación en el potencial y ≈ 1,4 ns en el campo), como corresponde a
núcleos sin dependencias entre puntos. El campo tarda alrededor de un 20 % más que el
potencial.

## Visualización

Se dibuja primero el mapa del potencial y, debajo, el del módulo del campo.
Si la malla es mayor que el área de dibujo (100 × 50 celdas) se promedian bloques
cuadrados de celdas, con el mismo factor en ambos ejes para no deformar la imagen.
Cada celda se pinta como dos espacios con color de fondo RGB (`\033[48;2;R;G;Bm`),
por lo que hace falta un terminal con soporte de color verdadero.
