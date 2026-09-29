# IC_campo_electrico

Cálculo secuencial del **potencial eléctrico** en una región plana del espacio.
Dadas *n* cargas puntuales colocadas en esa región, el programa evalúa el potencial
en cada punto de una malla y lo representa en el terminal con colores ANSI de 24 bits
(**rojo = potencial alto, azul = potencial bajo, verde oscuro ≈ potencial nulo**).
Al terminar imprime el tiempo de cálculo en segundos.

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
- Se trabaja en 2D: la malla es un corte plano de la región, y las cargas se
  consideran contenidas en ese mismo plano.

### Unidades

| Magnitud | Unidad en el fichero | Conversión interna |
|---|---|---|
| Tamaño del espacio | micrómetros (µm) | 1 celda de la malla = 1 µm |
| Posición de las cargas `(x, y)` | micrómetros (µm) | × 10⁻⁶ → metros |
| Valor de la carga `q` | nanoculombios (nC) | × 10⁻⁹ → culombios |
| Potencial calculado `V` | **voltios (V)** | — |

Todo el cálculo se hace en `double`.

### Singularidad en las cargas

`V → ∞` cuando `r → 0`. Para evitar la división por cero en las celdas que caen justo
sobre una carga se aplica un radio mínimo de `R_MIN = 0.5 µm` (la mitad del paso de la
malla): si `r < R_MIN`, se usa `R_MIN`. Es el equivalente al *softening* habitual en
simulaciones de N cuerpos.

### Escala de color

Los potenciales crecen sin límite junto a las cargas, así que una escala lineal
saturaría la imagen. Se normaliza con

```
t = tanh( V / V_ref )        V_ref = media de |V| sobre toda la malla
```

de modo que `t ∈ (−1, 1)`, y se mapea `t = +1 → rojo`, `t = −1 → azul`, `t = 0 → verde oscuro`.

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

El núcleo (`calcular_potencial()`) es un triple bucle de coste

```
O(ancho · alto · n_cargas)
```

con una raíz cuadrada y una división por evaluación. La memoria ocupada por la malla es
`ancho · alto · 8 bytes`.

### Casos de prueba incluidos

Se regeneran con `python3 generar_casos.py` (posiciones deterministas, semilla fija).
Tiempos medidos con `g++ -O2` en un único núcleo:

**Serie A — escalado del tamaño del espacio (100 cargas fijas)**

| Fichero | Espacio (µm) | Evaluaciones | Tiempo | Pico de RSS |
|---|---|---|---|---|
| `casos/espacio_2000.txt`  | 2000 × 2000   | 4,0 · 10⁸ | 0,95 s | — |
| `casos/espacio_3000.txt`  | 3000 × 3000   | 9,0 · 10⁸ | 2,16 s | — |
| `casos/espacio_5000.txt`  | 5000 × 5000   | 2,5 · 10⁹ | 5,98 s | 200 MB |
| `casos/espacio_7000.txt`  | 7000 × 7000   | 4,9 · 10⁹ | 11,73 s | 386 MB |
| `casos/espacio_9000.txt`  | 9000 × 9000   | 8,1 · 10⁹ | 19,60 s | 636 MB |
| `casos/espacio_11000.txt` | 11000 × 11000 | 1,21 · 10¹⁰ | 31,01 s | 949 MB |

**Serie B — escalado del número de cargas (espacio fijo de 4000 × 4000 µm)**

| Fichero | Cargas | Evaluaciones | Tiempo |
|---|---|---|---|
| `casos/cargas_50.txt`  | 50  | 8,0 · 10⁸ | 1,91 s |
| `casos/cargas_100.txt` | 100 | 1,6 · 10⁹ | 3,83 s |
| `casos/cargas_200.txt` | 200 | 3,2 · 10⁹ | 7,67 s |
| `casos/cargas_400.txt` | 400 | 6,4 · 10⁹ | 15,29 s |

El **caso de referencia** es `espacio.txt` (5000 × 5000 µm, 100 cargas, ≈ 6 s), que es
el que ejecuta `make run`. `ejemplo.txt` es el caso pequeño del formato de arriba,
útil para ver el dibujo con detalle.

En ambas series el tiempo escala linealmente con el número de evaluaciones
(≈ 2,4 ns por evaluación), como corresponde a un núcleo sin dependencias entre puntos.

## Visualización

Si la malla es mayor que el área de dibujo (100 × 50 celdas) se promedian bloques
cuadrados de celdas, con el mismo factor en ambos ejes para no deformar la imagen.
Cada celda se pinta como dos espacios con color de fondo RGB (`\033[48;2;R;G;Bm`),
por lo que hace falta un terminal con soporte de color verdadero.
