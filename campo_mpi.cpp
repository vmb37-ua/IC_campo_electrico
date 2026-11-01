// ---------------------------------------------------------------------------
// Practica 3 - Ingenieria de los Computadores (2026-27)
// Potencial y campo electrico con MPI - estrategia A (datos).
//
// ESTRATEGIA A: paralelismo de DATOS.  Todos los procesos ejecutan las mismas
// funciones (potencial y campo) sobre un bloque de filas distinto de la malla.
//
//   rank 0: lee el fichero
//      |  MPI_Bcast (ancho, alto, n) y MPI_Bcast (cargas)  -> todos
//      v
//   cada proceso calcula V y E de SUS filas               <- parte paralela
//      |  MPI_Gatherv (V) y MPI_Gatherv (E)                -> rank 0
//      v
//   rank 0: dibuja e imprime los tiempos (MPI_Reduce con el maximo)
//
// Compilar:  make campo_mpi
// Ejecutar:  mpirun -np 4 ./campo_mpi casos/espacio_5000.txt
// ---------------------------------------------------------------------------

#include <mpi.h>        // biblioteca MPI
#include <algorithm>    // std::min
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

// --- Constantes fisicas -----------------------------------------------------
// Las posiciones del fichero estan en micrometros y las cargas en nanoculombios.
static const double K_COULOMB = 8.9875517923e9;  // N*m^2/C^2
static const double NANO      = 1e-9;            // nC -> C
static const double MICRO     = 1e-6;            // um -> m
static const double R_MIN     = 0.5;             // radio minimo (um), evita 1/0

struct Carga {
    double q;  // nC
    double x;  // um
    double y;  // um
};

// ---------------------------------------------------------------------------
// Lectura del fichero de configuracion (todo lo no computacional, lo minimo)
// Formato:   AnchoxAlto
//            (carga, x, y)
// ---------------------------------------------------------------------------
static bool leer_espacio(const char* ruta, int& ancho, int& alto,
                         std::vector<Carga>& cargas) {
    std::ifstream f(ruta);
    if (!f) return false;

    std::string linea;
    bool tam_leido = false;

    while (std::getline(f, linea)) {
        Carga c;
        if (!tam_leido) {
            if (std::sscanf(linea.c_str(), " %dx%d", &ancho, &alto) == 2) {
                tam_leido = true;
            }
        } else if (std::sscanf(linea.c_str(), " ( %lf , %lf , %lf )",
                               &c.q, &c.x, &c.y) == 3) {
            cargas.push_back(c);
        }
    }
    return tam_leido && !cargas.empty();
}

// ---------------------------------------------------------------------------
// Reparto de la malla por BLOQUES DE FILAS CONTIGUOS entre nprocs procesos.
// Rellena los dos arrays que necesita MPI_Gatherv:
//   counts[p] = cuantos doubles calcula el proceso p (filas * ancho)
//   displs[p] = donde empieza su trozo dentro de la malla completa
// Si alto no es divisible entre nprocs, los 'resto' primeros procesos se
// llevan una fila mas: la diferencia maxima es de 1 fila (carga equilibrada).
// ---------------------------------------------------------------------------
static void repartir_filas(int alto, int ancho, int nprocs,
                           std::vector<int>& counts, std::vector<int>& displs) {
    counts.resize(nprocs);
    displs.resize(nprocs);
    const int base  = alto / nprocs;
    const int resto = alto % nprocs;
    for (int p = 0; p < nprocs; ++p) {
        const int y0    = p * base + std::min(p, resto);
        const int filas = base + (p < resto ? 1 : 0);
        counts[p] = filas * ancho;
        displs[p] = y0 * ancho;
    }
}

// ---------------------------------------------------------------------------
// NUCLEO DE CALCULO 1: potencial de las filas [y0, y0 + filas).
// V es el buffer LOCAL del proceso: solo contiene sus filas.
// ---------------------------------------------------------------------------
static void calcular_potencial(int ancho, int y0, int filas,
                               const std::vector<Carga>& cargas,
                               std::vector<double>& V) {
    const int n = static_cast<int>(cargas.size());

    for (int f = 0; f < filas; ++f) {
        const double py = static_cast<double>(y0 + f);       // fila GLOBAL
        double* fila = &V[static_cast<size_t>(f) * ancho];   // fila LOCAL

        for (int x = 0; x < ancho; ++x) {
            const double px = static_cast<double>(x);
            double v = 0.0;

            for (int i = 0; i < n; ++i) {
                const double dx = px - cargas[i].x;
                const double dy = py - cargas[i].y;
                double r = std::sqrt(dx * dx + dy * dy);
                if (r < R_MIN) r = R_MIN;
                v += cargas[i].q / r;
            }

            fila[x] = v * (K_COULOMB * NANO / MICRO);  // voltios
        }
    }
}

// ---------------------------------------------------------------------------
// NUCLEO DE CALCULO 2: modulo del campo de las filas [y0, y0 + filas).
// ---------------------------------------------------------------------------
static void calcular_campo(int ancho, int y0, int filas,
                           const std::vector<Carga>& cargas,
                           std::vector<double>& E) {
    const int n = static_cast<int>(cargas.size());

    for (int f = 0; f < filas; ++f) {
        const double py = static_cast<double>(y0 + f);       // fila GLOBAL
        double* fila = &E[static_cast<size_t>(f) * ancho];   // fila LOCAL

        for (int x = 0; x < ancho; ++x) {
            const double px = static_cast<double>(x);
            double ex = 0.0, ey = 0.0;

            for (int i = 0; i < n; ++i) {
                const double dx = px - cargas[i].x;
                const double dy = py - cargas[i].y;
                double r = std::sqrt(dx * dx + dy * dy);
                if (r < R_MIN) r = R_MIN;
                const double f3 = cargas[i].q / (r * r * r);
                ex += f3 * dx;
                ey += f3 * dy;
            }

            fila[x] = std::sqrt(ex * ex + ey * ey)
                    * (K_COULOMB * NANO / (MICRO * MICRO));  // V/m
        }
    }
}

// ---------------------------------------------------------------------------
// Representacion en el terminal con colores ANSI de 24 bits.
//   Escala DIVERGENTE (potencial): rojo = alto, azul = bajo, verde ~ 0.
//   Escala SECUENCIAL (|E| >= 0): azul oscuro = debil, amarillo = intenso.
// Si la malla es mayor que el terminal se promedian bloques de celdas.
// ---------------------------------------------------------------------------
enum Escala { DIVERGENTE, SECUENCIAL };

static void dibujar(int ancho, int alto, const std::vector<double>& D,
                    int cols_max, int filas_max, Escala escala) {
    // Factor de reduccion (el mismo en ambos ejes para no deformar la imagen).
    int paso = 1;
    while (ancho / paso > cols_max || alto / paso > filas_max) ++paso;

    const int cols  = ancho / paso;
    const int filas = alto / paso;

    // Escala: valor de referencia = media de |D|.  Con tanh() se comprimen
    // los picos enormes que aparecen junto a las cargas.
    double suma = 0.0;
    for (size_t i = 0; i < D.size(); ++i) suma += std::fabs(D[i]);
    const double ref = (suma / D.size()) + 1e-12;

    std::string salida;
    salida.reserve(static_cast<size_t>(cols) * filas * 24);
    char buf[48];

    for (int f = 0; f < filas; ++f) {
        for (int c = 0; c < cols; ++c) {
            // Media del bloque paso x paso.
            double m = 0.0;
            for (int dy = 0; dy < paso; ++dy)
                for (int dx = 0; dx < paso; ++dx)
                    m += D[static_cast<size_t>(f * paso + dy) * ancho + c * paso + dx];
            m /= static_cast<double>(paso) * paso;

            const double t = std::tanh(m / ref);  // -1 .. 1  (0 .. 1 si D >= 0)
            int rojo, verde, azul;
            if (escala == DIVERGENTE) {
                rojo  = static_cast<int>(255.0 * (t > 0 ?  t : 0.0));
                azul  = static_cast<int>(255.0 * (t < 0 ? -t : 0.0));
                verde = static_cast<int>(60.0 * (1.0 - std::fabs(t)));
            } else {
                rojo  = static_cast<int>(255.0 * t);
                verde = static_cast<int>(220.0 * t * t);
                azul  = static_cast<int>(80.0 * (1.0 - t));
            }

            std::snprintf(buf, sizeof(buf), "\033[48;2;%d;%d;%dm  ", rojo, verde, azul);
            salida += buf;
        }
        salida += "\033[0m\n";
    }
    std::fwrite(salida.data(), 1, salida.size(), stdout);
}

// ---------------------------------------------------------------------------
int main(int argc, char** argv) {
    // Arranca el entorno MPI. Debe ser la primera llamada MPI del programa;
    // recibe argc/argv porque mpirun puede anadir argumentos propios.
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    // rank = identificador de este proceso (0 .. nprocs-1) dentro de
    // MPI_COMM_WORLD, el comunicador que agrupa a TODOS los procesos.
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // nprocs = numero total de procesos lanzados (el -np de mpirun).
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    const char* ruta = (argc > 1) ? argv[1] : "espacio.txt";

    // --- Lectura en el rank 0 y difusion de los datos ----------------------
    // ancho, alto y n se empaquetan en un array para enviarlos en UN solo
    // mensaje (cada mensaje paga una latencia fija).
    int dims[3] = {0, 0, 0};
    std::vector<Carga> cargas;
    if (rank == 0) {
        int a = 0, h = 0;
        if (leer_espacio(ruta, a, h, cargas)) {
            dims[0] = a;
            dims[1] = h;
            dims[2] = static_cast<int>(cargas.size());
        } else {
            std::fprintf(stderr, "Error: no se pudo leer '%s'\n", ruta);
        }
    }
    // Difusion (colectiva): el rank 0 (raiz) envia los 3 enteros de dims y
    // el resto los recibe en su propio dims. La llaman TODOS los procesos.
    MPI_Bcast(dims, 3, MPI_INT, 0, MPI_COMM_WORLD);
    // Si la lectura fallo, dims llega a cero a todos y terminan todos a la
    // vez (si solo terminase el rank 0, el resto se quedaria bloqueado).
    if (dims[2] == 0) { MPI_Finalize(); return 1; }
    const int ancho = dims[0], alto = dims[1], n = dims[2];

    // MPI no reserva memoria: el buffer de recepcion debe tener ya su tamano.
    if (rank != 0) cargas.resize(n);
    // Difusion de las cargas. Carga son 3 double contiguos y std::vector
    // guarda los elementos seguidos, asi que el vector es un array de 3*n
    // doubles que se envia en un solo mensaje.
    MPI_Bcast(cargas.data(), 3 * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // --- Reparto de filas ---------------------------------------------------
    // Cada proceso calcula el reparto de TODOS (es barato y asi no hay que
    // comunicarlo) y se queda con su primera fila y su numero de filas.
    std::vector<int> counts, displs;
    repartir_filas(alto, ancho, nprocs, counts, displs);
    const int y0    = displs[rank] / ancho;
    const int filas = counts[rank] / ancho;

    // --- Memoria -------------------------------------------------------------
    // Buffers locales: cada proceso guarda solo sus filas (1/nprocs de la
    // malla). Las mallas completas solo las necesita el rank 0, que dibuja.
    std::vector<double> Vloc(static_cast<size_t>(filas) * ancho);
    std::vector<double> Eloc(static_cast<size_t>(filas) * ancho);
    std::vector<double> V, E;
    if (rank == 0) {
        V.resize(static_cast<size_t>(ancho) * alto);
        E.resize(static_cast<size_t>(ancho) * alto);
    }

    // --- Calculo -------------------------------------------------------------
    // Barrera: nadie sigue hasta que todos han llegado. Asi todos empiezan a
    // medir a la vez (el rank 0, que leyo el fichero, iria por detras).
    MPI_Barrier(MPI_COMM_WORLD);
    // Reloj de MPI (segundos, double); sustituye a std::chrono.
    const double t0 = MPI_Wtime();
    calcular_potencial(ancho, y0, filas, cargas, Vloc);
    calcular_campo(ancho, y0, filas, cargas, Eloc);
    const double t1 = MPI_Wtime();

    // --- Reunir los resultados en el rank 0 ----------------------------------
    // Recogida (colectiva): cada proceso envia sus counts[rank] doubles y el
    // rank 0 coloca el trozo del proceso p en V a partir de displs[p].
    // Version "v" porque los trozos pueden tener tamanos distintos.
    // En los procesos que no son raiz, el buffer de recepcion se ignora.
    MPI_Gatherv(Vloc.data(), counts[rank], MPI_DOUBLE,
                V.data(), counts.data(), displs.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    // Lo mismo para el campo.
    MPI_Gatherv(Eloc.data(), counts[rank], MPI_DOUBLE,
                E.data(), counts.data(), displs.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    // t2 - t1 = comunicacion + espera al proceso mas lento (el Gatherv
    // sincroniza implicitamente en el raiz).
    const double t2 = MPI_Wtime();

    // --- Tiempo global -------------------------------------------------------
    // Reduccion (colectiva): combina los tloc de todos con MPI_MAX, elemento
    // a elemento, y deja el resultado en tmax del rank 0. Se toma el maximo
    // porque el programa no acaba hasta que acaba el proceso mas lento.
    double tloc[2] = {t1 - t0, t2 - t0};
    double tmax[2];
    MPI_Reduce(tloc, tmax, 2, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // --- Salida (solo el rank 0, el unico con las mallas completas) ----------
    if (rank == 0) {
        std::printf("Espacio: %d x %d um   Cargas: %d   Procesos: %d\n",
                    ancho, alto, n, nprocs);
        std::printf("\nPotencial electrico V (rojo = alto, azul = bajo)\n");
        dibujar(ancho, alto, V, 100, 50, DIVERGENTE);
        std::printf("\nModulo del campo electrico |E| (amarillo = intenso, azul oscuro = debil)\n");
        dibujar(ancho, alto, E, 100, 50, SECUENCIAL);

        std::printf("\nTiempo de calculo (max):         %.3f s\n", tmax[0]);
        std::printf("Tiempo de comunicacion:          %.3f s\n", tmax[1] - tmax[0]);
        std::printf("Tiempo total:                    %.3f s\n", tmax[1]);
    }

    // Cierra MPI. Ultima llamada MPI; la hacen todos los procesos.
    MPI_Finalize();
    return 0;
}
