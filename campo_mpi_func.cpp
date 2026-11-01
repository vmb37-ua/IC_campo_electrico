// ---------------------------------------------------------------------------
// Practica 3 - Ingenieria de los Computadores (2026-27)
// Potencial y campo electrico con MPI - estrategia B (funcional).
//
// ESTRATEGIA B: paralelismo FUNCIONAL.  Los procesos se dividen en dos grupos
// que ejecutan funciones distintas a la vez: uno calcula el potencial V y el
// otro el campo |E| (son independientes, solo dependen de las cargas).
// Dentro de cada grupo las filas se reparten como en la estrategia A.
//
//   rank 0: lee el fichero;  MPI_Bcast (dims, cargas) -> todos
//   MPI_Comm_split:  grupo V = ranks 0 .. nV-1   |  grupo E = ranks nV .. N-1
//   grupo V: calcula sus filas de V              |  grupo E: sus filas de E
//   MPI_Gatherv en el grupo -> rank 0 (V)        |  MPI_Gatherv -> rank nV (E)
//                                                |  MPI_Send (E) -> rank 0
//   rank 0: MPI_Recv (E), dibuja e imprime los tiempos
//
// Compilar:  make campo_mpi_func
// Ejecutar:  mpirun -np 4 ./campo_mpi_func casos/espacio_5000.txt
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

    // Con un solo proceso no se pueden formar dos grupos.
    if (nprocs < 2) {
        if (rank == 0) std::fprintf(stderr, "La estrategia B necesita al menos 2 procesos\n");
        MPI_Finalize();
        return 1;
    }

    // --- Division en dos grupos ----------------------------------------------
    // nV procesos para el potencial; el resto para el campo. Con un numero
    // impar de procesos el campo se lleva uno mas (es ~20 % mas caro).
    const int nV    = nprocs / 2;
    const int color = (rank < nV) ? 0 : 1;   // 0 = potencial, 1 = campo
    MPI_Comm grupo;
    // Divide MPI_COMM_WORLD en comunicadores nuevos: los procesos con el
    // mismo color quedan en el mismo. El tercer argumento (rank) fija el
    // orden dentro del grupo. Es colectiva: la llaman todos.
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &grupo);
    int grank, gsize;
    // Rank y tamano DENTRO del grupo: el world rank nV es el grank 0 del
    // grupo del campo.
    MPI_Comm_rank(grupo, &grank);
    MPI_Comm_size(grupo, &gsize);

    // --- Reparto de filas dentro del grupo ------------------------------------
    std::vector<int> counts, displs;
    repartir_filas(alto, ancho, gsize, counts, displs);
    const int y0    = displs[grank] / ancho;
    const int filas = counts[grank] / ancho;

    // --- Memoria -------------------------------------------------------------
    // loc: las filas de este proceso (de V o de E segun el grupo).
    // malla: la malla completa del grupo, solo en su raiz (grank 0), es decir,
    // world rank 0 para V y world rank nV para E.
    std::vector<double> loc(static_cast<size_t>(filas) * ancho);
    std::vector<double> malla;
    if (grank == 0) malla.resize(static_cast<size_t>(ancho) * alto);

    // --- Calculo -------------------------------------------------------------
    // Barrera sobre TODOS los procesos para que los dos grupos empiecen a
    // medir a la vez.
    MPI_Barrier(MPI_COMM_WORLD);
    const double t0 = MPI_Wtime();
    // Cada grupo ejecuta una FUNCION distinta: paralelismo funcional.
    if (color == 0) calcular_potencial(ancho, y0, filas, cargas, loc);
    else            calcular_campo(ancho, y0, filas, cargas, loc);
    const double t1 = MPI_Wtime();

    // --- Reunir resultados, fase 1: dentro de cada grupo ---------------------
    // Igual que en A, pero el comunicador es 'grupo': cada grupo junta sus
    // filas en SU raiz y los dos Gatherv ocurren a la vez.
    MPI_Gatherv(loc.data(), counts[grank], MPI_DOUBLE,
                malla.data(), counts.data(), displs.data(), MPI_DOUBLE,
                0, grupo);

    // --- Reunir resultados, fase 2: E del rank nV al rank 0 ------------------
    // Comunicacion PUNTO A PUNTO (un emisor y un receptor).
    std::vector<double>& V = malla;   // en el rank 0, 'malla' es el potencial
    std::vector<double> E;
    const int total = ancho * alto;   // cabe en int incluso con 11000 x 11000
    if (rank == nV) {
        // Envia la malla de E completa al world rank 0 con etiqueta 0.
        // Bloquea hasta que el buffer se puede reutilizar.
        MPI_Send(malla.data(), total, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    } else if (rank == 0) {
        E.resize(total);
        // Bloquea hasta recibir del rank nV un mensaje con etiqueta 0.
        // MPI_STATUS_IGNORE: no se necesita la informacion del mensaje.
        MPI_Recv(E.data(), total, MPI_DOUBLE, nV, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    const double t2 = MPI_Wtime();

    // --- Tiempo global -------------------------------------------------------
    // Se separa el calculo de cada grupo para ver el desequilibrio entre V y
    // E: cada proceso pone su tiempo en la posicion de su grupo y 0 en la
    // otra. La reduccion con MPI_MAX deja en el rank 0 el maximo de cada
    // posicion (calculo de V, calculo de E y total).
    double tloc[3] = {color == 0 ? t1 - t0 : 0.0,
                      color == 1 ? t1 - t0 : 0.0,
                      t2 - t0};
    double tmax[3];
    MPI_Reduce(tloc, tmax, 3, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // --- Salida (solo el rank 0) ----------------------------------------------
    if (rank == 0) {
        std::printf("Espacio: %d x %d um   Cargas: %d   Procesos: %d (V: %d, E: %d)\n",
                    ancho, alto, n, nprocs, nV, nprocs - nV);
        std::printf("\nPotencial electrico V (rojo = alto, azul = bajo)\n");
        dibujar(ancho, alto, V, 100, 50, DIVERGENTE);
        std::printf("\nModulo del campo electrico |E| (amarillo = intenso, azul oscuro = debil)\n");
        dibujar(ancho, alto, E, 100, 50, SECUENCIAL);

        const double calc = (tmax[0] > tmax[1]) ? tmax[0] : tmax[1];
        std::printf("\nTiempo de calculo grupo V (max): %.3f s\n", tmax[0]);
        std::printf("Tiempo de calculo grupo E (max): %.3f s\n", tmax[1]);
        std::printf("Tiempo de comunicacion:          %.3f s\n", tmax[2] - calc);
        std::printf("Tiempo total:                    %.3f s\n", tmax[2]);
    }

    // Libera el comunicador creado con MPI_Comm_split.
    MPI_Comm_free(&grupo);
    // Cierra MPI. Ultima llamada MPI; la hacen todos los procesos.
    MPI_Finalize();
    return 0;
}
