// ---------------------------------------------------------------------------
// Practica 2 - Ingenieria de los Computadores (2026-27)
// Calculo secuencial del potencial y del campo electrico en una region del
// espacio.
//
// Version SECUENCIAL de referencia: no usa hebras, procesos ni ningun tipo
// de paralelismo explicito.  Hay dos nucleos de calculo independientes
// (potencial y campo), ambos con un doble bucle sobre los puntos de la malla
// que para cada punto acumula la contribucion de todas las cargas.
// ---------------------------------------------------------------------------

#include <chrono>
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
// NUCLEO DE CALCULO 1: potencial.  O(ancho * alto * n_cargas)
// Para cada punto (x, y) de la malla:  V = sum_i  k * q_i / r_i
// ---------------------------------------------------------------------------
static void calcular_potencial(int ancho, int alto,
                               const std::vector<Carga>& cargas,
                               std::vector<double>& V) {
    const int n = static_cast<int>(cargas.size());

    for (int y = 0; y < alto; ++y) {
        const double py = static_cast<double>(y);
        double* fila = &V[static_cast<size_t>(y) * ancho];

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
// NUCLEO DE CALCULO 2: modulo del campo electrico.  O(ancho * alto * n_cargas)
// Para cada punto (x, y) de la malla:  E = sum_i  k * q_i * (r - r_i) / r_i^3
// Se acumulan las componentes Ex, Ey y se guarda |E|.
// Es independiente del potencial: solo necesita las cargas.
// ---------------------------------------------------------------------------
static void calcular_campo(int ancho, int alto,
                           const std::vector<Carga>& cargas,
                           std::vector<double>& E) {
    const int n = static_cast<int>(cargas.size());

    for (int y = 0; y < alto; ++y) {
        const double py = static_cast<double>(y);
        double* fila = &E[static_cast<size_t>(y) * ancho];

        for (int x = 0; x < ancho; ++x) {
            const double px = static_cast<double>(x);
            double ex = 0.0, ey = 0.0;

            for (int i = 0; i < n; ++i) {
                const double dx = px - cargas[i].x;
                const double dy = py - cargas[i].y;
                double r = std::sqrt(dx * dx + dy * dy);
                if (r < R_MIN) r = R_MIN;
                const double f = cargas[i].q / (r * r * r);
                ex += f * dx;
                ey += f * dy;
            }

            // um en el numerador y um^3 en el denominador -> factor 1/um^2
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
    const char* ruta = (argc > 1) ? argv[1] : "espacio.txt";

    int ancho = 0, alto = 0;
    std::vector<Carga> cargas;
    if (!leer_espacio(ruta, ancho, alto, cargas)) {
        std::fprintf(stderr, "Error: no se pudo leer '%s'\n", ruta);
        return 1;
    }

    std::printf("Espacio: %d x %d um   Cargas: %zu\n", ancho, alto, cargas.size());

    std::vector<double> V(static_cast<size_t>(ancho) * alto);
    std::vector<double> E(static_cast<size_t>(ancho) * alto);

    const auto t0 = std::chrono::steady_clock::now();
    calcular_potencial(ancho, alto, cargas, V);
    const auto t1 = std::chrono::steady_clock::now();
    calcular_campo(ancho, alto, cargas, E);
    const auto t2 = std::chrono::steady_clock::now();

    std::printf("\nPotencial electrico V (rojo = alto, azul = bajo)\n");
    dibujar(ancho, alto, V, 100, 50, DIVERGENTE);
    std::printf("\nModulo del campo electrico |E| (amarillo = intenso, azul oscuro = debil)\n");
    dibujar(ancho, alto, E, 100, 50, SECUENCIAL);

    const double seg_v = std::chrono::duration<double>(t1 - t0).count();
    const double seg_e = std::chrono::duration<double>(t2 - t1).count();
    std::printf("\nTiempo de calculo del potencial: %.3f s\n", seg_v);
    std::printf("Tiempo de calculo del campo:     %.3f s\n", seg_e);
    std::printf("Tiempo de calculo total:         %.3f s\n", seg_v + seg_e);
    return 0;
}
