// ---------------------------------------------------------------------------
// Practica 2 - Ingenieria de los Computadores (2026-27)
// Calculo secuencial del potencial electrico en una region del espacio.
//
// Version SECUENCIAL de referencia: no usa hebras, procesos ni ningun tipo
// de paralelismo explicito.  El nucleo de calculo es el doble bucle sobre
// los puntos de la malla, que para cada punto acumula la contribucion de
// todas las cargas.
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
// NUCLEO DE CALCULO.  O(ancho * alto * n_cargas)
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
// Representacion en el terminal con colores ANSI de 24 bits.
// rojo = potencial alto, azul = potencial bajo.
// Si la malla es mayor que el terminal se promedian bloques de celdas.
// ---------------------------------------------------------------------------
static void dibujar(int ancho, int alto, const std::vector<double>& V,
                    int cols_max, int filas_max) {
    // Factor de reduccion (el mismo en ambos ejes para no deformar la imagen).
    int paso = 1;
    while (ancho / paso > cols_max || alto / paso > filas_max) ++paso;

    const int cols  = ancho / paso;
    const int filas = alto / paso;

    // Escala: valor de referencia = media de |V|.  Con tanh() se comprimen
    // los picos enormes que aparecen junto a las cargas.
    double suma = 0.0;
    for (size_t i = 0; i < V.size(); ++i) suma += std::fabs(V[i]);
    const double vref = (suma / V.size()) + 1e-12;

    std::string salida;
    salida.reserve(static_cast<size_t>(cols) * filas * 24);
    char buf[48];

    for (int f = 0; f < filas; ++f) {
        for (int c = 0; c < cols; ++c) {
            // Media del bloque paso x paso.
            double m = 0.0;
            for (int dy = 0; dy < paso; ++dy)
                for (int dx = 0; dx < paso; ++dx)
                    m += V[static_cast<size_t>(f * paso + dy) * ancho + c * paso + dx];
            m /= static_cast<double>(paso) * paso;

            const double t = std::tanh(m / vref);          // -1 .. 1
            const int rojo  = static_cast<int>(255.0 * (t > 0 ?  t : 0.0));
            const int azul  = static_cast<int>(255.0 * (t < 0 ? -t : 0.0));
            const int verde = static_cast<int>(60.0 * (1.0 - std::fabs(t)));

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

    const auto t0 = std::chrono::steady_clock::now();
    calcular_potencial(ancho, alto, cargas, V);
    const auto t1 = std::chrono::steady_clock::now();

    dibujar(ancho, alto, V, 100, 50);

    const double seg = std::chrono::duration<double>(t1 - t0).count();
    std::printf("Tiempo de calculo: %.3f s\n", seg);
    return 0;
}
