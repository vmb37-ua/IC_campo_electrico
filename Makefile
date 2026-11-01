# Practicas 2 y 3 - Ingenieria de los Computadores
# campo:          version secuencial (sin paralelismo).
# campo_mpi:      version MPI, estrategia A (paralelismo de datos).
# campo_mpi_func: version MPI, estrategia B (paralelismo funcional).

CXX      = g++
CXXFLAGS = -O2 -std=c++11
TARGET   = campo
SRC      = campo.cpp
CASO     = espacio.txt
MPICXX   = mpic++
NP       = 4

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

campo_mpi: campo_mpi.cpp
	$(MPICXX) $(CXXFLAGS) -o campo_mpi campo_mpi.cpp

campo_mpi_func: campo_mpi_func.cpp
	$(MPICXX) $(CXXFLAGS) -o campo_mpi_func campo_mpi_func.cpp

mpi: campo_mpi campo_mpi_func

# Ejecuta el caso de prueba de referencia
run: $(TARGET)
	./$(TARGET) $(CASO)

# Ejecucion local con MPI:  make run_mpi NP=8 CASO=casos/espacio_9000.txt
run_mpi: campo_mpi
	mpirun -np $(NP) ./campo_mpi $(CASO)

run_mpi_func: campo_mpi_func
	mpirun -np $(NP) ./campo_mpi_func $(CASO)

clean:
	rm -f $(TARGET) campo_mpi campo_mpi_func

.PHONY: all mpi run run_mpi run_mpi_func clean casos

# Regenera los ficheros de casos de prueba
casos:
	python3 generar_casos.py
