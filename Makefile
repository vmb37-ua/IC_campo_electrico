# Practica 2 - Ingenieria de los Computadores
# Version secuencial (sin paralelismo).

CXX      = g++
CXXFLAGS = -O2 -std=c++11
TARGET   = campo
SRC      = campo.cpp
CASO     = espacio.txt

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

# Ejecuta el caso de prueba de referencia
run: $(TARGET)
	./$(TARGET) $(CASO)

clean:
	rm -f $(TARGET)

.PHONY: all run clean casos

# Regenera los ficheros de casos de prueba
casos:
	python3 generar_casos.py
