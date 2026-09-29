# IC_campo_electrico

Programa que calcula el potencial eléctrico en una región del espacio (m x m, en micrómetros) dado en C++. Con n cargas predefinidas en el espacio se calcula el valor del potencial eléctrico en cada punto. Se representa en la terminal dicho espacio con colores de códigos de escape ANSI con: rojo=alto potencial; azul=bajo potencial.

Después, imprimirá el tiempo que ha tardado en ejecutarse el programa en segundos.

Modo de uso:
En el fichero "espacio.txt" se define el espacio y las cargas de la siquiente forma:
100x100
(5, 34, 5)
(-4, 20, 60)

El valor de carga es el primer argumento de los paréntesis que representan las cargas, la posición después.
La coordenada 0,0 corresponde a la esquina superior izquierda.
