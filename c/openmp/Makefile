# Makefile de la version OpenMP
#   make          compila sim_openmp y prueba_arbol_omp
#   make clean    borra binarios
# Reutiliza el nucleo de ../common (se compila aqui con -fopenmp).
#
# .RECIPEPREFIX: las ordenes empiezan con '>' en vez de TAB, para que el
# archivo no se dane si el editor convierte los tabs en espacios.
.RECIPEPREFIX = >

CC      = gcc
CFLAGS  = -O2 -Wall -Wextra -std=c11 -D_POSIX_C_SOURCE=199309L -fopenmp -I../common
LDLIBS  = -lm

COMUN   = ../common/particulas.c ../common/directa.c ../common/morton.c \
          ../common/arbol.c ../common/recorrido.c

all: sim_openmp prueba_arbol_omp

sim_openmp: sim_openmp.c arbol_omp.c arbol_omp.h $(COMUN)
> $(CC) $(CFLAGS) -o $@ sim_openmp.c arbol_omp.c $(COMUN) $(LDLIBS)

prueba_arbol_omp: prueba_arbol_omp.c arbol_omp.c arbol_omp.h $(COMUN)
> $(CC) $(CFLAGS) -o $@ prueba_arbol_omp.c arbol_omp.c $(COMUN) $(LDLIBS)

clean:
> rm -f sim_openmp prueba_arbol_omp

.PHONY: all clean