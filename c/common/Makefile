# Makefile del nucleo C
#   make          compila sim_directa y prueba_morton
#   make clean    borra binarios y objetos
CC      = gcc
CFLAGS  = -O2 -Wall -Wextra -std=c11 -D_POSIX_C_SOURCE=199309L
LDLIBS  = -lm

NUCLEO  = particulas.o directa.o leapfrog.o morton.o

all: sim_directa prueba_morton

sim_directa: sim_directa.o $(NUCLEO)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

prueba_morton: prueba_morton.o particulas.o morton.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o sim_directa prueba_morton

.PHONY: all clean
