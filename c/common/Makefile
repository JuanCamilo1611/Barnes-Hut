# Makefile del nucleo C
#   make          compila sim_directa
#   make clean    borra binarios y objetos
CC      = gcc
CFLAGS  = -O2 -Wall -Wextra -std=c11 -D_POSIX_C_SOURCE=199309L
LDLIBS  = -lm

NUCLEO  = particulas.o directa.o leapfrog.o

all: sim_directa

sim_directa: sim_directa.o $(NUCLEO)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o sim_directa

.PHONY: all clean
