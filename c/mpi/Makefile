# Makefile de la version MPI
#   make            compila sim_mpi (reutiliza el nucleo de c/common)
#   make clean      borra binarios y objetos
#   make prueba     verifica con 4 procesos contra la suma directa
CC      = mpicc
COMUN   = ../common
CFLAGS  = -O2 -Wall -Wextra -std=c11 -D_POSIX_C_SOURCE=199309L -I$(COMUN)
LDLIBS  = -lm

NUCLEO  = particulas.o directa.o morton.o arbol.o
PROPIOS = locales.o reparto.o let.o sim_mpi.o

all: sim_mpi

sim_mpi: $(PROPIOS) $(NUCLEO)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

# Los objetos del nucleo se compilan aqui para no tocar c/common.
%.o: $(COMUN)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $<

prueba: sim_mpi
	cd ../.. && python3 tools/verificar.py --programa "mpirun -np 4 c/mpi/sim_mpi" --hilos 4

clean:
	rm -f *.o sim_mpi

.PHONY: all clean prueba
