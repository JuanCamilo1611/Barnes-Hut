# Makefile de la version CUDA (para Colab, GPU T4 = sm_75)
#   make         compila sim_cuda
#   make clean
.RECIPEPREFIX = >
NVCC    = nvcc
FLAGS   = -O2 -arch=sm_75 -I../c/common -Xcompiler -D_POSIX_C_SOURCE=199309L
COMUN   = ../c/common/particulas.c ../c/common/directa.c \
          ../c/common/morton.c ../c/common/arbol.c

sim_cuda: sim_cuda.cu $(COMUN)
> $(NVCC) $(FLAGS) -o $@ sim_cuda.cu $(COMUN) -lm

clean:
> rm -f sim_cuda

.PHONY: clean