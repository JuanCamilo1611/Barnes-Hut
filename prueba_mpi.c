
#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[]) {

    // Inicializar MPI
    MPI_Init(&argc, &argv);

    // Obtener identificador del proceso
    int id;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    // Obtener cantidad total de procesos
    int total;
    MPI_Comm_size(MPI_COMM_WORLD, &total);

    printf("Hola desde el proceso %d de %d\n",
           id, total);

    // Finalizar MPI
    MPI_Finalize();

    return 0;
}


