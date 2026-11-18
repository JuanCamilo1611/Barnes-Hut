/*
 * particulas.h
 * Estructura comun de las particulas y lectura/escritura de archivos.
 *
 * Los datos se guardan en arreglos planos (no un struct por particula):
 *   pos[3*i + k], vel[3*i + k], acc[3*i + k]   con k = 0 (x), 1 (y), 2 (z)
 *   masa[i]
 * Asi se pueden pasar tal cual a OpenMP, MPI y CUDA.
 */
#ifndef PARTICULAS_H
#define PARTICULAS_H

typedef struct {
    int n;          /* numero de particulas */
    double *pos;    /* 3*n posiciones  */
    double *vel;    /* 3*n velocidades */
    double *acc;    /* 3*n aceleraciones del paso actual */
    double *masa;   /* n masas */
} Sistema;

/* Reserva memoria para n particulas (todo en cero). */
Sistema crear_sistema(int n);

/* Libera la memoria del sistema. */
void liberar_sistema(Sistema *s);

/* Lee un archivo de n filas con 7 columnas: x y z vx vy vz masa.
 * Devuelve 0 si todo salio bien, -1 si hubo error. */
int leer_particulas(const char *ruta, int n, Sistema *s);

/* Escribe una linea "ax ay az" por particula, en el orden de entrada. */
int escribir_aceleraciones(const char *ruta, const Sistema *s);

#endif
