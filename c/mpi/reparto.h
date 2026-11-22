/*
 * reparto.h
 * Particionado del espacio entre procesos MPI.
 *
 * Cada particula recibe una CLAVE de 64 bits (clave_particula). Las
 * claves se ordenan y se cortan en P tramos de igual PESO, y cada
 * proceso se queda con un tramo. El peso es el costo medido en el paso
 * anterior (interacciones del recorrido) o 1 por particula.
 *
 *   PARTICION_MORTON     clave = codigo Morton -> cada proceso recibe un
 *                        tramo de la curva Z: una region compacta.
 *   PARTICION_BISECCION  clave = coordenada x -> cada proceso recibe una
 *                        rebanada del cubo a lo largo de x.
 *
 * Para cambiar la estrategia basta con la variable de entorno
 * BH_PARTICION (morton | biseccion) o con agregar un caso nuevo en
 * clave_particula() de reparto.c.
 */
#ifndef REPARTO_H
#define REPARTO_H

#include <mpi.h>
#include "locales.h"

typedef enum { PARTICION_MORTON, PARTICION_BISECCION } Particion;

/* Cubo raiz COMUN a todos los procesos (MPI_Allreduce del min y max),
 * con el mismo margen que caja_envolvente() del nucleo C. */
void caja_global(const Locales *L, double minimo[3], double *lado,
                 MPI_Comm comm);

/* Redistribuye las particulas: al volver, el proceso r tiene las del
 * tramo r de claves.  usar_costo = 0 reparte por numero de particulas. */
void repartir(Locales *L, const double minimo[3], double lado,
              Particion modo, int usar_costo, MPI_Comm comm);

#endif
