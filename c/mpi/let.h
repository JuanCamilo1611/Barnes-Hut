/*
 * let.h
 * Locally Essential Tree: la parte del arbol de un proceso que OTRO
 * proceso necesita para calcular sus fuerzas.
 *
 * Para cada proceso destino se recorre el arbol local en PREORDEN contra
 * la caja (AABB) de las particulas del destino:
 *   - si la celda esta lejos de TODA la caja (lado < theta * distancia
 *     minima a la caja), se envia como un solo nodo y se SALTA su rama
 *     entera: i = fin[i] + 1 (posible porque en preorden la rama es un
 *     rango contiguo, contrato seccion 3);
 *   - si no, se envia el nodo y se sigue con sus hijos.
 * Lo que se recibe de cada proceso es un ArbolPlano del contrato, mas
 * pequeno: sus hojas (part = PART_REMOTA) son particulas remotas o celdas
 * resumidas. Ningun proceso tiene el arbol global completo.
 */
#ifndef LET_H
#define LET_H

#include <mpi.h>
#include "arbol.h"
#include "locales.h"

/* Marca de hoja remota: nunca coincide con un indice local, asi que el
 * recorrido no la confunde con la propia particula. */
#define PART_REMOTA INT32_MAX

typedef struct {
    int P;
    ArbolPlano *arbol;      /* P arboles recibidos (el propio queda vacio) */
    long nodos_enviados;    /* estadisticas del ultimo intercambio */
    long nodos_recibidos;
} LETRecibidos;

/* Arma y envia los LET de este proceso a todos los demas (MPI_Alltoallv)
 * y recibe los que le corresponden. */
void intercambiar_let(const ArbolPlano *local, const Locales *L,
                      double theta, LETRecibidos *rec, MPI_Comm comm);

void liberar_let(LETRecibidos *rec);

/* Suma en acc[3] la aceleracion que el arbol t ejerce sobre el punto p,
 * excluyendo la hoja cuya part == excluir. Devuelve las interacciones
 * evaluadas (costo).  Mismo criterio que c/common/recorrido.c. */
long acumular_arbol(const ArbolPlano *t, const double p[3], int excluir,
                    double theta, double eps, double acc[3]);

#endif
