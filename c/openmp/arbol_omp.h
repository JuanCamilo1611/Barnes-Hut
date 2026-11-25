/*
 * arbol_omp.h
 * Construccion PARALELA del arbol plano a partir del orden Morton.
 *
 * Idea (enunciado 4.3: "a partir de un ordenamiento por Morton"):
 *   1. Esqueleto (secuencial, barato): se baja por el arbol solo mientras
 *      un nodo tenga muchas particulas. Cada rango "pequeno" se vuelve una
 *      TAREA: un subarbol completo e independiente de los demas.
 *   2. Cada hilo construye sus tareas en un arbol LOCAL (sin candados,
 *      porque los rangos Morton no se tocan entre si).
 *   3. Se calcula donde empieza cada subarbol en el arbol final (preorden)
 *      y se copian en paralelo, sumando ese desplazamiento a hijos[].
 *   4. Se llenan los pocos nodos del esqueleto y su masa / centro de masa.
 * El resultado es IDENTICO al arbol secuencial (mismo preorden).
 */
#ifndef ARBOL_OMP_H
#define ARBOL_OMP_H

#include <stdint.h>
#include "particulas.h"
#include "arbol.h"

/* Un elemento del esqueleto: un nodo de arriba o una tarea (subarbol). */
typedef struct {
    int es_tarea;
    int ini, fin, prof;      /* rango de particulas en orden Morton */
    double centro[3], lado;  /* cubo del elemento */
    int hijos[8];            /* elementos hijos (solo nodos), -1 si no hay */
    int tarea;               /* numero de tarea (solo tareas) */
    int inicio;              /* indice de su primer nodo en el arbol final */
} Elemento;

/* Memoria que se reutiliza entre pasos (para no reservar cada vez). */
typedef struct {
    Elemento *elem;  int n_elem, cap_elem;
    int *tareas;     int n_tareas, cap_tareas;   /* elemento de cada tarea */
    ArbolPlano *locales; int n_locales;          /* un arbol por tarea */
    int tareas_por_hilo;     /* granularidad: cuantas tareas por hilo */
} ConstructorOMP;

ConstructorOMP crear_constructor(int tareas_por_hilo);
void liberar_constructor(ConstructorOMP *k);

/* Igual que construir_arbol() del nucleo, pero en paralelo. */
void construir_arbol_omp(ArbolPlano *a, ConstructorOMP *k, const Sistema *s,
                         const int *orden, const uint64_t *codigos,
                         const double minimo[3], double lado);

#endif