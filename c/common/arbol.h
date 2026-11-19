/*
 * arbol.h
 * Octree plano (SoA) construido desde el orden Morton.
 * Sigue docs/contrato_arbol_plano.md:
 *   - nodo 0 = raiz;
 *   - numeracion en PREORDEN: cada nodo va antes que sus descendientes y
 *     el subarbol de un nodo es un rango contiguo de indices;
 *   - hijos[8*n + k] = -1 si el octante k esta vacio;
 *   - part[n] >= 0 solo en hojas; -1 en nodos internos;
 *   - max_profundidad = 20 con bucketing en el tope.
 */
#ifndef ARBOL_H
#define ARBOL_H

#include <stdint.h>
#include "particulas.h"

#define MAX_PROFUNDIDAD 20   /* contrato, seccion 5 */

typedef struct {
    int cap;          /* nodos reservados */
    int n_nodos;      /* nodos usados */
    /* --- los 6 arreglos del contrato --- */
    double  *centro;  /* 3*cap: centro geometrico del cubo */
    double  *lado;    /* cap:   arista del cubo */
    double  *masa;    /* cap:   masa total de la celda */
    double  *com;     /* 3*cap: centro de masa */
    int32_t *hijos;   /* 8*cap: hijo por octante, -1 si no existe */
    int32_t *part;    /* cap:   particula de la hoja, -1 si interno */
    /* --- auxiliares SOLO de construccion (no se comparten) --- */
    double  *hojas_masa;      /* cap */
    double  *hojas_momento;   /* 3*cap: masa * posicion */
} ArbolPlano;

/* Reserva un arbol vacio con capacidad inicial (crece solo si hace falta). */
ArbolPlano crear_arbol(int cap_inicial);
void liberar_arbol(ArbolPlano *a);

/* Construye el arbol a partir de las particulas ya ordenadas por Morton.
 *   orden[]   indices de particula en orden Morton (ordenar_por_morton)
 *   codigos[] codigo Morton de cada particula (indexado por particula)
 *   minimo, lado: el cubo raiz (caja_envolvente)
 * Deja masa[] y com[] ya calculados. */
void construir_arbol(ArbolPlano *a, const Sistema *s,
                     const int *orden, const uint64_t *codigos,
                     const double minimo[3], double lado);

/* Profundidad maxima alcanzada (para el informe y las pruebas). */
int profundidad_arbol(const ArbolPlano *a);

#endif
