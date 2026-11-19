/*
 * recorrido.h
 * Recorrido del arbol plano SIN recursion (pila explicita) para calcular
 * aceleraciones con el criterio de apertura theta (contrato, seccion 4).
 */
#ifndef RECORRIDO_H
#define RECORRIDO_H

#include "particulas.h"
#include "arbol.h"

/* Pila maxima: por cada nivel quedan como mucho 7 hermanos pendientes. */
#define TAM_PILA (8 * (MAX_PROFUNDIDAD + 2))

/* Aceleracion sobre la particula i (escribe s->acc[3i..3i+2]). */
void aceleracion_particula(const ArbolPlano *a, Sistema *s, int i,
                           double theta, double eps);

/* Aceleraciones de todas las particulas, en un hilo. */
void aceleraciones_arbol(const ArbolPlano *a, Sistema *s,
                         double theta, double eps);

#endif
