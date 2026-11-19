/*
 * directa.h
 * Suma directa O(N^2) y energia total. Es la referencia de correctitud
 * en C: Barnes-Hut debe dar casi las mismas aceleraciones.
 *
 * Constante G = 1. El suavizado eps evita fuerzas infinitas:
 *   a_i = sum_j  G m_j (r_j - r_i) / (|r_j - r_i|^2 + eps^2)^(3/2)
 */
#ifndef DIRECTA_H
#define DIRECTA_H

#include "particulas.h"

#define G_GRAVEDAD 1.0

/* Calcula s->acc para todas las particulas sumando todos los pares. */
void aceleraciones_directas(Sistema *s, double eps);

/* Energia cinetica + potencial (cada pareja se cuenta una vez). O(N^2). */
double energia_total(const Sistema *s, double eps);

#endif
