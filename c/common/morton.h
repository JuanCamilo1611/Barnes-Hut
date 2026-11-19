/*
 * morton.h
 * Codigos Morton (curva Z) de 63 bits y ordenamiento de particulas.
 *
 * Cada eje se convierte en un entero de 21 bits (0 .. 2^21 - 1) y los
 * bits se intercalan como x y z x y z ... con el bit de x como el MAS
 * significativo. Asi los 3 bits mas altos del codigo son exactamente el
 * octante del contrato (docs/contrato_arbol_plano.md, seccion 2):
 *     octante = 4*(x >= cx) + 2*(y >= cy) + 1*(z >= cz)
 * y en general los bits 3d..3d+2 (contando desde arriba) dan el octante
 * en la profundidad d. 21 bits por eje alcanzan para max_profundidad = 20.
 */
#ifndef MORTON_H
#define MORTON_H

#include <stdint.h>
#include "particulas.h"

#define MORTON_BITS 21   /* bits por eje: 3 * 21 = 63 bits por codigo */

/* Cubo que contiene a todas las particulas: esquina minima y arista.
 * Se le suma un margen pequeno para que ninguna quede en el borde. */
void caja_envolvente(const Sistema *s, double minimo[3], double *lado);

/* Intercala tres enteros de 21 bits en un codigo de 63 bits (x arriba). */
uint64_t morton_codificar(uint32_t ix, uint32_t iy, uint32_t iz);

/* Calcula el codigo Morton de cada particula dentro del cubo dado. */
void calcular_morton(const Sistema *s, const double minimo[3], double lado,
                     uint64_t *codigos);

/* Llena orden[0..n-1] con los indices de las particulas ordenados por
 * codigo Morton (en empate, por indice, para que sea reproducible).
 * No mueve las particulas: solo devuelve la permutacion. */
void ordenar_por_morton(int n, const uint64_t *codigos, int *orden);

#endif
