/*
 * locales.h
 * Particulas que posee UN proceso MPI. Cambian de dueno en cada paso
 * (reparto.c), por eso cada una lleva su id global (fila del archivo).
 */
#ifndef LOCALES_H
#define LOCALES_H

#include <stdint.h>
#include "particulas.h"

typedef struct {
    int n, cap;
    double  *pos;    /* 3*n */
    double  *vel;    /* 3*n */
    double  *acc;    /* 3*n */
    double  *masa;   /* n   */
    double  *costo;  /* n: interacciones del ultimo recorrido (peso del reparto) */
    int32_t *id;     /* n: indice global de la particula */
} Locales;

/* Particula empaquetada para enviarla por MPI (como MPI_BYTE). */
typedef struct {
    double pos[3], vel[3], masa, costo;
    int32_t id, relleno;
} ParticulaMPI;

void locales_reservar(Locales *L, int n);   /* cap >= n, conserva datos */
void locales_liberar(Locales *L);
void locales_empacar(const Locales *L, int i, ParticulaMPI *p);
void locales_desempacar(Locales *L, int i, const ParticulaMPI *p);

/* Vista Sistema (sin copiar) para usar las funciones del nucleo C. */
static inline Sistema locales_sistema(const Locales *L)
{
    Sistema s = { L->n, L->pos, L->vel, L->acc, L->masa };
    return s;
}

#endif
