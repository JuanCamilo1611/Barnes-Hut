/*
 * leapfrog.h
 * Integrador leapfrog (patada - desplazamiento - patada) que REUTILIZA
 * las fuerzas: las aceleraciones del final de un paso son las del
 * inicio del siguiente, asi que solo se calculan una vez por paso.
 *
 * El calculo de fuerzas se pasa como funcion, para usar el mismo
 * integrador con la suma directa, con Barnes-Hut o con OpenMP.
 */
#ifndef LEAPFROG_H
#define LEAPFROG_H

#include "particulas.h"

/* Firma de cualquier funcion que llene s->acc. */
typedef void (*FuncionFuerzas)(Sistema *s, void *config);

/* Ejecuta 'pasos' pasos de tamano dt.
 * PRECONDICION: s->acc ya debe tener las aceleraciones iniciales.
 * Si tiempo_fuerzas no es NULL, acumula ahi los segundos gastados
 * calculando fuerzas; el resto del tiempo es integracion. */
void simular_leapfrog(Sistema *s, int pasos, double dt,
                      FuncionFuerzas calcular, void *config,
                      double *tiempo_fuerzas, double *tiempo_integrar);

#endif
