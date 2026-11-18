/* tiempo.h - reloj de pared en segundos (clock_gettime monotonico). */
#ifndef TIEMPO_H
#define TIEMPO_H

#include <time.h>

static inline double ahora(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + 1e-9 * (double)t.tv_nsec;
}

#endif
