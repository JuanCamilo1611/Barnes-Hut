/*
 * prueba_arbol_omp.c
 * Comprueba que el arbol construido en paralelo es IDENTICO, nodo por
 * nodo, al arbol secuencial del nucleo (mismo preorden, mismos valores).
 *
 * Uso:  prueba_arbol_omp <archivo> <N> <hilos>
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "particulas.h"
#include "morton.h"
#include "arbol.h"
#include "arbol_omp.h"
#include "tiempo.h"

static int iguales(const void *x, const void *y, size_t bytes)
{
    return memcmp(x, y, bytes) == 0;
}

int main(int argc, char **argv)
{
    if (argc < 4) {
        fprintf(stderr, "uso: %s <archivo> <N> <hilos>\n", argv[0]);
        return 1;
    }
    int n = atoi(argv[2]);
    omp_set_num_threads(atoi(argv[3]));

    Sistema s;
    if (leer_particulas(argv[1], n, &s) != 0)
        return 1;

    double minimo[3], lado;
    uint64_t *codigos = malloc((size_t)n * sizeof(uint64_t));
    int *orden = malloc((size_t)n * sizeof(int));
    caja_envolvente(&s, minimo, &lado);
    calcular_morton(&s, minimo, lado, codigos);
    ordenar_por_morton(n, codigos, orden);

    ArbolPlano sec = crear_arbol(2 * n);
    double t0 = ahora();
    construir_arbol(&sec, &s, orden, codigos, minimo, lado);
    double t1 = ahora();

    ArbolPlano par = crear_arbol(16);
    ConstructorOMP k = crear_constructor(8);
    construir_arbol_omp(&par, &k, &s, orden, codigos, minimo, lado); /* calentar */
    double t2 = ahora();
    construir_arbol_omp(&par, &k, &s, orden, codigos, minimo, lado);
    double t3 = ahora();

    int m = sec.n_nodos;
    int ok = (par.n_nodos == m)
        && iguales(sec.centro, par.centro, 3 * (size_t)m * sizeof(double))
        && iguales(sec.lado,   par.lado,       (size_t)m * sizeof(double))
        && iguales(sec.masa,   par.masa,       (size_t)m * sizeof(double))
        && iguales(sec.com,    par.com,    3 * (size_t)m * sizeof(double))
        && iguales(sec.hijos,  par.hijos,  8 * (size_t)m * sizeof(int32_t))
        && iguales(sec.part,   par.part,       (size_t)m * sizeof(int32_t));

    printf("hilos=%d  tareas=%d  nodos secuencial=%d  nodos paralelo=%d\n",
           omp_get_max_threads(), k.n_tareas, sec.n_nodos, par.n_nodos);
    printf("tiempo_build_secuencial=%.6f\n", t1 - t0);
    printf("tiempo_build_paralelo=%.6f\n", t3 - t2);
    printf("%s\n", ok ? "OK: el arbol paralelo es identico al secuencial."
                      : "FALLA: los arboles son distintos.");

    liberar_constructor(&k);
    liberar_arbol(&sec);
    liberar_arbol(&par);
    liberar_sistema(&s);
    free(codigos);
    free(orden);
    return ok ? 0 : 1;
}