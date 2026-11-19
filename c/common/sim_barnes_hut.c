/*
 * sim_barnes_hut.c
 * Simulacion Barnes-Hut en C, un hilo. Base de OpenMP, MPI y CUDA.
 *
 * Uso (contrato comun del README):
 *   sim_barnes_hut <archivo> <N> <pasos> <dt> <theta> <eps> <hilos> [salida_acc]
 * hilos se acepta pero no se usa en esta version.
 *
 * En cada calculo de fuerzas: Morton -> ordenar -> arbol -> recorrido.
 * El arbol se reconstruye en cada paso porque las particulas se mueven.
 */
#include <stdio.h>
#include <stdlib.h>
#include "particulas.h"
#include "directa.h"
#include "leapfrog.h"
#include "morton.h"
#include "arbol.h"
#include "recorrido.h"
#include "tiempo.h"

#define N_MAXIMO_ENERGIA 5000   /* la energia es O(N^2): solo para N pequeno */

typedef struct {
    double theta, eps;
    ArbolPlano arbol;
    uint64_t *codigos;
    int *orden;
    double t_sort, t_build, t_walk;   /* acumulados por fase */
} ConfigBH;

/* Calcula s->acc con Barnes-Hut y acumula el tiempo de cada fase. */
static void fuerzas_barnes_hut(Sistema *s, void *config)
{
    ConfigBH *c = config;
    double minimo[3], lado;

    double t0 = ahora();
    caja_envolvente(s, minimo, &lado);
    calcular_morton(s, minimo, lado, c->codigos);
    ordenar_por_morton(s->n, c->codigos, c->orden);
    double t1 = ahora();
    construir_arbol(&c->arbol, s, c->orden, c->codigos, minimo, lado);
    double t2 = ahora();
    aceleraciones_arbol(&c->arbol, s, c->theta, c->eps);
    double t3 = ahora();

    c->t_sort  += t1 - t0;
    c->t_build += t2 - t1;
    c->t_walk  += t3 - t2;
}

int main(int argc, char **argv)
{
    if (argc < 8) {
        fprintf(stderr,
                "uso: %s <archivo> <N> <pasos> <dt> <theta> <eps> <hilos> [salida_acc]\n",
                argv[0]);
        return 1;
    }

    const char *archivo = argv[1];
    int n        = atoi(argv[2]);
    int pasos    = atoi(argv[3]);
    double dt    = atof(argv[4]);
    double theta = atof(argv[5]);
    double eps   = atof(argv[6]);
    /* argv[7] = hilos (esta version es secuencial) */
    const char *salida_acc = (argc > 8) ? argv[8] : "aceleraciones.txt";

    Sistema s;
    if (leer_particulas(archivo, n, &s) != 0)
        return 1;

    ConfigBH c = { .theta = theta, .eps = eps };
    c.arbol   = crear_arbol(2 * n);
    c.codigos = malloc((size_t)n * sizeof(uint64_t));
    c.orden   = malloc((size_t)n * sizeof(int));
    if (!c.codigos || !c.orden) {
        fprintf(stderr, "error: no hay memoria\n");
        return 1;
    }

    /* Aceleraciones iniciales: se escriben para verificar.py
     * y son las que usa el primer paso del leapfrog. */
    fuerzas_barnes_hut(&s, &c);
    if (escribir_aceleraciones(salida_acc, &s) != 0)
        return 1;

    int con_energia = (n <= N_MAXIMO_ENERGIA);
    double e_inicial = con_energia ? energia_total(&s, eps) : 0.0;

    double t_fuerzas = 0.0, t_integrar = 0.0;
    simular_leapfrog(&s, pasos, dt, fuerzas_barnes_hut, &c,
                     &t_fuerzas, &t_integrar);

    double t_total = c.t_sort + c.t_build + c.t_walk + t_integrar;
    printf("tiempo_build=%.6f\n", c.t_build);
    printf("tiempo_sort=%.6f\n", c.t_sort);
    printf("tiempo_walk=%.6f\n", c.t_walk);
    printf("tiempo_integrate=%.6f\n", t_integrar);
    printf("tiempo_total=%.6f\n", t_total);
    if (con_energia) {
        double e_final = energia_total(&s, eps);
        printf("energia_inicial=%.12f\n", e_inicial);
        printf("energia_final=%.12f\n", e_final);
    }

    free(c.codigos);
    free(c.orden);
    liberar_arbol(&c.arbol);
    liberar_sistema(&s);
    return 0;
}
