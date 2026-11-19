/*
 * sim_directa.c
 * Simulacion por suma directa O(N^2) en C (un hilo).
 * Sirve para comprobar que el nucleo C (lectura, leapfrog, energia)
 * coincide con la referencia en Python, antes de meter el octree.
 *
 * Uso (contrato comun del README):
 *   sim_directa <archivo> <N> <pasos> <dt> <theta> <eps> <hilos> [salida_acc]
 * theta e hilos se aceptan pero no se usan en esta version.
 */
#include <stdio.h>
#include <stdlib.h>
#include "particulas.h"
#include "directa.h"
#include "leapfrog.h"
#include "tiempo.h"

#define N_MAXIMO_ENERGIA 5000   /* la energia es O(N^2): solo para N pequeno */

/* Adaptador: el integrador pide una funcion (Sistema*, void*). */
static void fuerzas_directas(Sistema *s, void *config)
{
    double eps = *(double *)config;
    aceleraciones_directas(s, eps);
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
    /* argv[5] = theta (no se usa en suma directa) */
    double eps   = atof(argv[6]);
    /* argv[7] = hilos (esta version es secuencial) */
    const char *salida_acc = (argc > 8) ? argv[8] : "aceleraciones.txt";

    Sistema s;
    if (leer_particulas(archivo, n, &s) != 0)
        return 1;

    double t_inicio = ahora();

    /* Aceleraciones iniciales: se escriben para verificar.py
     * y son las que usa el primer paso del leapfrog. */
    aceleraciones_directas(&s, eps);
    double t_fuerzas = ahora() - t_inicio;
    double t_integrar = 0.0;

    if (escribir_aceleraciones(salida_acc, &s) != 0)
        return 1;

    int con_energia = (n <= N_MAXIMO_ENERGIA);
    double e_inicial = con_energia ? energia_total(&s, eps) : 0.0;

    simular_leapfrog(&s, pasos, dt, fuerzas_directas, &eps,
                     &t_fuerzas, &t_integrar);
    double t_total = t_fuerzas + t_integrar;

    /* Salida clave=valor (contrato del README). En suma directa no hay
     * arbol ni ordenamiento, asi que build y sort valen 0. */
    printf("tiempo_build=%.6f\n", 0.0);
    printf("tiempo_sort=%.6f\n", 0.0);
    printf("tiempo_walk=%.6f\n", t_fuerzas);
    printf("tiempo_integrate=%.6f\n", t_integrar);
    printf("tiempo_total=%.6f\n", t_total);
    if (con_energia) {
        double e_final = energia_total(&s, eps);
        printf("energia_inicial=%.12f\n", e_inicial);
        printf("energia_final=%.12f\n", e_final);
    }

    liberar_sistema(&s);
    return 0;
}
