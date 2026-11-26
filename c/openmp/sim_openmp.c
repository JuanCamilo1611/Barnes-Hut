/*
 * sim_openmp.c
 * Simulacion Barnes-Hut en C con OpenMP (memoria compartida, hilos).
 *
 * Uso (contrato comun del README):
 *   sim_openmp <archivo> <N> <pasos> <dt> <theta> <eps> <hilos> [salida_acc]
 *
 * Que se paraleliza en cada calculo de fuerzas:
 *   - caja envolvente (reduccion min/max) y codigos Morton (parallel for);
 *   - construccion del arbol desde el orden Morton (arbol_omp.c);
 *   - recorrido: una particula por iteracion, en orden Morton, con
 *     schedule(runtime) -> por defecto dynamic,64 (ver REPARTO abajo);
 *   - integracion leapfrog (parallel for).
 * El ordenamiento (qsort) sigue siendo secuencial.
 *
 * Variables de entorno opcionales (para medir y para la sustentacion):
 *   OMP_SCHEDULE="static" | "dynamic,64" | "guided"   reparto del recorrido
 *   BH_ARBOL_SECUENCIAL=1   construye el arbol con la version de un hilo
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "particulas.h"
#include "directa.h"
#include "morton.h"
#include "arbol.h"
#include "recorrido.h"
#include "tiempo.h"
#include "arbol_omp.h"

#define N_MAXIMO_ENERGIA 5000   /* la energia es O(N^2): solo para N pequeno */

/* REPARTO del recorrido si no se define OMP_SCHEDULE: bloques de 64
 * particulas que los hilos van tomando a medida que terminan. Las
 * particulas en zonas densas abren mas nodos (cuestan mas), por eso un
 * reparto estatico deja hilos esperando. */
#define REPARTO_POR_DEFECTO omp_sched_dynamic
#define BLOQUE_POR_DEFECTO  64

/* Tareas de construccion por hilo (granularidad del arbol paralelo). */
#define TAREAS_POR_HILO 8

typedef struct {
    double theta, eps;
    int arbol_secuencial;
    ArbolPlano arbol;
    ConstructorOMP constructor;
    uint64_t *codigos;
    int *orden;
    double t_sort, t_build, t_walk;   /* acumulados por fase */
} ConfigOMP;

/* Caja envolvente en paralelo (reduccion de minimos y maximos). */
static void caja_envolvente_omp(const Sistema *s, double minimo[3],
                                double *lado)
{
    double xmin = s->pos[0], ymin = s->pos[1], zmin = s->pos[2];
    double xmax = xmin, ymax = ymin, zmax = zmin;

    #pragma omp parallel for reduction(min: xmin, ymin, zmin) \
                             reduction(max: xmax, ymax, zmax)
    for (int i = 0; i < s->n; i++) {
        double x = s->pos[3 * i], y = s->pos[3 * i + 1], z = s->pos[3 * i + 2];
        if (x < xmin) xmin = x;
        if (y < ymin) ymin = y;
        if (z < zmin) zmin = z;
        if (x > xmax) xmax = x;
        if (y > ymax) ymax = y;
        if (z > zmax) zmax = z;
    }

    /* Igual que caja_envolvente() del nucleo: cubo + margen. */
    double arista = fmax(xmax - xmin, fmax(ymax - ymin, zmax - zmin));
    if (arista == 0.0)
        arista = 1.0;
    double margen = 1e-6 * arista;
    minimo[0] = xmin - margen;
    minimo[1] = ymin - margen;
    minimo[2] = zmin - margen;
    *lado = arista + 2.0 * margen;
}

/* Igual que a_entero() de morton.c. */
static uint32_t a_entero(double p, double minimo, double lado)
{
    const double escala = (double)(1u << MORTON_BITS);
    double t = (p - minimo) / lado * escala;
    if (t < 0.0) t = 0.0;
    if (t > escala - 1.0) t = escala - 1.0;
    return (uint32_t)t;
}

/* Codigos Morton en paralelo: cada particula es independiente. */
static void calcular_morton_omp(const Sistema *s, const double minimo[3],
                                double lado, uint64_t *codigos)
{
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < s->n; i++) {
        uint32_t ix = a_entero(s->pos[3 * i + 0], minimo[0], lado);
        uint32_t iy = a_entero(s->pos[3 * i + 1], minimo[1], lado);
        uint32_t iz = a_entero(s->pos[3 * i + 2], minimo[2], lado);
        codigos[i] = morton_codificar(ix, iy, iz);
    }
}

/* Calcula s->acc con Barnes-Hut y acumula el tiempo de cada fase. */
static void fuerzas_omp(Sistema *s, ConfigOMP *c)
{
    double minimo[3], lado;

    double t0 = ahora();
    caja_envolvente_omp(s, minimo, &lado);
    calcular_morton_omp(s, minimo, lado, c->codigos);
    ordenar_por_morton(s->n, c->codigos, c->orden);
    double t1 = ahora();

    if (c->arbol_secuencial)
        construir_arbol(&c->arbol, s, c->orden, c->codigos, minimo, lado);
    else
        construir_arbol_omp(&c->arbol, &c->constructor, s, c->orden,
                            c->codigos, minimo, lado);
    double t2 = ahora();

    /* Recorrido en orden Morton: particulas vecinas visitan casi los
     * mismos nodos (mejor uso de cache). Cada hilo escribe solo
     * s->acc de sus particulas -> no hay condicion de carrera. */
    const int *orden = c->orden;
    #pragma omp parallel for schedule(runtime)
    for (int j = 0; j < s->n; j++)
        aceleracion_particula(&c->arbol, s, orden[j], c->theta, c->eps);
    double t3 = ahora();

    c->t_sort  += t1 - t0;
    c->t_build += t2 - t1;
    c->t_walk  += t3 - t2;
}

/* v = v + (dt/2) * a */
static void media_patada(Sistema *s, double dt)
{
    const double h = 0.5 * dt;
    #pragma omp parallel for schedule(static)
    for (int k = 0; k < 3 * s->n; k++)
        s->vel[k] += h * s->acc[k];
}

/* x = x + dt * v */
static void desplazar(Sistema *s, double dt)
{
    #pragma omp parallel for schedule(static)
    for (int k = 0; k < 3 * s->n; k++)
        s->pos[k] += dt * s->vel[k];
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
    int hilos    = atoi(argv[7]);
    const char *salida_acc = (argc > 8) ? argv[8] : "aceleraciones.txt";

    if (hilos > 0)
        omp_set_num_threads(hilos);
    if (getenv("OMP_SCHEDULE") == NULL)
        omp_set_schedule(REPARTO_POR_DEFECTO, BLOQUE_POR_DEFECTO);

    Sistema s;
    if (leer_particulas(archivo, n, &s) != 0)
        return 1;

    ConfigOMP c = { .theta = theta, .eps = eps };
    c.arbol_secuencial = (getenv("BH_ARBOL_SECUENCIAL") != NULL);
    c.arbol       = crear_arbol(2 * n);
    c.constructor = crear_constructor(TAREAS_POR_HILO);
    c.codigos = malloc((size_t)n * sizeof(uint64_t));
    c.orden   = malloc((size_t)n * sizeof(int));
    if (!c.codigos || !c.orden) {
        fprintf(stderr, "error: no hay memoria\n");
        return 1;
    }

    /* Aceleraciones iniciales: para verificar.py y para el primer paso. */
    fuerzas_omp(&s, &c);
    if (escribir_aceleraciones(salida_acc, &s) != 0)
        return 1;

    int con_energia = (n <= N_MAXIMO_ENERGIA);
    double e_inicial = con_energia ? energia_total(&s, eps) : 0.0;

    /* Leapfrog con reuso de fuerzas (igual que leapfrog.c del nucleo). */
    double t_integrar = 0.0;
    for (int paso = 0; paso < pasos; paso++) {
        double t0 = ahora();
        media_patada(&s, dt);
        desplazar(&s, dt);
        double t1 = ahora();
        fuerzas_omp(&s, &c);
        double t2 = ahora();
        media_patada(&s, dt);
        double t3 = ahora();
        t_integrar += (t1 - t0) + (t3 - t2);
    }

    double t_total = c.t_sort + c.t_build + c.t_walk + t_integrar;
    printf("hilos=%d\n", omp_get_max_threads());
    printf("arbol=%s\n", c.arbol_secuencial ? "secuencial" : "paralelo");
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
    liberar_constructor(&c.constructor);
    liberar_arbol(&c.arbol);
    liberar_sistema(&s);
    return 0;
}