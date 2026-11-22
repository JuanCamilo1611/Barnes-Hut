/*
 * sim_mpi.c
 * Barnes-Hut distribuido con MPI.
 *
 * Uso (contrato comun del README):
 *   mpirun -np P sim_mpi <archivo> <N> <pasos> <dt> <theta> <eps> <hilos> [salida_acc]
 * El numero de procesos lo fija mpirun; <hilos> solo se compara con P.
 *
 * Variables de entorno (para cambiar la estrategia sin recompilar):
 *   BH_PARTICION = morton (defecto) | biseccion      -> reparto.c
 *   BH_BALANCE   = costo  (defecto) | particulas      -> peso del reparto
 *
 * En cada calculo de fuerzas, cada proceso:
 *   1. caja global (Allreduce) y reparto por claves (Alltoallv)   [reparto]
 *   2. Morton y orden de SUS particulas                            [sort]
 *   3. arbol local solo con SUS particulas (nucleo C)             [build]
 *   4. intercambio de LET con los demas procesos (Alltoallv)       [let]
 *   5. recorrido: arbol local + LET recibidos                      [walk]
 * El leapfrog reutiliza las fuerzas: una sola evaluacion por paso.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "particulas.h"
#include "directa.h"
#include "morton.h"
#include "arbol.h"
#include "locales.h"
#include "reparto.h"
#include "let.h"

#define N_MAXIMO_ENERGIA 5000   /* la energia es O(N^2): solo para N pequeno */

typedef struct {
    double theta, eps;
    Particion modo;
    int usar_costo;
    ArbolPlano arbol;
    LETRecibidos let;
    uint64_t *codigos;
    int *orden, cap;
    /* tiempos acumulados de ESTE proceso */
    double t_reparto, t_sort, t_build, t_let, t_walk, t_integrar;
    long let_recibidos_max, nodos_locales_max;
} Config;

/* ------------------------------------------------------------------ */

static void fuerzas(Locales *L, Config *c, MPI_Comm comm)
{
    double minimo[3], lado;

    double t0 = MPI_Wtime();
    caja_global(L, minimo, &lado, comm);
    repartir(L, minimo, lado, c->modo, c->usar_costo, comm);
    double t1 = MPI_Wtime();

    Sistema s = locales_sistema(L);
    if (L->n > c->cap) {
        c->cap = 2 * L->n;
        c->codigos = realloc(c->codigos, (size_t)c->cap * sizeof(uint64_t));
        c->orden = realloc(c->orden, (size_t)c->cap * sizeof(int));
    }
    calcular_morton(&s, minimo, lado, c->codigos);
    ordenar_por_morton(L->n, c->codigos, c->orden);
    double t2 = MPI_Wtime();

    /* Mismo cubo raiz en todos: las celdas locales son celdas globales. */
    construir_arbol(&c->arbol, &s, c->orden, c->codigos, minimo, lado);
    double t3 = MPI_Wtime();

    intercambiar_let(&c->arbol, L, c->theta, &c->let, comm);
    double t4 = MPI_Wtime();

    for (int i = 0; i < L->n; i++) {
        double a[3] = { 0.0, 0.0, 0.0 };
        const double *p = &L->pos[3 * i];
        long costo = acumular_arbol(&c->arbol, p, i, c->theta, c->eps, a);
        for (int r = 0; r < c->let.P; r++)
            costo += acumular_arbol(&c->let.arbol[r], p, -1,
                                    c->theta, c->eps, a);
        for (int k = 0; k < 3; k++)
            L->acc[3 * i + k] = a[k];
        L->costo[i] = (double)costo;
    }
    double t5 = MPI_Wtime();

    c->t_reparto += t1 - t0;
    c->t_sort    += t2 - t1;
    c->t_build   += t3 - t2;
    c->t_let     += t4 - t3;
    c->t_walk    += t5 - t4;
    if (c->let.nodos_recibidos > c->let_recibidos_max)
        c->let_recibidos_max = c->let.nodos_recibidos;
    if (c->arbol.n_nodos > c->nodos_locales_max)
        c->nodos_locales_max = c->arbol.n_nodos;
}

/* v += (dt/2) a   */
static void media_patada(Locales *L, double dt)
{
    for (int k = 0; k < 3 * L->n; k++)
        L->vel[k] += 0.5 * dt * L->acc[k];
}

/* x += dt v */
static void desplazar(Locales *L, double dt)
{
    for (int k = 0; k < 3 * L->n; k++)
        L->pos[k] += dt * L->vel[k];
}

/* ------------------------------------------------------------------ */
/* Entrada/salida: solo el proceso 0 toca archivos.                   */

typedef struct {
    double pos[3], vel[3], acc[3], masa;
    int32_t id, relleno;
} Completa;

/* Junta todas las particulas en el proceso 0, en el orden del archivo. */
static int reunir(const Locales *L, int n_total, Sistema *g, MPI_Comm comm)
{
    int P, yo;
    MPI_Comm_size(comm, &P);
    MPI_Comm_rank(comm, &yo);

    Completa *mias = malloc((size_t)(L->n ? L->n : 1) * sizeof(Completa));
    for (int i = 0; i < L->n; i++) {
        for (int k = 0; k < 3; k++) {
            mias[i].pos[k] = L->pos[3 * i + k];
            mias[i].vel[k] = L->vel[3 * i + k];
            mias[i].acc[k] = L->acc[3 * i + k];
        }
        mias[i].masa = L->masa[i];
        mias[i].id = L->id[i];
    }
    int bytes = L->n * (int)sizeof(Completa);
    int *cuantos = NULL, *desp = NULL;
    Completa *todas = NULL;
    if (yo == 0) {
        cuantos = malloc((size_t)P * sizeof(int));
        desp = malloc((size_t)P * sizeof(int));
    }
    MPI_Gather(&bytes, 1, MPI_INT, cuantos, 1, MPI_INT, 0, comm);
    if (yo == 0) {
        int t = 0;
        for (int r = 0; r < P; r++) {
            desp[r] = t;
            t += cuantos[r];
        }
        todas = malloc((size_t)t);
    }
    MPI_Gatherv(mias, bytes, MPI_BYTE, todas, cuantos, desp, MPI_BYTE, 0, comm);
    if (yo == 0) {
        for (int i = 0; i < n_total; i++) {
            const Completa *q = &todas[i];
            int j = q->id;
            for (int k = 0; k < 3; k++) {
                g->pos[3 * j + k] = q->pos[k];
                g->vel[3 * j + k] = q->vel[k];
                g->acc[3 * j + k] = q->acc[k];
            }
            g->masa[j] = q->masa;
        }
        free(cuantos); free(desp); free(todas);
    }
    free(mias);
    return 0;
}

/* El proceso 0 lee el archivo y reparte bloques consecutivos; el primer
 * calculo de fuerzas ya los redistribuye por la estrategia elegida. */
static int distribuir_inicial(const char *archivo, int n, Locales *L,
                              MPI_Comm comm)
{
    int P, yo;
    MPI_Comm_size(comm, &P);
    MPI_Comm_rank(comm, &yo);

    int ok = 1;
    Sistema s = { 0 };
    ParticulaMPI *todas = NULL;
    if (yo == 0) {
        ok = (leer_particulas(archivo, n, &s) == 0);
        if (ok) {
            todas = malloc((size_t)n * sizeof(ParticulaMPI));
            for (int i = 0; i < n; i++) {
                for (int k = 0; k < 3; k++) {
                    todas[i].pos[k] = s.pos[3 * i + k];
                    todas[i].vel[k] = s.vel[3 * i + k];
                }
                todas[i].masa = s.masa[i];
                todas[i].costo = 0.0;
                todas[i].id = i;
                todas[i].relleno = 0;
            }
            liberar_sistema(&s);
        }
    }
    MPI_Bcast(&ok, 1, MPI_INT, 0, comm);
    if (!ok)
        return -1;

    int *bytes = malloc((size_t)P * sizeof(int));
    int *desp = malloc((size_t)P * sizeof(int));
    for (int r = 0, t = 0; r < P; r++) {
        int cuantas = n / P + (r < n % P ? 1 : 0);
        bytes[r] = cuantas * (int)sizeof(ParticulaMPI);
        desp[r] = t;
        t += bytes[r];
    }
    int mias = bytes[yo] / (int)sizeof(ParticulaMPI);
    ParticulaMPI *buf = malloc((size_t)(mias ? mias : 1) * sizeof(ParticulaMPI));
    MPI_Scatterv(todas, bytes, desp, MPI_BYTE, buf, bytes[yo], MPI_BYTE, 0, comm);

    locales_reservar(L, mias);
    L->n = mias;
    for (int i = 0; i < mias; i++)
        locales_desempacar(L, i, &buf[i]);

    free(bytes); free(desp); free(buf); free(todas);
    return 0;
}

static double energia_global(const Locales *L, int n, double eps,
                             MPI_Comm comm)
{
    int yo;
    MPI_Comm_rank(comm, &yo);
    Sistema g = { 0 };
    if (yo == 0)
        g = crear_sistema(n);
    reunir(L, n, &g, comm);
    double e = 0.0;
    if (yo == 0) {
        e = energia_total(&g, eps);
        liberar_sistema(&g);
    }
    return e;
}

/* ------------------------------------------------------------------ */

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    int P, yo;
    MPI_Comm_size(comm, &P);
    MPI_Comm_rank(comm, &yo);

    if (argc < 8) {
        if (yo == 0)
            fprintf(stderr, "uso: mpirun -np P %s <archivo> <N> <pasos> <dt> "
                            "<theta> <eps> <hilos> [salida_acc]\n", argv[0]);
        MPI_Finalize();
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

    Config c;
    memset(&c, 0, sizeof(c));
    c.theta = theta;
    c.eps = eps;
    const char *ep = getenv("BH_PARTICION");
    c.modo = (ep && strcmp(ep, "biseccion") == 0) ? PARTICION_BISECCION
                                                  : PARTICION_MORTON;
    const char *eb = getenv("BH_BALANCE");
    c.usar_costo = !(eb && strcmp(eb, "particulas") == 0);

    if (yo == 0 && hilos != P)
        fprintf(stderr, "aviso: hilos=%d pero mpirun lanzo %d procesos; "
                        "se usan %d\n", hilos, P, P);

    Locales L;
    memset(&L, 0, sizeof(L));
    if (distribuir_inicial(archivo, n, &L, comm) != 0) {
        MPI_Finalize();
        return 1;
    }
    c.arbol = crear_arbol(2 * (n / P + 1));

    MPI_Barrier(comm);
    double t_inicio = MPI_Wtime();

    /* Aceleraciones iniciales: se escriben para verificar.py y son las
     * del primer medio paso del leapfrog. */
    fuerzas(&L, &c, comm);
    double t_pausa = MPI_Wtime();
    {
        Sistema g = { 0 };
        if (yo == 0)
            g = crear_sistema(n);
        reunir(&L, n, &g, comm);
        int ok = 1;
        if (yo == 0) {
            ok = (escribir_aceleraciones(salida_acc, &g) == 0);
            liberar_sistema(&g);
        }
        MPI_Bcast(&ok, 1, MPI_INT, 0, comm);
        if (!ok) {
            MPI_Finalize();
            return 1;
        }
    }
    int con_energia = (n <= N_MAXIMO_ENERGIA);
    double e_inicial = con_energia ? energia_global(&L, n, eps, comm) : 0.0;
    MPI_Barrier(comm);
    double t_fuera = MPI_Wtime() - t_pausa;   /* E/S: no cuenta en el total */

    for (int paso = 0; paso < pasos; paso++) {
        double a = MPI_Wtime();
        media_patada(&L, dt);
        desplazar(&L, dt);
        double b = MPI_Wtime();
        fuerzas(&L, &c, comm);
        double d = MPI_Wtime();
        media_patada(&L, dt);
        c.t_integrar += (b - a) + (MPI_Wtime() - d);
    }
    MPI_Barrier(comm);
    double t_total = MPI_Wtime() - t_inicio - t_fuera;

    /* Cada fase se reporta como el MAXIMO entre procesos (el mas lento
     * marca el ritmo); el walk tambien con su promedio (desbalance). */
    double mio[7] = { c.t_build, c.t_sort, c.t_walk, c.t_integrar,
                      c.t_reparto, c.t_let, c.t_walk };
    double mx[7], suma[7];
    MPI_Reduce(mio, mx, 7, MPI_DOUBLE, MPI_MAX, 0, comm);
    MPI_Reduce(mio, suma, 7, MPI_DOUBLE, MPI_SUM, 0, comm);
    long nodos[2] = { c.let_recibidos_max, c.nodos_locales_max }, nodos_mx[2];
    MPI_Reduce(nodos, nodos_mx, 2, MPI_LONG, MPI_MAX, 0, comm);

    double e_final = con_energia ? energia_global(&L, n, eps, comm) : 0.0;

    if (yo == 0) {
        printf("tiempo_build=%.6f\n", mx[0]);
        printf("tiempo_sort=%.6f\n", mx[1]);
        printf("tiempo_walk=%.6f\n", mx[2]);
        printf("tiempo_integrate=%.6f\n", mx[3]);
        printf("tiempo_total=%.6f\n", t_total);
        printf("tiempo_reparto=%.6f\n", mx[4]);
        printf("tiempo_let=%.6f\n", mx[5]);
        printf("desbalance_walk=%.3f\n",
               suma[6] > 0.0 ? mx[6] / (suma[6] / P) : 1.0);
        printf("procesos=%d\n", P);
        printf("nodos_arbol_local_max=%ld\n", nodos_mx[1]);
        printf("nodos_let_recibidos_max=%ld\n", nodos_mx[0]);
        if (con_energia) {
            printf("energia_inicial=%.12f\n", e_inicial);
            printf("energia_final=%.12f\n", e_final);
        }
    }

    liberar_let(&c.let);
    liberar_arbol(&c.arbol);
    free(c.codigos);
    free(c.orden);
    locales_liberar(&L);
    MPI_Finalize();
    return 0;
}
