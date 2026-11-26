#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <stdint.h>
#include "reparto.h"
#include "morton.h"

#define MARGEN   1e-6   /* igual que caja_envolvente() en c/common/morton.c */
#define MUESTRAS 256    /* muestras por proceso para elegir los cortes      */

void caja_global(const Locales *L, double minimo[3], double *lado,
                 MPI_Comm comm)
{
    double mn[3] = { DBL_MAX, DBL_MAX, DBL_MAX };
    double mx[3] = { -DBL_MAX, -DBL_MAX, -DBL_MAX };
    for (int i = 0; i < L->n; i++)
        for (int k = 0; k < 3; k++) {
            double p = L->pos[3 * i + k];
            if (p < mn[k]) mn[k] = p;
            if (p > mx[k]) mx[k] = p;
        }
    MPI_Allreduce(MPI_IN_PLACE, mn, 3, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(MPI_IN_PLACE, mx, 3, MPI_DOUBLE, MPI_MAX, comm);

    double arista = 0.0;
    for (int k = 0; k < 3; k++)
        if (mx[k] - mn[k] > arista)
            arista = mx[k] - mn[k];
    if (arista == 0.0)
        arista = 1.0;
    double margen = MARGEN * arista;
    for (int k = 0; k < 3; k++)
        minimo[k] = mn[k] - margen;
    *lado = arista + 2.0 * margen;
}

/* Entero de 21 bits de una coordenada (como a_entero() en morton.c). */
static uint32_t a_entero(double p, double minimo, double lado)
{
    const double escala = (double)(1u << MORTON_BITS);
    double t = (p - minimo) / lado * escala;
    if (t < 0.0) t = 0.0;
    if (t > escala - 1.0) t = escala - 1.0;
    return (uint32_t)t;
}

/* ===== ESTRATEGIA DE PARTICIONADO: la clave que se ordena y se corta ===== */
static uint64_t clave_particula(const Locales *L, int i,
                                const double minimo[3], double lado,
                                Particion modo)
{
    const double *p = &L->pos[3 * i];
    switch (modo) {
    case PARTICION_BISECCION: {
        /* Solo la coordenada x: rebanadas perpendiculares al eje x. */
        double t = (p[0] - minimo[0]) / lado;
        if (t < 0.0) t = 0.0;
        if (t > 1.0) t = 1.0;
        return (uint64_t)(t * (double)(UINT64_C(1) << 62));
    }
    case PARTICION_MORTON:
    default:
        return morton_codificar(a_entero(p[0], minimo[0], lado),
                                a_entero(p[1], minimo[1], lado),
                                a_entero(p[2], minimo[2], lado));
    }
}

typedef struct { uint64_t clave; double peso; } Muestra;

static int comparar_muestras(const void *a, const void *b)
{
    uint64_t x = ((const Muestra *)a)->clave, y = ((const Muestra *)b)->clave;
    return (x > y) - (x < y);
}

static double peso(const Locales *L, int i, int usar_costo)
{
    return (usar_costo && L->costo[i] > 0.0) ? L->costo[i] : 1.0;
}

void repartir(Locales *L, const double minimo[3], double lado,
              Particion modo, int usar_costo, MPI_Comm comm)
{
    int P, yo;
    MPI_Comm_size(comm, &P);
    MPI_Comm_rank(comm, &yo);
    int n = L->n;

    /* 1. Claves locales, ordenadas (ordenar_por_morton ordena cualquier
     *    clave de 64 bits, con empates por indice). */
    uint64_t *clave = malloc((size_t)(n ? n : 1) * sizeof(uint64_t));
    int *orden = malloc((size_t)(n ? n : 1) * sizeof(int));
    for (int i = 0; i < n; i++)
        clave[i] = clave_particula(L, i, minimo, lado, modo);
    ordenar_por_morton(n, clave, orden);

    /* 2. Muestras locales en cuantiles de PESO: cada una representa
     *    W_local / MUESTRAS de trabajo. */
    double W = 0.0;
    for (int i = 0; i < n; i++)
        W += peso(L, i, usar_costo);
    int m = n > 0 ? MUESTRAS : 0;
    Muestra mias[MUESTRAS];
    double acum = 0.0;
    int j = 0;
    for (int k = 0; k < m; k++) {
        double objetivo = (k + 0.5) * W / m;
        while (j < n - 1 && acum + peso(L, orden[j], usar_costo) < objetivo) {
            acum += peso(L, orden[j], usar_costo);
            j++;
        }
        mias[k].clave = clave[orden[j]];
        mias[k].peso = W / m;
    }

    /* 3. Todos reciben todas las muestras y calculan los MISMOS cortes. */
    int *cuantas = malloc((size_t)P * sizeof(int));
    int *desp = malloc((size_t)P * sizeof(int));
    int bytes = m * (int)sizeof(Muestra);
    MPI_Allgather(&bytes, 1, MPI_INT, cuantas, 1, MPI_INT, comm);
    int total = 0;
    for (int r = 0; r < P; r++) {
        desp[r] = total;
        total += cuantas[r];
    }
    Muestra *todas = malloc((size_t)(total ? total : 1));
    MPI_Allgatherv(mias, bytes, MPI_BYTE, todas, cuantas, desp, MPI_BYTE, comm);
    int nm = total / (int)sizeof(Muestra);
    qsort(todas, (size_t)nm, sizeof(Muestra), comparar_muestras);

    double Wtot = 0.0;
    for (int k = 0; k < nm; k++)
        Wtot += todas[k].peso;
    uint64_t *corte = malloc((size_t)(P > 1 ? P - 1 : 1) * sizeof(uint64_t));
    int t = 1;
    acum = 0.0;
    for (int k = 0; k < nm && t < P; k++) {
        acum += todas[k].peso;
        while (t < P && acum >= t * Wtot / P)
            corte[t - 1] = todas[k].clave, t++;
    }
    for (; t < P; t++)
        corte[t - 1] = UINT64_MAX;

    /* 4. Destino de cada particula: el tramo de claves donde cae. Como
     *    estan ordenadas, los destinos son crecientes: un solo barrido. */
    ParticulaMPI *envio = malloc((size_t)(n ? n : 1) * sizeof(ParticulaMPI));
    int *env_bytes = calloc((size_t)P, sizeof(int));
    int *env_desp = malloc((size_t)P * sizeof(int));
    int d = 0;
    for (int q = 0; q < n; q++) {
        int i = orden[q];
        while (d < P - 1 && clave[i] >= corte[d])
            d++;
        locales_empacar(L, i, &envio[q]);
        env_bytes[d] += (int)sizeof(ParticulaMPI);
    }
    for (int r = 0, s = 0; r < P; r++) {
        env_desp[r] = s;
        s += env_bytes[r];
    }

    /* 5. Intercambio: primero cuantos bytes, luego las particulas. */
    int *rec_bytes = malloc((size_t)P * sizeof(int));
    int *rec_desp = malloc((size_t)P * sizeof(int));
    MPI_Alltoall(env_bytes, 1, MPI_INT, rec_bytes, 1, MPI_INT, comm);
    int rtotal = 0;
    for (int r = 0; r < P; r++) {
        rec_desp[r] = rtotal;
        rtotal += rec_bytes[r];
    }
    ParticulaMPI *recibo = malloc((size_t)(rtotal ? rtotal : 1));
    MPI_Alltoallv(envio, env_bytes, env_desp, MPI_BYTE,
                  recibo, rec_bytes, rec_desp, MPI_BYTE, comm);

    int nuevo_n = rtotal / (int)sizeof(ParticulaMPI);
    locales_reservar(L, nuevo_n);
    L->n = nuevo_n;
    for (int i = 0; i < nuevo_n; i++)
        locales_desempacar(L, i, &recibo[i]);

    free(clave); free(orden); free(cuantas); free(desp); free(todas);
    free(corte); free(envio); free(env_bytes); free(env_desp);
    free(rec_bytes); free(rec_desp); free(recibo);
}
