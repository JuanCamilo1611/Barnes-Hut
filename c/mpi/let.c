#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "let.h"
#include "directa.h"     /* G_GRAVEDAD */
#include "recorrido.h"   /* TAM_PILA   */

/* Nodo empaquetado para enviarlo por MPI (como MPI_BYTE). */
typedef struct {
    double centro[3], lado, masa, com[3];
    int32_t hijos[8], part, relleno;
} NodoLET;

typedef struct { double lo[3], hi[3]; int vacia, relleno; } Caja;

/* ---------- criterio de apertura contra una CAJA (no un punto) ---------- */

/* Distancia al cuadrado del centro de masa a la caja (0 si esta dentro). */
static double dist2_caja(const double c[3], const Caja *b)
{
    double d2 = 0.0;
    for (int k = 0; k < 3; k++) {
        double d = 0.0;
        if (c[k] < b->lo[k]) d = b->lo[k] - c[k];
        else if (c[k] > b->hi[k]) d = c[k] - b->hi[k];
        d2 += d * d;
    }
    return d2;
}

/* La caja toca el cubo de la celda? (equivale a la prueba de
 * "la particula esta dentro de la celda" del recorrido secuencial). */
static int toca_cubo(const ArbolPlano *a, int n, const Caja *b)
{
    double h = 0.5 * a->lado[n];
    for (int k = 0; k < 3; k++) {
        double c = a->centro[3 * n + k];
        if (b->hi[k] < c - h || b->lo[k] > c + h)
            return 0;
    }
    return 1;
}

/* Se debe abrir la celda n para ALGUN punto de la caja? */
static int abrir_para_caja(const ArbolPlano *a, int n, const Caja *b,
                           double theta2)
{
    double s = a->lado[n];
    double d2 = dist2_caja(&a->com[3 * n], b);
    return d2 == 0.0 || s * s >= theta2 * d2 || toca_cubo(a, n, b);
}

/* ---------- armado del LET para una caja destino ---------- */

static void copiar_nodo(const ArbolPlano *a, int n, NodoLET *o)
{
    for (int k = 0; k < 3; k++) {
        o->centro[k] = a->centro[3 * n + k];
        o->com[k] = a->com[3 * n + k];
    }
    o->lado = a->lado[n];
    o->masa = a->masa[n];
    o->relleno = 0;
}

/* Escribe en out[] el LET para la caja b; devuelve cuantos nodos.
 * fin[n] = ultimo descendiente de n (rango contiguo en preorden).
 * nuevo[] es memoria de trabajo de tamano n_nodos. */
static int armar_let(const ArbolPlano *a, const int *fin, const Caja *b,
                     double theta, NodoLET *out, int *nuevo)
{
    const double theta2 = theta * theta;
    int m = 0, i = 0;

    while (i < a->n_nodos) {
        NodoLET *o = &out[m];
        copiar_nodo(a, i, o);
        nuevo[i] = m++;

        if (a->part[i] >= 0) {                     /* hoja: va tal cual   */
            o->part = PART_REMOTA;
            for (int k = 0; k < 8; k++) o->hijos[k] = -1;
            i++;
        } else if (!abrir_para_caja(a, i, b, theta2)) {
            o->part = PART_REMOTA;                 /* resumen: hoja       */
            for (int k = 0; k < 8; k++) o->hijos[k] = -1;
            i = fin[i] + 1;                        /* salta la rama entera */
        } else {
            o->part = -1;                          /* se abre: sus hijos  */
            for (int k = 0; k < 8; k++)            /* vienen a continuacion */
                o->hijos[k] = a->hijos[8 * i + k]; /* (indices viejos)    */
            i++;
        }
    }

    /* Reindexado: los hijos de los nodos abiertos siempre se enviaron. */
    for (int q = 0; q < m; q++)
        if (out[q].part == -1)
            for (int k = 0; k < 8; k++)
                if (out[q].hijos[k] >= 0)
                    out[q].hijos[k] = nuevo[out[q].hijos[k]];
    return m;
}

/* fin[n]: como un hijo siempre tiene indice mayor que su padre, basta
 * una pasada de atras hacia adelante. */
static void calcular_fin(const ArbolPlano *a, int *fin)
{
    for (int n = a->n_nodos - 1; n >= 0; n--) {
        fin[n] = n;
        for (int k = 0; k < 8; k++) {
            int h = a->hijos[8 * n + k];
            if (h >= 0 && fin[h] > fin[n])
                fin[n] = fin[h];
        }
    }
}

static void caja_local(const Locales *L, Caja *b)
{
    memset(b, 0, sizeof(*b));
    b->vacia = (L->n == 0);
    for (int k = 0; k < 3; k++) {
        b->lo[k] = DBL_MAX;
        b->hi[k] = -DBL_MAX;
    }
    for (int i = 0; i < L->n; i++)
        for (int k = 0; k < 3; k++) {
            double p = L->pos[3 * i + k];
            if (p < b->lo[k]) b->lo[k] = p;
            if (p > b->hi[k]) b->hi[k] = p;
        }
}

void liberar_let(LETRecibidos *rec)
{
    if (!rec->arbol)
        return;
    for (int r = 0; r < rec->P; r++)
        liberar_arbol(&rec->arbol[r]);
    free(rec->arbol);
    rec->arbol = NULL;
}

void intercambiar_let(const ArbolPlano *local, const Locales *L,
                      double theta, LETRecibidos *rec, MPI_Comm comm)
{
    int P, yo;
    MPI_Comm_size(comm, &P);
    MPI_Comm_rank(comm, &yo);

    /* 1. Todos conocen la caja de las particulas de todos. */
    Caja mia, *cajas = malloc((size_t)P * sizeof(Caja));
    caja_local(L, &mia);
    MPI_Allgather(&mia, sizeof(Caja), MPI_BYTE, cajas, sizeof(Caja),
                  MPI_BYTE, comm);

    /* 2. Un LET por destino, uno detras de otro en el mismo buffer. */
    int nn = local->n_nodos;
    int *fin = malloc((size_t)(nn ? nn : 1) * sizeof(int));
    int *nuevo = malloc((size_t)(nn ? nn : 1) * sizeof(int));
    calcular_fin(local, fin);

    NodoLET *envio = malloc((size_t)(P * (nn ? nn : 1)) * sizeof(NodoLET));
    int *env_bytes = calloc((size_t)P, sizeof(int));
    int *env_desp = malloc((size_t)P * sizeof(int));
    int usados = 0;
    for (int r = 0; r < P; r++) {
        env_desp[r] = usados * (int)sizeof(NodoLET);
        if (r == yo || cajas[r].vacia || nn == 0)
            continue;
        int m = armar_let(local, fin, &cajas[r], theta, &envio[usados], nuevo);
        env_bytes[r] = m * (int)sizeof(NodoLET);
        usados += m;
    }

    /* 3. Intercambio: cuantos bytes y luego los nodos. */
    int *rec_bytes = malloc((size_t)P * sizeof(int));
    int *rec_desp = malloc((size_t)P * sizeof(int));
    MPI_Alltoall(env_bytes, 1, MPI_INT, rec_bytes, 1, MPI_INT, comm);
    int rtotal = 0;
    for (int r = 0; r < P; r++) {
        rec_desp[r] = rtotal;
        rtotal += rec_bytes[r];
    }
    NodoLET *recibo = malloc((size_t)(rtotal ? rtotal : 1));
    MPI_Alltoallv(envio, env_bytes, env_desp, MPI_BYTE,
                  recibo, rec_bytes, rec_desp, MPI_BYTE, comm);

    /* 4. Cada LET recibido se vuelve un ArbolPlano del contrato. */
    liberar_let(rec);
    rec->P = P;
    rec->arbol = calloc((size_t)P, sizeof(ArbolPlano));
    rec->nodos_enviados = usados;
    rec->nodos_recibidos = rtotal / (long)sizeof(NodoLET);
    for (int r = 0; r < P; r++) {
        int m = rec_bytes[r] / (int)sizeof(NodoLET);
        const NodoLET *q = (const NodoLET *)((const char *)recibo + rec_desp[r]);
        ArbolPlano *t = &rec->arbol[r];
        *t = crear_arbol(m);
        t->n_nodos = m;
        for (int n = 0; n < m; n++) {
            for (int k = 0; k < 3; k++) {
                t->centro[3 * n + k] = q[n].centro[k];
                t->com[3 * n + k] = q[n].com[k];
            }
            t->lado[n] = q[n].lado;
            t->masa[n] = q[n].masa;
            t->part[n] = q[n].part;
            for (int k = 0; k < 8; k++)
                t->hijos[8 * n + k] = q[n].hijos[k];
        }
    }

    free(cajas); free(fin); free(nuevo); free(envio);
    free(env_bytes); free(env_desp); free(rec_bytes); free(rec_desp);
    free(recibo);
}

/* ---------- recorrido (igual que c/common/recorrido.c, pero acumula) ---------- */

long acumular_arbol(const ArbolPlano *t, const double p[3], int excluir,
                    double theta, double eps, double acc[3])
{
    if (t->n_nodos == 0)
        return 0;
    const double theta2 = theta * theta, eps2 = eps * eps;
    double ax = 0.0, ay = 0.0, az = 0.0;
    long interacciones = 0;
    int pila[TAM_PILA];
    int tope = 0;
    pila[tope++] = 0;

    while (tope > 0) {
        int n = pila[--tope];
        double m = t->masa[n];
        if (m == 0.0)
            continue;
        double dx = t->com[3 * n + 0] - p[0];
        double dy = t->com[3 * n + 1] - p[1];
        double dz = t->com[3 * n + 2] - p[2];
        double d2 = dx * dx + dy * dy + dz * dz;
        int es_hoja = (t->part[n] >= 0);
        if (es_hoja && t->part[n] == excluir)
            continue;

        int abrir = 0;
        if (!es_hoja) {
            double s = t->lado[n];
            abrir = (s * s >= theta2 * d2);
            if (!abrir) {
                double h = 0.5 * s;
                abrir = fabs(p[0] - t->centro[3 * n + 0]) <= h &&
                        fabs(p[1] - t->centro[3 * n + 1]) <= h &&
                        fabs(p[2] - t->centro[3 * n + 2]) <= h;
            }
        }

        if (abrir) {
            for (int o = 7; o >= 0; o--) {
                int h = t->hijos[8 * n + o];
                if (h >= 0)
                    pila[tope++] = h;
            }
        } else {
            double r2 = d2 + eps2;
            double f = G_GRAVEDAD * m / (r2 * sqrt(r2));
            ax += f * dx;
            ay += f * dy;
            az += f * dz;
            interacciones++;
        }
    }
    acc[0] += ax;
    acc[1] += ay;
    acc[2] += az;
    return interacciones;
}
