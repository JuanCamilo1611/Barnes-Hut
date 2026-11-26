#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "arbol_omp.h"

static void *reservar(void *p, size_t bytes)
{
    void *q = realloc(p, bytes);
    if (!q) {
        fprintf(stderr, "error: no hay memoria (arbol_omp)\n");
        exit(1);
    }
    return q;
}

ConstructorOMP crear_constructor(int tareas_por_hilo)
{
    ConstructorOMP k;
    memset(&k, 0, sizeof(k));
    k.tareas_por_hilo = (tareas_por_hilo > 0) ? tareas_por_hilo : 8;
    return k;
}

void liberar_constructor(ConstructorOMP *k)
{
    for (int j = 0; j < k->n_locales; j++)
        liberar_arbol(&k->locales[j]);
    free(k->locales);
    free(k->elem);
    free(k->tareas);
    memset(k, 0, sizeof(*k));
}

/* Igual que en arbol.c: octante de un codigo en la profundidad d. */
static inline int octante_morton(uint64_t codigo, int d)
{
    return (int)((codigo >> (60 - 3 * d)) & 7);
}

static int nuevo_elemento(ConstructorOMP *k)
{
    if (k->n_elem == k->cap_elem) {
        k->cap_elem = k->cap_elem ? 2 * k->cap_elem : 256;
        k->elem = reservar(k->elem, (size_t)k->cap_elem * sizeof(Elemento));
    }
    return k->n_elem++;
}

static void agregar_tarea(ConstructorOMP *k, int e)
{
    if (k->n_tareas == k->cap_tareas) {
        k->cap_tareas = k->cap_tareas ? 2 * k->cap_tareas : 256;
        k->tareas = reservar(k->tareas, (size_t)k->cap_tareas * sizeof(int));
    }
    k->elem[e].tarea = k->n_tareas;
    k->tareas[k->n_tareas++] = e;
}

/* Paso 1: esqueleto en preorden (mismo recorrido que construir_rango). */
static int esqueleto(ConstructorOMP *k, const int *orden,
                     const uint64_t *codigos, int ini, int fin,
                     const double centro[3], double lado, int prof,
                     int umbral)
{
    int e = nuevo_elemento(k);
    Elemento *el = &k->elem[e];
    el->ini = ini;  el->fin = fin;  el->prof = prof;  el->lado = lado;
    for (int d = 0; d < 3; d++)
        el->centro[d] = centro[d];
    for (int o = 0; o < 8; o++)
        el->hijos[o] = -1;

    /* Rango pequeno (o tope de profundidad): se vuelve una tarea. */
    if (fin - ini <= umbral || prof == MAX_PROFUNDIDAD) {
        el->es_tarea = 1;
        agregar_tarea(k, e);
        return e;
    }

    el->es_tarea = 0;
    int j = ini;
    double mitad = lado / 4.0;
    for (int oct = 0; oct < 8 && j < fin; oct++) {
        int m = j;
        while (m < fin && octante_morton(codigos[orden[m]], prof) == oct)
            m++;
        if (m == j)
            continue;
        double c[3] = {
            centro[0] + ((oct & 4) ? mitad : -mitad),
            centro[1] + ((oct & 2) ? mitad : -mitad),
            centro[2] + ((oct & 1) ? mitad : -mitad),
        };
        int hijo = esqueleto(k, orden, codigos, j, m, c, lado / 2.0,
                             prof + 1, umbral);
        k->elem[e].hijos[oct] = hijo;   /* elem pudo moverse: por indice */
        j = m;
    }
    return e;
}

/* Copia el subarbol local L al arbol final a partir del nodo 'off'. */
static void copiar_subarbol(ArbolPlano *a, const ArbolPlano *L, int off)
{
    int n = L->n_nodos;
    memcpy(&a->centro[3 * off], L->centro, 3 * (size_t)n * sizeof(double));
    memcpy(&a->com[3 * off],    L->com,    3 * (size_t)n * sizeof(double));
    memcpy(&a->lado[off],       L->lado,       (size_t)n * sizeof(double));
    memcpy(&a->masa[off],       L->masa,       (size_t)n * sizeof(double));
    memcpy(&a->part[off],       L->part,       (size_t)n * sizeof(int32_t));
    for (int q = 0; q < 8 * n; q++) {
        int h = L->hijos[q];
        a->hijos[8 * off + q] = (h >= 0) ? h + off : -1;  /* reindexar */
    }
}

void construir_arbol_omp(ArbolPlano *a, ConstructorOMP *k, const Sistema *s,
                         const int *orden, const uint64_t *codigos,
                         const double minimo[3], double lado)
{
    /* Granularidad: unas 'tareas_por_hilo' tareas por hilo, para que
     * schedule(dynamic) pueda equilibrar subarboles de distinto tamano. */
    int hilos = omp_get_max_threads();
    int umbral = s->n / (k->tareas_por_hilo * hilos);
    if (umbral < 1)
        umbral = 1;

    /* 1. Esqueleto (secuencial). */
    k->n_elem = 0;
    k->n_tareas = 0;
    double centro[3];
    for (int d = 0; d < 3; d++)
        centro[d] = minimo[d] + lado / 2.0;
    esqueleto(k, orden, codigos, 0, s->n, centro, lado, 0, umbral);

    if (k->n_locales < k->n_tareas) {
        k->locales = reservar(k->locales,
                              (size_t)k->n_tareas * sizeof(ArbolPlano));
        for (int j = k->n_locales; j < k->n_tareas; j++)
            k->locales[j] = crear_arbol(64);
        k->n_locales = k->n_tareas;
    }

    /* 2. Cada tarea es un subarbol independiente: en paralelo y sin
     *    candados (cada hilo escribe solo en su arbol local). */
    #pragma omp parallel for schedule(dynamic, 1)
    for (int j = 0; j < k->n_tareas; j++) {
        const Elemento *el = &k->elem[k->tareas[j]];
        construir_subarbol(&k->locales[j], s, orden, codigos, el->ini,
                           el->fin, el->centro, el->lado, el->prof);
    }

    /* 3. Donde empieza cada elemento en el arbol final (suma prefija en
     *    preorden: un nodo ocupa 1, una tarea ocupa su subarbol). */
    int pos = 0;
    for (int e = 0; e < k->n_elem; e++) {
        k->elem[e].inicio = pos;
        pos += k->elem[e].es_tarea ? k->locales[k->elem[e].tarea].n_nodos : 1;
    }
    arbol_reservar(a, pos);
    a->n_nodos = pos;

    #pragma omp parallel for schedule(dynamic, 1)
    for (int j = 0; j < k->n_tareas; j++)
        copiar_subarbol(a, &k->locales[j], k->elem[k->tareas[j]].inicio);

    /* 4. Nodos del esqueleto, de atras hacia adelante: sus hijos (que van
     *    despues en preorden) ya tienen masa y centro de masa. */
    for (int e = k->n_elem - 1; e >= 0; e--) {
        const Elemento *el = &k->elem[e];
        if (el->es_tarea)
            continue;
        int n = el->inicio;
        double m = 0.0, mom[3] = {0.0, 0.0, 0.0};
        for (int d = 0; d < 3; d++)
            a->centro[3 * n + d] = el->centro[d];
        a->lado[n] = el->lado;
        a->part[n] = -1;
        for (int o = 0; o < 8; o++) {
            int h = (el->hijos[o] >= 0) ? k->elem[el->hijos[o]].inicio : -1;
            a->hijos[8 * n + o] = h;
            if (h < 0)
                continue;
            m += a->masa[h];
            for (int d = 0; d < 3; d++)
                mom[d] += a->masa[h] * a->com[3 * h + d];
        }
        a->masa[n] = m;
        for (int d = 0; d < 3; d++)
            a->com[3 * n + d] = (m > 0.0) ? mom[d] / m : el->centro[d];
    }
}