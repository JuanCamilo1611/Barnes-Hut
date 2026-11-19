/*
 * prueba_arbol.c
 * Comprueba que el arbol cumple las invariantes del contrato.
 *
 * Uso:  prueba_arbol <archivo> <N>
 *
 * Pruebas:
 *   1. Preorden: todo hijo tiene indice mayor que su padre, y el
 *      subarbol de cada nodo es un rango contiguo de indices.
 *   2. Hojas e internos: part >= 0 solo en hojas (sin hijos).
 *   3. Cada particula esta en exactamente una hoja (o en un bucket).
 *   4. La masa de la raiz es la masa total y su centro de masa es el
 *      centro de masa del sistema.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "particulas.h"
#include "morton.h"
#include "arbol.h"
#include "tiempo.h"

static int fallos = 0;

static void revisar(int ok, const char *que)
{
    printf("[%s] %s\n", ok ? " OK " : "FALLA", que);
    if (!ok) fallos++;
}

/* Ultimo indice del subarbol de n (recorre el hijo de mayor indice). */
static int ultimo_descendiente(const ArbolPlano *a, int n)
{
    for (;;) {
        int mayor = -1;
        for (int o = 0; o < 8; o++)
            if (a->hijos[8 * n + o] > mayor)
                mayor = a->hijos[8 * n + o];
        if (mayor < 0)
            return n;
        n = mayor;
    }
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        fprintf(stderr, "uso: %s <archivo> <N>\n", argv[0]);
        return 1;
    }
    int n = atoi(argv[2]);
    Sistema s;
    if (leer_particulas(argv[1], n, &s) != 0)
        return 1;

    double minimo[3], lado;
    uint64_t *codigos = malloc((size_t)n * sizeof(uint64_t));
    int *orden = malloc((size_t)n * sizeof(int));
    ArbolPlano a = crear_arbol(2 * n);

    double t0 = ahora();
    caja_envolvente(&s, minimo, &lado);
    calcular_morton(&s, minimo, lado, codigos);
    double t1 = ahora();
    ordenar_por_morton(n, codigos, orden);
    double t2 = ahora();
    construir_arbol(&a, &s, orden, codigos, minimo, lado);
    double t3 = ahora();

    /* 1. Preorden */
    int hijo_mayor = 1, contiguo = 1;
    for (int v = 0; v < a.n_nodos; v++) {
        int esperado = v + 1;   /* el primer hijo va justo despues */
        for (int o = 0; o < 8; o++) {
            int h = a.hijos[8 * v + o];
            if (h < 0) continue;
            if (h <= v) hijo_mayor = 0;
            if (h != esperado) contiguo = 0;
            esperado = ultimo_descendiente(&a, h) + 1;
        }
    }
    revisar(hijo_mayor, "todo hijo tiene indice mayor que su padre");
    revisar(contiguo, "cada subarbol es un rango contiguo (preorden)");

    /* 2. Hojas e internos */
    int bien_formado = 1;
    for (int v = 0; v < a.n_nodos; v++) {
        int tiene_hijos = 0;
        for (int o = 0; o < 8; o++)
            if (a.hijos[8 * v + o] >= 0) tiene_hijos = 1;
        if ((a.part[v] >= 0) == tiene_hijos) bien_formado = 0;
    }
    revisar(bien_formado, "part >= 0 solo en hojas; los internos tienen hijos");

    /* 3. Cada particula en una sola hoja */
    int *veces = calloc((size_t)n, sizeof(int));
    int hojas = 0, buckets = 0;
    for (int v = 0; v < a.n_nodos; v++)
        if (a.part[v] >= 0) {
            hojas++;
            veces[a.part[v]]++;
            if (fabs(a.masa[v] - s.masa[a.part[v]]) > 1e-15) buckets++;
        }
    int una_vez = 1;
    for (int i = 0; i < n; i++)
        if (veces[i] > 1) una_vez = 0;
    revisar(una_vez && hojas + 0 <= n, "ninguna particula aparece en dos hojas");
    revisar(hojas == n || buckets > 0,
            "hojas = N (o hay buckets en la profundidad maxima)");

    /* 4. Masa y centro de masa de la raiz */
    double mt = 0.0, cm[3] = {0, 0, 0};
    for (int i = 0; i < n; i++) {
        mt += s.masa[i];
        for (int k = 0; k < 3; k++) cm[k] += s.masa[i] * s.pos[3 * i + k];
    }
    double err_cm = 0.0;
    for (int k = 0; k < 3; k++)
        err_cm += fabs(cm[k] / mt - a.com[k]);
    revisar(fabs(mt - a.masa[0]) < 1e-9, "masa de la raiz = masa total");
    revisar(err_cm < 1e-9, "centro de masa de la raiz = centro de masa total");

    printf("nodos=%d  hojas=%d  buckets=%d  profundidad=%d\n",
           a.n_nodos, hojas, buckets, profundidad_arbol(&a));
    printf("tiempo_morton=%.6f\ntiempo_sort=%.6f\ntiempo_build=%.6f\n",
           t1 - t0, t2 - t1, t3 - t2);

    free(veces); free(codigos); free(orden);
    liberar_arbol(&a);
    liberar_sistema(&s);
    printf(fallos == 0 ? "\nTodas las pruebas pasaron.\n"
                       : "\n%d prueba(s) fallaron.\n", fallos);
    return fallos != 0;
}
