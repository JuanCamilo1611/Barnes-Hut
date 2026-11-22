#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "arbol.h"
#include "morton.h"

static void *reservar(void *p, size_t bytes)
{
    void *q = realloc(p, bytes);
    if (!q) {
        fprintf(stderr, "error: no hay memoria para el arbol\n");
        exit(1);
    }
    return q;
}

static void ampliar(ArbolPlano *a, int nueva)
{
    size_t c = (size_t)nueva;
    a->centro        = reservar(a->centro,        3 * c * sizeof(double));
    a->lado          = reservar(a->lado,              c * sizeof(double));
    a->masa          = reservar(a->masa,              c * sizeof(double));
    a->com           = reservar(a->com,           3 * c * sizeof(double));
    a->hijos         = reservar(a->hijos,         8 * c * sizeof(int32_t));
    a->part          = reservar(a->part,              c * sizeof(int32_t));
    a->hojas_masa    = reservar(a->hojas_masa,        c * sizeof(double));
    a->hojas_momento = reservar(a->hojas_momento, 3 * c * sizeof(double));
    a->cap = nueva;
}

ArbolPlano crear_arbol(int cap_inicial)
{
    ArbolPlano a;
    memset(&a, 0, sizeof(a));
    ampliar(&a, cap_inicial > 16 ? cap_inicial : 16);
    return a;
}

void liberar_arbol(ArbolPlano *a)
{
    free(a->centro); free(a->lado); free(a->masa); free(a->com);
    free(a->hijos);  free(a->part);
    free(a->hojas_masa); free(a->hojas_momento);
    memset(a, 0, sizeof(*a));
}

/* Crea un nodo vacio al final del arreglo y devuelve su indice. */
static int nuevo_nodo(ArbolPlano *a, const double centro[3], double lado)
{
    if (a->n_nodos == a->cap)
        ampliar(a, 2 * a->cap);
    int n = a->n_nodos++;
    for (int k = 0; k < 3; k++) {
        a->centro[3 * n + k] = centro[k];
        a->com[3 * n + k] = 0.0;
        a->hojas_momento[3 * n + k] = 0.0;
    }
    a->lado[n] = lado;
    a->masa[n] = 0.0;
    a->hojas_masa[n] = 0.0;
    a->part[n] = -1;
    for (int k = 0; k < 8; k++)
        a->hijos[8 * n + k] = -1;
    return n;
}

/* Octante de un codigo Morton en la profundidad d (0 = hijos de la raiz):
 * son los 3 bits que estan 3*d posiciones por debajo de los mas altos. */
static inline int octante_morton(uint64_t codigo, int d)
{
    return (int)((codigo >> (60 - 3 * d)) & 7);
}

/* Construye el subarbol con las particulas orden[ini..fin-1] (todas caen
 * en este cubo). Crea el nodo ANTES que sus hijos -> preorden. */
static int construir_rango(ArbolPlano *a, const Sistema *s,
                           const int *orden, const uint64_t *codigos,
                           int ini, int fin, const double centro[3],
                           double lado, int prof)
{
    int nodo = nuevo_nodo(a, centro, lado);

    /* Hoja: una sola particula, o profundidad maxima (bucketing). */
    if (fin - ini == 1 || prof == MAX_PROFUNDIDAD) {
        a->part[nodo] = orden[ini];
        for (int j = ini; j < fin; j++) {
            int i = orden[j];
            double m = s->masa[i];
            a->hojas_masa[nodo] += m;
            for (int k = 0; k < 3; k++)
                a->hojas_momento[3 * nodo + k] += m * s->pos[3 * i + k];
        }
        return nodo;
    }

    /* Nodo interno: como las particulas estan ordenadas por Morton, las
     * de cada octante forman un bloque seguido. Se recorre una vez. */
    int j = ini;
    double mitad = lado / 4.0;
    for (int oct = 0; oct < 8 && j < fin; oct++) {
        int k = j;
        while (k < fin && octante_morton(codigos[orden[k]], prof) == oct)
            k++;
        if (k == j)
            continue;            /* octante vacio: hijos[] queda en -1 */

        double c[3] = {
            centro[0] + ((oct & 4) ? mitad : -mitad),
            centro[1] + ((oct & 2) ? mitad : -mitad),
            centro[2] + ((oct & 1) ? mitad : -mitad),
        };
        int hijo = construir_rango(a, s, orden, codigos, j, k, c,
                                   lado / 2.0, prof + 1);
        a->hijos[8 * nodo + oct] = hijo;   /* 'a' pudo crecer: se usa indice */
        j = k;
    }
    return nodo;
}

/* Masa y centro de masa: como todo hijo tiene indice mayor que su padre,
 * basta una pasada de atras hacia adelante, sin pila ni recursion. */
static void calcular_masa_com(ArbolPlano *a)
{
    for (int n = a->n_nodos - 1; n >= 0; n--) {
        double m = 0.0, mom[3] = {0.0, 0.0, 0.0};
        if (a->part[n] >= 0) {               /* hoja */
            m = a->hojas_masa[n];
            for (int k = 0; k < 3; k++)
                mom[k] = a->hojas_momento[3 * n + k];
        } else {                             /* interno */
            for (int o = 0; o < 8; o++) {
                int h = a->hijos[8 * n + o];
                if (h < 0) continue;
                m += a->masa[h];
                for (int k = 0; k < 3; k++)
                    mom[k] += a->masa[h] * a->com[3 * h + k];
            }
        }
        a->masa[n] = m;
        for (int k = 0; k < 3; k++)
            a->com[3 * n + k] = (m > 0.0) ? mom[k] / m : a->centro[3 * n + k];
    }
}

void construir_arbol(ArbolPlano *a, const Sistema *s,
                     const int *orden, const uint64_t *codigos,
                     const double minimo[3], double lado)
{
    a->n_nodos = 0;
    double centro[3];
    for (int k = 0; k < 3; k++)
        centro[k] = minimo[k] + lado / 2.0;
    construir_rango(a, s, orden, codigos, 0, s->n, centro, lado, 0);
    calcular_masa_com(a);
}

int profundidad_arbol(const ArbolPlano *a)
{
    /* La arista se divide a la mitad en cada nivel. */
    int maxp = 0;
    for (int n = 0; n < a->n_nodos; n++) {
        int p = 0;
        double l = a->lado[n];
        while (l < a->lado[0] * 0.75) { l *= 2.0; p++; }
        if (p > maxp) maxp = p;
    }
    return maxp;
}
