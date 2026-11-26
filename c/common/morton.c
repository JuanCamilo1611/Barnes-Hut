#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "morton.h"
#define MARGEN 1e-6   /* fraccion de la arista que se agrega al cubo */

void caja_envolvente(const Sistema *s, double minimo[3], double *lado)
{
    double maximo[3];
    for (int k = 0; k < 3; k++) {
        minimo[k] = s->pos[k];
        maximo[k] = s->pos[k];
    }
    for (int i = 1; i < s->n; i++) {
        for (int k = 0; k < 3; k++) {
            double p = s->pos[3 * i + k];
            if (p < minimo[k]) minimo[k] = p;
            if (p > maximo[k]) maximo[k] = p;
        }
    }

    /* Un cubo (no una caja): la arista es la mayor de las tres. */
    double arista = 0.0;
    for (int k = 0; k < 3; k++)
        if (maximo[k] - minimo[k] > arista)
            arista = maximo[k] - minimo[k];
    if (arista == 0.0)
        arista = 1.0;   /* caso raro: todas las particulas en el mismo punto */

    double margen = MARGEN * arista;
    for (int k = 0; k < 3; k++)
        minimo[k] -= margen;
    *lado = arista + 2.0 * margen;
}

/* Separa los 21 bits de v dejando dos ceros entre cada uno:
 * ...b2 b1 b0  ->  ...b2 0 0 b1 0 0 b0   (truco de mascaras) */
static uint64_t separar_bits(uint32_t v)
{
    uint64_t x = v & 0x1fffff;                     /* 21 bits */
    x = (x | (x << 32)) & 0x1f00000000ffffULL;
    x = (x | (x << 16)) & 0x1f0000ff0000ffULL;
    x = (x | (x << 8))  & 0x100f00f00f00f00fULL;
    x = (x | (x << 4))  & 0x10c30c30c30c30c3ULL;
    x = (x | (x << 2))  & 0x1249249249249249ULL;
    return x;
}

uint64_t morton_codificar(uint32_t ix, uint32_t iy, uint32_t iz)
{
    return (separar_bits(ix) << 2) | (separar_bits(iy) << 1) | separar_bits(iz);
}

/* Convierte una coordenada a entero de 21 bits dentro del cubo. */
static uint32_t a_entero(double p, double minimo, double lado)
{
    const double escala = (double)(1u << MORTON_BITS);   /* 2^21 */
    double t = (p - minimo) / lado * escala;
    if (t < 0.0) t = 0.0;
    if (t > escala - 1.0) t = escala - 1.0;
    return (uint32_t)t;
}

void calcular_morton(const Sistema *s, const double minimo[3], double lado,
                     uint64_t *codigos)
{
    for (int i = 0; i < s->n; i++) {
        uint32_t ix = a_entero(s->pos[3 * i + 0], minimo[0], lado);
        uint32_t iy = a_entero(s->pos[3 * i + 1], minimo[1], lado);
        uint32_t iz = a_entero(s->pos[3 * i + 2], minimo[2], lado);
        codigos[i] = morton_codificar(ix, iy, iz);
    }
}

/* Par (codigo, indice): se ordenan los pares juntos, asi la funcion de
 * comparacion no necesita ninguna variable global y es segura con hilos. */
typedef struct {
    uint64_t codigo;
    int indice;
} ParMorton;

static int comparar_pares(const void *a, const void *b)
{
    const ParMorton *p = (const ParMorton *)a;
    const ParMorton *q = (const ParMorton *)b;
    if (p->codigo < q->codigo) return -1;
    if (p->codigo > q->codigo) return 1;
    return (p->indice > q->indice) - (p->indice < q->indice);  /* empate: por indice */
}

void ordenar_por_morton(int n, const uint64_t *codigos, int *orden)
{
    ParMorton *pares = malloc((size_t)n * sizeof(ParMorton));
    if (pares == NULL) {
        fprintf(stderr, "ordenar_por_morton: sin memoria\n");
        exit(1);
    }
    for (int i = 0; i < n; i++) {
        pares[i].codigo = codigos[i];
        pares[i].indice = i;
    }
    qsort(pares, (size_t)n, sizeof(ParMorton), comparar_pares);
    for (int i = 0; i < n; i++)
        orden[i] = pares[i].indice;
    free(pares);
}