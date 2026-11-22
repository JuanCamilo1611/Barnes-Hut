#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "locales.h"

static void *crecer(void *p, size_t bytes)
{
    void *q = realloc(p, bytes ? bytes : 1);
    if (!q) {
        fprintf(stderr, "error: no hay memoria para las particulas locales\n");
        exit(1);
    }
    return q;
}

void locales_reservar(Locales *L, int n)
{
    if (n <= L->cap)
        return;
    int c = n > 2 * L->cap ? n : 2 * L->cap;
    size_t z = (size_t)c;
    L->pos   = crecer(L->pos,   3 * z * sizeof(double));
    L->vel   = crecer(L->vel,   3 * z * sizeof(double));
    L->acc   = crecer(L->acc,   3 * z * sizeof(double));
    L->masa  = crecer(L->masa,      z * sizeof(double));
    L->costo = crecer(L->costo,     z * sizeof(double));
    L->id    = crecer(L->id,        z * sizeof(int32_t));
    L->cap = c;
}

void locales_liberar(Locales *L)
{
    free(L->pos); free(L->vel); free(L->acc);
    free(L->masa); free(L->costo); free(L->id);
    memset(L, 0, sizeof(*L));
}

void locales_empacar(const Locales *L, int i, ParticulaMPI *p)
{
    for (int k = 0; k < 3; k++) {
        p->pos[k] = L->pos[3 * i + k];
        p->vel[k] = L->vel[3 * i + k];
    }
    p->masa = L->masa[i];
    p->costo = L->costo[i];
    p->id = L->id[i];
    p->relleno = 0;
}

void locales_desempacar(Locales *L, int i, const ParticulaMPI *p)
{
    for (int k = 0; k < 3; k++) {
        L->pos[3 * i + k] = p->pos[k];
        L->vel[3 * i + k] = p->vel[k];
        L->acc[3 * i + k] = 0.0;
    }
    L->masa[i] = p->masa;
    L->costo[i] = p->costo;
    L->id[i] = p->id;
}
