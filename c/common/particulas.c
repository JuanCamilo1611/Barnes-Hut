#include <stdio.h>
#include <stdlib.h>
#include "particulas.h"

Sistema crear_sistema(int n)
{
    Sistema s;
    s.n = n;
    s.pos  = calloc(3 * (size_t)n, sizeof(double));
    s.vel  = calloc(3 * (size_t)n, sizeof(double));
    s.acc  = calloc(3 * (size_t)n, sizeof(double));
    s.masa = calloc((size_t)n, sizeof(double));
    if (!s.pos || !s.vel || !s.acc || !s.masa) {
        fprintf(stderr, "error: no hay memoria para %d particulas\n", n);
        exit(1);
    }
    return s;
}

void liberar_sistema(Sistema *s)
{
    free(s->pos);
    free(s->vel);
    free(s->acc);
    free(s->masa);
    s->pos = s->vel = s->acc = s->masa = NULL;
    s->n = 0;
}

int leer_particulas(const char *ruta, int n, Sistema *s)
{
    FILE *f = fopen(ruta, "r");
    if (!f) {
        fprintf(stderr, "error: no se pudo abrir %s\n", ruta);
        return -1;
    }

    *s = crear_sistema(n);

    for (int i = 0; i < n; i++) {
        double *p = &s->pos[3 * i];
        double *v = &s->vel[3 * i];
        int leidos = fscanf(f, "%lf %lf %lf %lf %lf %lf %lf",
                            &p[0], &p[1], &p[2],
                            &v[0], &v[1], &v[2],
                            &s->masa[i]);
        if (leidos != 7) {
            fprintf(stderr, "error: %s fila %d no tiene 7 valores "
                            "(el archivo tiene menos de N=%d filas?)\n",
                    ruta, i + 1, n);
            fclose(f);
            return -1;
        }
    }

    fclose(f);
    return 0;
}

int escribir_aceleraciones(const char *ruta, const Sistema *s)
{
    FILE *f = fopen(ruta, "w");
    if (!f) {
        fprintf(stderr, "error: no se pudo escribir %s\n", ruta);
        return -1;
    }
    for (int i = 0; i < s->n; i++) {
        fprintf(f, "%.15e %.15e %.15e\n",
                s->acc[3 * i], s->acc[3 * i + 1], s->acc[3 * i + 2]);
    }
    fclose(f);
    return 0;
}
