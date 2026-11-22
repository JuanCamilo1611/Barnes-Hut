/*
 * prueba_morton.c
 * Comprueba los codigos Morton antes de construir el octree.
 *
 * Uso:  prueba_morton <archivo> <N>
 *
 * Pruebas:
 *   1. Codigos conocidos: x, y, z por separado y la esquina maxima.
 *   2. Tras ordenar, los codigos quedan de menor a mayor.
 *   3. Los 3 bits mas altos de cada codigo son el octante del contrato
 *      respecto al centro del cubo raiz.
 */
#include <stdio.h>
#include <stdlib.h>
#include "particulas.h"
#include "morton.h"
#include "tiempo.h"

static int fallos = 0;

static void revisar(int ok, const char *que)
{
    printf("[%s] %s\n", ok ? " OK " : "FALLA", que);
    if (!ok) fallos++;
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        fprintf(stderr, "uso: %s <archivo> <N>\n", argv[0]);
        return 1;
    }
    int n = atoi(argv[2]);

    /* 1. Codigos conocidos */
    uint32_t m = (1u << MORTON_BITS) - 1;
    revisar(morton_codificar(1, 0, 0) == 4, "x=1 -> codigo 4 (bit de x arriba)");
    revisar(morton_codificar(0, 1, 0) == 2, "y=1 -> codigo 2");
    revisar(morton_codificar(0, 0, 1) == 1, "z=1 -> codigo 1");
    revisar(morton_codificar(m, m, m) == (1ULL << 63) - 1,
            "esquina maxima -> 2^63 - 1 (usa los 63 bits)");

    Sistema s;
    if (leer_particulas(argv[1], n, &s) != 0)
        return 1;

    double minimo[3], lado;
    uint64_t *codigos = malloc((size_t)n * sizeof(uint64_t));
    int *orden = malloc((size_t)n * sizeof(int));

    double t0 = ahora();
    caja_envolvente(&s, minimo, &lado);
    calcular_morton(&s, minimo, lado, codigos);
    double t1 = ahora();
    ordenar_por_morton(n, codigos, orden);
    double t2 = ahora();

    /* 2. Orden creciente */
    int creciente = 1;
    for (int i = 1; i < n; i++)
        if (codigos[orden[i - 1]] > codigos[orden[i]])
            creciente = 0;
    revisar(creciente, "codigos ordenados de menor a mayor");

    /* 3. Octante de la raiz = 3 bits mas altos */
    int coinciden = 1;
    for (int i = 0; i < n; i++) {
        double c[3];
        for (int k = 0; k < 3; k++)
            c[k] = minimo[k] + lado / 2.0;
        int octante = 4 * (s.pos[3 * i + 0] >= c[0])
                    + 2 * (s.pos[3 * i + 1] >= c[1])
                    + 1 * (s.pos[3 * i + 2] >= c[2]);
        if ((int)(codigos[i] >> 60) != octante)
            coinciden = 0;
    }
    revisar(coinciden, "3 bits mas altos = octante del contrato");

    printf("cubo: min=(%.4f, %.4f, %.4f) lado=%.6f\n",
           minimo[0], minimo[1], minimo[2], lado);
    printf("tiempo_morton=%.6f\n", t1 - t0);
    printf("tiempo_sort=%.6f\n", t2 - t1);
    printf("primeras 5 en orden Morton:\n");
    for (int i = 0; i < 5 && i < n; i++)
        printf("  particula %6d  codigo %20llu\n",
               orden[i], (unsigned long long)codigos[orden[i]]);

    free(codigos);
    free(orden);
    liberar_sistema(&s);

    printf(fallos == 0 ? "\nTodas las pruebas pasaron.\n"
                       : "\n%d prueba(s) fallaron.\n", fallos);
    return fallos != 0;
}
