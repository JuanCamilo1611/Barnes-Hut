#include <math.h>
#include "directa.h"

void aceleraciones_directas(Sistema *s, double eps)
{
    const int n = s->n;
    const double eps2 = eps * eps;

    for (int i = 0; i < n; i++) {
        double xi = s->pos[3 * i];
        double yi = s->pos[3 * i + 1];
        double zi = s->pos[3 * i + 2];
        double ax = 0.0, ay = 0.0, az = 0.0;

        for (int j = 0; j < n; j++) {
            if (j == i) continue;               /* una particula no se atrae a si misma */
            double dx = s->pos[3 * j]     - xi;
            double dy = s->pos[3 * j + 1] - yi;
            double dz = s->pos[3 * j + 2] - zi;
            double r2 = dx * dx + dy * dy + dz * dz + eps2;
            double inv_r3 = 1.0 / (r2 * sqrt(r2)); /* 1 / r^3 */
            double f = G_GRAVEDAD * s->masa[j] * inv_r3;
            ax += f * dx;
            ay += f * dy;
            az += f * dz;
        }

        s->acc[3 * i]     = ax;
        s->acc[3 * i + 1] = ay;
        s->acc[3 * i + 2] = az;
    }
}

double energia_total(const Sistema *s, double eps)
{
    const int n = s->n;
    const double eps2 = eps * eps;
    double cinetica = 0.0;
    double potencial = 0.0;

    for (int i = 0; i < n; i++) {
        const double *v = &s->vel[3 * i];
        cinetica += 0.5 * s->masa[i] * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {     /* j > i: cada pareja una vez */
            double dx = s->pos[3 * j]     - s->pos[3 * i];
            double dy = s->pos[3 * j + 1] - s->pos[3 * i + 1];
            double dz = s->pos[3 * j + 2] - s->pos[3 * i + 2];
            double r = sqrt(dx * dx + dy * dy + dz * dz + eps2);
            potencial -= G_GRAVEDAD * s->masa[i] * s->masa[j] / r;
        }
    }

    return cinetica + potencial;
}
