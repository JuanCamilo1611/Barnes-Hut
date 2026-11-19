#include <math.h>
#include "recorrido.h"
#include "directa.h"   /* G_GRAVEDAD */

void aceleracion_particula(const ArbolPlano *a, Sistema *s, int i,
                           double theta, double eps)
{
    const double *p = &s->pos[3 * i];
    const double theta2 = theta * theta, eps2 = eps * eps;
    double ax = 0.0, ay = 0.0, az = 0.0;

    int pila[TAM_PILA];
    int tope = 0;
    pila[tope++] = 0;                       /* empezar en la raiz */

    while (tope > 0) {
        int n = pila[--tope];
        double m = a->masa[n];
        if (m == 0.0)
            continue;                       /* celda vacia: no aporta */

        double dx = a->com[3 * n + 0] - p[0];
        double dy = a->com[3 * n + 1] - p[1];
        double dz = a->com[3 * n + 2] - p[2];
        double d2 = dx * dx + dy * dy + dz * dz;
        int es_hoja = (a->part[n] >= 0);

        if (es_hoja && a->part[n] == i)
            continue;                       /* no atraerse a si misma */

        /* Abrir si s/d >= theta, es decir s^2 >= theta^2 * d^2.
         * Tambien se abre si la particula esta dentro de la celda
         * (con theta grande el com podria quedar "lejos" igual). */
        int abrir = 0;
        if (!es_hoja) {
            double s_ = a->lado[n];
            abrir = (s_ * s_ >= theta2 * d2);
            if (!abrir) {
                double h = 0.5 * s_;
                abrir = fabs(p[0] - a->centro[3 * n + 0]) <= h &&
                        fabs(p[1] - a->centro[3 * n + 1]) <= h &&
                        fabs(p[2] - a->centro[3 * n + 2]) <= h;
            }
        }

        if (abrir) {
            for (int o = 7; o >= 0; o--) {  /* al reves: se visitan 0..7 */
                int h = a->hijos[8 * n + o];
                if (h >= 0)
                    pila[tope++] = h;
            }
        } else {
            /* La celda completa actua como una sola masa en su com. */
            double r2 = d2 + eps2;
            double f = G_GRAVEDAD * m / (r2 * sqrt(r2));
            ax += f * dx;
            ay += f * dy;
            az += f * dz;
        }
    }

    s->acc[3 * i + 0] = ax;
    s->acc[3 * i + 1] = ay;
    s->acc[3 * i + 2] = az;
}

void aceleraciones_arbol(const ArbolPlano *a, Sistema *s,
                         double theta, double eps)
{
    for (int i = 0; i < s->n; i++)
        aceleracion_particula(a, s, i, theta, eps);
}
