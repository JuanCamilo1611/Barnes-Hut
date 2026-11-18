#include "leapfrog.h"
#include "tiempo.h"

/* v = v + (dt/2) * a */
static void media_patada(Sistema *s, double dt)
{
    const double h = 0.5 * dt;
    for (int k = 0; k < 3 * s->n; k++)
        s->vel[k] += h * s->acc[k];
}

/* x = x + dt * v */
static void desplazar(Sistema *s, double dt)
{
    for (int k = 0; k < 3 * s->n; k++)
        s->pos[k] += dt * s->vel[k];
}

void simular_leapfrog(Sistema *s, int pasos, double dt,
                      FuncionFuerzas calcular, void *config,
                      double *tiempo_fuerzas, double *tiempo_integrar)
{
    double t_fuerzas = 0.0, t_integrar = 0.0;

    for (int paso = 0; paso < pasos; paso++) {
        double t0 = ahora();
        media_patada(s, dt);      /* 1. medio paso de velocidad con a(t)     */
        desplazar(s, dt);         /* 2. paso completo de posicion            */
        double t1 = ahora();

        calcular(s, config);      /* 3. a(t+dt): UNA sola vez por paso       */
        double t2 = ahora();

        media_patada(s, dt);      /* 4. otro medio paso de velocidad con a(t+dt) */
        double t3 = ahora();      /*    s->acc queda lista para el paso siguiente */

        t_integrar += (t1 - t0) + (t3 - t2);
        t_fuerzas  += (t2 - t1);
    }

    if (tiempo_fuerzas)  *tiempo_fuerzas  += t_fuerzas;
    if (tiempo_integrar) *tiempo_integrar += t_integrar;
}
