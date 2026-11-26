/*
 * sim_cuda.cu
 * Simulacion Barnes-Hut con CUDA (Colab, GPU T4).
 *
 * Uso (contrato comun del README):
 *   sim_cuda <archivo> <N> <pasos> <dt> <theta> <eps> <hilos> [salida_acc]
 * 'hilos' = hilos por bloque de CUDA (si es menor que 32 se usan 128;
 *           asi verificar.py, que pasa hilos = 1, funciona igual).
 *
 * Fases de cada calculo de fuerzas (tiempos con cudaEvent):
 *   sort  : caja envolvente, codigos Morton en GPU (1 hilo por particula)
 *           y ordenamiento en GPU con thrust::sort_by_key (radix sort,
 *           estable: empates por indice, igual que el nucleo C);
 *   build : arbol plano en preorden con el nucleo C (construir_arbol) a
 *           partir del orden Morton calculado en la GPU, y copia a la GPU
 *           de los 6 arreglos del contrato (SoA -> accesos coalescentes);
 *   walk  : recorrido SIN recursion (pila explicita en cada hilo). Un hilo
 *           por particula EN ORDEN MORTON: los 32 hilos de un warp son
 *           particulas vecinas, abren casi los mismos nodos -> poca
 *           divergencia;
 *   integrate: leapfrog en GPU (dos medias patadas y un desplazamiento).
 *
 * Todo en double (float64), igual que la referencia.
 */
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cuda_runtime.h>
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#include <thrust/sequence.h>

extern "C" {
#include "particulas.h"
#include "directa.h"
#include "morton.h"
#include "arbol.h"
}

#define N_MAXIMO_ENERGIA 5000
#define HILOS_POR_DEFECTO 128
#define PILA_GPU (8 * (MAX_PROFUNDIDAD + 2))   /* igual que recorrido.h */

/* Detiene el programa con un mensaje claro si una llamada de CUDA falla. */
#define CUDA_OK(llamada)                                                   \
    do {                                                                   \
        cudaError_t err_ = (llamada);                                      \
        if (err_ != cudaSuccess) {                                         \
            fprintf(stderr, "CUDA error %s en %s:%d\n",                    \
                    cudaGetErrorString(err_), __FILE__, __LINE__);         \
            exit(1);                                                       \
        }                                                                  \
    } while (0)

/* ------------------------------------------------------------------ */
/* Kernels                                                             */
/* ------------------------------------------------------------------ */

__device__ uint64_t separar_bits_gpu(uint32_t v)
{
    uint64_t x = v & 0x1fffff;
    x = (x | (x << 32)) & 0x1f00000000ffffULL;
    x = (x | (x << 16)) & 0x1f0000ff0000ffULL;
    x = (x | (x << 8))  & 0x100f00f00f00f00fULL;
    x = (x | (x << 4))  & 0x10c30c30c30c30c3ULL;
    x = (x | (x << 2))  & 0x1249249249249249ULL;
    return x;
}

__device__ uint32_t a_entero_gpu(double p, double minimo, double lado)
{
    const double escala = (double)(1u << MORTON_BITS);
    double t = (p - minimo) / lado * escala;
    if (t < 0.0) t = 0.0;
    if (t > escala - 1.0) t = escala - 1.0;
    return (uint32_t)t;
}

/* Un hilo por particula: codigo Morton (mismas operaciones que morton.c). */
__global__ void kernel_morton(int n, const double *pos, double mx, double my,
                              double mz, double lado, uint64_t *codigos)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    uint32_t ix = a_entero_gpu(pos[3 * i + 0], mx, lado);
    uint32_t iy = a_entero_gpu(pos[3 * i + 1], my, lado);
    uint32_t iz = a_entero_gpu(pos[3 * i + 2], mz, lado);
    codigos[i] = (separar_bits_gpu(ix) << 2) | (separar_bits_gpu(iy) << 1)
               | separar_bits_gpu(iz);
}

/* Recorrido sin recursion. El hilo j calcula la particula orden[j]. */
__global__ void kernel_recorrido(int n, const int *orden, const double *pos,
                                 double *acc, const double *centro,
                                 const double *lado, const double *masa,
                                 const double *com, const int32_t *hijos,
                                 const int32_t *part, double theta2,
                                 double eps2)
{
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j >= n) return;
    int i = orden[j];

    double px = pos[3 * i], py = pos[3 * i + 1], pz = pos[3 * i + 2];
    double ax = 0.0, ay = 0.0, az = 0.0;

    int pila[PILA_GPU];
    int tope = 0;
    pila[tope++] = 0;                       /* raiz */

    while (tope > 0) {
        int nodo = pila[--tope];
        double m = masa[nodo];
        if (m == 0.0) continue;             /* celda vacia */

        double dx = com[3 * nodo + 0] - px;
        double dy = com[3 * nodo + 1] - py;
        double dz = com[3 * nodo + 2] - pz;
        double d2 = dx * dx + dy * dy + dz * dz;
        int p = part[nodo];
        int es_hoja = (p >= 0);
        if (es_hoja && p == i) continue;    /* no atraerse a si misma */

        int abrir = 0;
        if (!es_hoja) {
            double s = lado[nodo];
            abrir = (s * s >= theta2 * d2);
            if (!abrir) {                   /* particula dentro de la celda */
                double h = 0.5 * s;
                abrir = fabs(px - centro[3 * nodo + 0]) <= h &&
                        fabs(py - centro[3 * nodo + 1]) <= h &&
                        fabs(pz - centro[3 * nodo + 2]) <= h;
            }
        }

        if (abrir) {
            for (int o = 7; o >= 0; o--) {  /* se visitan 0..7 */
                int h = hijos[8 * nodo + o];
                if (h >= 0) pila[tope++] = h;
            }
        } else {
            double r2 = d2 + eps2;
            double f = G_GRAVEDAD * m / (r2 * sqrt(r2));
            ax += f * dx;
            ay += f * dy;
            az += f * dz;
        }
    }
    acc[3 * i + 0] = ax;
    acc[3 * i + 1] = ay;
    acc[3 * i + 2] = az;
}

/* v += (dt/2) a      (una componente por hilo: 3n hilos) */
__global__ void kernel_media_patada(int n3, double *vel, const double *acc,
                                    double h)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < n3) vel[k] += h * acc[k];
}

/* x += dt v */
__global__ void kernel_desplazar(int n3, double *pos, const double *vel,
                                 double dt)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < n3) pos[k] += dt * vel[k];
}

/* ------------------------------------------------------------------ */
/* Datos en la GPU y tiempos                                           */
/* ------------------------------------------------------------------ */

typedef struct {
    int n, hilos;
    double theta, eps;
    /* particulas */
    double *pos, *vel, *acc;
    uint64_t *codigos;
    int *orden;
    /* arbol (los 6 arreglos del contrato) */
    int cap;
    double *centro, *lado, *masa, *com;
    int32_t *hijos, *part;
    /* tiempos acumulados en ms (cudaEvent) */
    float t_sort, t_build, t_walk, t_integrate;
    cudaEvent_t e0, e1, e2, e3;
} GPU;

static int bloques(int total, int hilos) { return (total + hilos - 1) / hilos; }

static void reservar_arbol_gpu(GPU *g, int nodos)
{
    if (nodos <= g->cap) return;
    if (g->cap > 0) {
        cudaFree(g->centro); cudaFree(g->lado); cudaFree(g->masa);
        cudaFree(g->com);    cudaFree(g->hijos); cudaFree(g->part);
    }
    g->cap = nodos + nodos / 4;   /* margen para pasos siguientes */
    size_t c = (size_t)g->cap;
    CUDA_OK(cudaMalloc(&g->centro, 3 * c * sizeof(double)));
    CUDA_OK(cudaMalloc(&g->lado,       c * sizeof(double)));
    CUDA_OK(cudaMalloc(&g->masa,       c * sizeof(double)));
    CUDA_OK(cudaMalloc(&g->com,    3 * c * sizeof(double)));
    CUDA_OK(cudaMalloc(&g->hijos,  8 * c * sizeof(int32_t)));
    CUDA_OK(cudaMalloc(&g->part,       c * sizeof(int32_t)));
}

/* Calcula g->acc (en la GPU). 's' es la copia en CPU (para el arbol). */
static void fuerzas_gpu(GPU *g, Sistema *s, ArbolPlano *arbol,
                        uint64_t *codigos, int *orden)
{
    int n = g->n, H = g->hilos;
    float ms;

    /* --- sort: caja, Morton y ordenamiento en GPU --- */
    CUDA_OK(cudaEventRecord(g->e0));
    CUDA_OK(cudaMemcpy(s->pos, g->pos, 3 * (size_t)n * sizeof(double),
                       cudaMemcpyDeviceToHost));
    double minimo[3], lado;
    caja_envolvente(s, minimo, &lado);
    kernel_morton<<<bloques(n, H), H>>>(n, g->pos, minimo[0], minimo[1],
                                        minimo[2], lado, g->codigos);
    CUDA_OK(cudaGetLastError());
    thrust::device_ptr<uint64_t> k(g->codigos);
    thrust::device_ptr<int> v(g->orden);
    thrust::sequence(v, v + n);
    /* Se ordena una COPIA de las claves: el arbol necesita los codigos
     * indexados por particula (codigos[i]), no ordenados. */
    CUDA_OK(cudaMemcpy(codigos, g->codigos, (size_t)n * sizeof(uint64_t),
                       cudaMemcpyDeviceToHost));
    thrust::sort_by_key(k, k + n, v);
    CUDA_OK(cudaMemcpy(orden, g->orden, (size_t)n * sizeof(int),
                       cudaMemcpyDeviceToHost));
    CUDA_OK(cudaEventRecord(g->e1));

    /* --- build: arbol plano (preorden) y copia a la GPU --- */
    construir_arbol(arbol, s, orden, codigos, minimo, lado);
    int m = arbol->n_nodos;
    reservar_arbol_gpu(g, m);
    CUDA_OK(cudaMemcpy(g->centro, arbol->centro, 3 * (size_t)m * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_OK(cudaMemcpy(g->lado,   arbol->lado,       (size_t)m * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_OK(cudaMemcpy(g->masa,   arbol->masa,       (size_t)m * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_OK(cudaMemcpy(g->com,    arbol->com,    3 * (size_t)m * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_OK(cudaMemcpy(g->hijos,  arbol->hijos,  8 * (size_t)m * sizeof(int32_t), cudaMemcpyHostToDevice));
    CUDA_OK(cudaMemcpy(g->part,   arbol->part,       (size_t)m * sizeof(int32_t), cudaMemcpyHostToDevice));
    CUDA_OK(cudaEventRecord(g->e2));

    /* --- walk: un hilo por particula en orden Morton --- */
    kernel_recorrido<<<bloques(n, H), H>>>(n, g->orden, g->pos, g->acc,
        g->centro, g->lado, g->masa, g->com, g->hijos, g->part,
        g->theta * g->theta, g->eps * g->eps);
    CUDA_OK(cudaGetLastError());
    CUDA_OK(cudaEventRecord(g->e3));
    CUDA_OK(cudaEventSynchronize(g->e3));

    CUDA_OK(cudaEventElapsedTime(&ms, g->e0, g->e1)); g->t_sort  += ms;
    CUDA_OK(cudaEventElapsedTime(&ms, g->e1, g->e2)); g->t_build += ms;
    CUDA_OK(cudaEventElapsedTime(&ms, g->e2, g->e3)); g->t_walk  += ms;
}

int main(int argc, char **argv)
{
    if (argc < 8) {
        fprintf(stderr,
                "uso: %s <archivo> <N> <pasos> <dt> <theta> <eps> <hilos> [salida_acc]\n",
                argv[0]);
        return 1;
    }
    const char *archivo = argv[1];
    int n        = atoi(argv[2]);
    int pasos    = atoi(argv[3]);
    double dt    = atof(argv[4]);
    double theta = atof(argv[5]);
    double eps   = atof(argv[6]);
    int hilos    = atoi(argv[7]);
    const char *salida_acc = (argc > 8) ? argv[8] : "aceleraciones.txt";

    Sistema s;
    if (leer_particulas(archivo, n, &s) != 0)
        return 1;

    GPU g = {0};
    g.n = n;
    g.hilos = (hilos >= 32) ? hilos : HILOS_POR_DEFECTO;
    g.theta = theta;
    g.eps = eps;
    size_t b3 = 3 * (size_t)n * sizeof(double);
    CUDA_OK(cudaMalloc(&g.pos, b3));
    CUDA_OK(cudaMalloc(&g.vel, b3));
    CUDA_OK(cudaMalloc(&g.acc, b3));
    CUDA_OK(cudaMalloc(&g.codigos, (size_t)n * sizeof(uint64_t)));
    CUDA_OK(cudaMalloc(&g.orden, (size_t)n * sizeof(int)));
    CUDA_OK(cudaMemcpy(g.pos, s.pos, b3, cudaMemcpyHostToDevice));
    CUDA_OK(cudaMemcpy(g.vel, s.vel, b3, cudaMemcpyHostToDevice));
    CUDA_OK(cudaEventCreate(&g.e0)); CUDA_OK(cudaEventCreate(&g.e1));
    CUDA_OK(cudaEventCreate(&g.e2)); CUDA_OK(cudaEventCreate(&g.e3));

    ArbolPlano arbol = crear_arbol(2 * n);
    uint64_t *codigos = (uint64_t *)malloc((size_t)n * sizeof(uint64_t));
    int *orden = (int *)malloc((size_t)n * sizeof(int));

    /* Aceleraciones iniciales: para verificar.py y para el primer paso. */
    fuerzas_gpu(&g, &s, &arbol, codigos, orden);
    CUDA_OK(cudaMemcpy(s.acc, g.acc, b3, cudaMemcpyDeviceToHost));
    if (escribir_aceleraciones(salida_acc, &s) != 0)
        return 1;

    int con_energia = (n <= N_MAXIMO_ENERGIA);
    double e_inicial = con_energia ? energia_total(&s, eps) : 0.0;

    /* Leapfrog con reuso de fuerzas, todo en la GPU. */
    int n3 = 3 * n, B = bloques(n3, g.hilos);
    cudaEvent_t i0, i1;
    CUDA_OK(cudaEventCreate(&i0)); CUDA_OK(cudaEventCreate(&i1));
    float ms;
    for (int paso = 0; paso < pasos; paso++) {
        CUDA_OK(cudaEventRecord(i0));
        kernel_media_patada<<<B, g.hilos>>>(n3, g.vel, g.acc, 0.5 * dt);
        kernel_desplazar<<<B, g.hilos>>>(n3, g.pos, g.vel, dt);
        CUDA_OK(cudaEventRecord(i1));
        CUDA_OK(cudaEventSynchronize(i1));
        CUDA_OK(cudaEventElapsedTime(&ms, i0, i1)); g.t_integrate += ms;

        fuerzas_gpu(&g, &s, &arbol, codigos, orden);

        CUDA_OK(cudaEventRecord(i0));
        kernel_media_patada<<<B, g.hilos>>>(n3, g.vel, g.acc, 0.5 * dt);
        CUDA_OK(cudaEventRecord(i1));
        CUDA_OK(cudaEventSynchronize(i1));
        CUDA_OK(cudaEventElapsedTime(&ms, i0, i1)); g.t_integrate += ms;
    }
    CUDA_OK(cudaGetLastError());

    /* Estado final a la CPU (para la energia). */
    CUDA_OK(cudaMemcpy(s.pos, g.pos, b3, cudaMemcpyDeviceToHost));
    CUDA_OK(cudaMemcpy(s.vel, g.vel, b3, cudaMemcpyDeviceToHost));

    double t_total = (g.t_sort + g.t_build + g.t_walk + g.t_integrate) / 1000.0;
    printf("hilos_por_bloque=%d\n", g.hilos);
    printf("tiempo_build=%.6f\n", g.t_build / 1000.0);
    printf("tiempo_sort=%.6f\n", g.t_sort / 1000.0);
    printf("tiempo_walk=%.6f\n", g.t_walk / 1000.0);
    printf("tiempo_integrate=%.6f\n", g.t_integrate / 1000.0);
    printf("tiempo_total=%.6f\n", t_total);
    if (con_energia) {
        printf("energia_inicial=%.12f\n", e_inicial);
        printf("energia_final=%.12f\n", energia_total(&s, eps));
    }

    cudaFree(g.pos); cudaFree(g.vel); cudaFree(g.acc);
    cudaFree(g.codigos); cudaFree(g.orden);
    if (g.cap > 0) {
        cudaFree(g.centro); cudaFree(g.lado); cudaFree(g.masa);
        cudaFree(g.com); cudaFree(g.hijos); cudaFree(g.part);
    }
    liberar_arbol(&arbol);
    free(codigos); free(orden);
    liberar_sistema(&s);
    return 0;
}