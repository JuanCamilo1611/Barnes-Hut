#!/usr/bin/env python3
"""
Compara la suma directa contra el Barnes-Hut de arbol plano.

    python3 tools/benchmark.py --N 5000
    python3 tools/benchmark.py --N 20000 --sin-directa

La suma directa es O(N^2): con N grande tarda demasiado, por eso se
puede omitir con --sin-directa y dejar solo el tiempo del arbol plano.
"""

import argparse
import os
import random
import sys
import time

import numpy as np

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, RAIZ)
sys.path.insert(0, os.path.join(RAIZ, "python_mp"))

from arbol import ArbolPlano


def generar(n, seed=42):
    rng = random.Random(seed)
    posiciones = np.array(
        [[rng.uniform(0, 1), rng.uniform(0, 1), rng.uniform(0, 1)]
         for _ in range(n)],
        dtype=np.float64,
    )
    masas = np.full(n, 1.0 / n, dtype=np.float64)
    return posiciones, masas


def directa(posiciones, masas, epsilon=0.01, G=1.0):
    """Suma directa vectorizada por bloques (solo para referencia)."""
    n = len(posiciones)
    epsilon_cuadrado = epsilon * epsilon
    aceleraciones = np.zeros((n, 3), dtype=np.float64)
    bloque = 256

    for inicio in range(0, n, bloque):
        fin = min(inicio + bloque, n)

        diferencia = (
            posiciones[None, :, :] - posiciones[inicio:fin, None, :]
        )
        distancia_cuadrada = (
            np.einsum("ijk,ijk->ij", diferencia, diferencia)
            + epsilon_cuadrado
        )

        factor = (
            G * masas[None, :]
            / distancia_cuadrada ** 1.5
        )
        # la propia particula no aporta: la diferencia es cero y el
        # factor sigue finito, asi que se anula la diagonal
        filas = np.arange(fin - inicio)
        factor[filas, np.arange(inicio, fin)] = 0.0

        aceleraciones[inicio:fin] = np.einsum(
            "ij,ijk->ik", factor, diferencia
        )

    return aceleraciones


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, default=5000)
    parser.add_argument("--theta", type=float, default=0.3)
    parser.add_argument("--eps", type=float, default=0.01)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--sin-directa", action="store_true")
    parser.add_argument("--repeticiones", type=int, default=1)
    args = parser.parse_args()

    posiciones, masas = generar(args.N, args.seed)

    print("BENCHMARK  N=%d  theta=%.2f  eps=%.3f"
          % (args.N, args.theta, args.eps))

    if not args.sin_directa:
        inicio = time.perf_counter()
        directas = directa(posiciones, masas, args.eps)
        tiempo_directo = time.perf_counter() - inicio
        print("  suma directa      : %8.3f s" % tiempo_directo)
    else:
        directas = None
        tiempo_directo = None

    mejor_build = mejor_walk = float("inf")

    for _ in range(args.repeticiones):
        arbol = ArbolPlano()

        inicio = time.perf_counter()
        arbol.construir(posiciones, masas)
        build = time.perf_counter() - inicio

        inicio = time.perf_counter()
        aceleraciones = arbol.aceleraciones(
            posiciones, theta=args.theta, epsilon=args.eps
        )
        walk = time.perf_counter() - inicio

        mejor_build = min(mejor_build, build)
        mejor_walk = min(mejor_walk, walk)

    print("  construccion      : %8.3f s" % mejor_build)
    print("  recorrido         : %8.3f s" % mejor_walk)
    print("  arbol plano total : %8.3f s" % (mejor_build + mejor_walk))
    print("  nodos             : %d" % arbol.n_nodos)

    if directas is not None:
        maximo, promedio = 0.0, 0.0
        usadas = 0
        for i in range(args.N):
            ref = directas[i]
            tes = aceleraciones[i]
            diferencia = np.linalg.norm(tes - ref)
            magnitud = np.linalg.norm(ref)
            if magnitud < 1e-300:
                continue
            error = 100.0 * diferencia / magnitud
            maximo = max(maximo, error)
            promedio += error
            usadas += 1
        promedio /= max(usadas, 1)

        print("  error vs directa  : max %.6f %%  promedio %.6f %%"
              % (maximo, promedio))

        total_bh = mejor_build + mejor_walk
        print("  relacion BH/directa: %.3fx" % (total_bh / tiempo_directo))

        if total_bh < tiempo_directo:
            print("  OK: el arbol plano es mas rapido que la suma directa")
        else:
            print("  FALLO: el arbol plano aun no gana a la suma directa")
            return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
