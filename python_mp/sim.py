#!/usr/bin/env python3
"""
Version Python de Barnes-Hut con arbol plano en arreglos.

Contrato comun (ver README.md):

    python3 python_mp/sim.py <archivo> <N> <pasos> <dt> <theta> <eps> \
        <hilos> [salida_acc]

Imprime en stdout las lineas clave=valor con los tiempos por fase y la
energia, y escribe las aceleraciones iniciales en el archivo pedido.
"""

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from arbol import ArbolPlano

# Tamano del lote de particulas del recorrido vectorizado.  Se puede
# cambiar con BH_LOTE para medir el efecto sobre el uso de cache.
LOTE = int(os.environ.get("BH_LOTE", "2048"))


def leer_particulas(ruta):
    """Lee el archivo de N filas y 7 columnas."""

    datos = np.loadtxt(ruta, dtype=np.float64)

    if datos.ndim == 1:
        datos = datos.reshape(1, -1)

    if datos.shape[1] != 7:
        raise SystemExit(
            "se esperaban 7 columnas (x y z vx vy vz masa) en %s" % ruta
        )

    return (
        np.ascontiguousarray(datos[:, 0:3]),
        np.ascontiguousarray(datos[:, 3:6]),
        np.ascontiguousarray(datos[:, 6]),
    )


def energia_cinetica(velocidades, masas):
    rapidez_cuadrada = np.einsum("ij,ij->i", velocidades, velocidades)
    return 0.5 * float(np.sum(masas * rapidez_cuadrada))


def energia_potencial(posiciones, masas, epsilon, G=1.0, bloque=512):
    """
    Energia potencial gravitacional total.

    Se procesa en bloques para no crear una matriz de NxN completa y
    cada pareja se cuenta una sola vez.
    """
    n = len(posiciones)
    epsilon_cuadrado = epsilon * epsilon
    indices = np.arange(n)
    total = 0.0

    for inicio in range(0, n, bloque):
        fin = min(inicio + bloque, n)

        diferencia = (
            posiciones[inicio:fin, None, :] - posiciones[None, :, :]
        )
        distancia_cuadrada = (
            np.einsum("ijk,ijk->ij", diferencia, diferencia)
            + epsilon_cuadrado
        )

        # Solo las parejas (i, j) con j > i.
        mascara = indices[None, :] > indices[inicio:fin, None]

        productos = (
            masas[inicio:fin, None] * masas[None, :]
        ) / np.sqrt(distancia_cuadrada)

        total += float(np.sum(productos[mascara]))

    return -G * total


def calcular_energia(posiciones, velocidades, masas, epsilon, G=1.0):
    cinetica = energia_cinetica(velocidades, masas)
    potencial = energia_potencial(posiciones, masas, epsilon, G)
    return cinetica + potencial


def calcular_fuerzas(arbol, paralelo, posiciones, masas, theta, G,
                     epsilon):
    """
    Construye el arbol y calcula las aceleraciones.

    Devuelve (aceleraciones, segundos_construir, segundos_recorrer).
    """

    if paralelo is None:
        momento = time.perf_counter()
        arbol.construir(posiciones, masas)
        tiempo_construir = time.perf_counter() - momento

        momento = time.perf_counter()
        aceleraciones = arbol.aceleraciones(
            posiciones, theta=theta, G=G, epsilon=epsilon, lote=LOTE
        )
        tiempo_recorrer = time.perf_counter() - momento

        return aceleraciones, tiempo_construir, tiempo_recorrer

    tiempo_construir = paralelo.construir(posiciones)
    aceleraciones, tiempos, tiempo_pared = paralelo.recorrer(posiciones)

    if paralelo.estrategia == "compartido":
        # El arbol ya estaba construido: todo el tiempo de la
        # llamada es recorrido.
        return aceleraciones, tiempo_construir, tiempo_pared

    # En la estrategia local cada trabajador construye el suyo.
    tiempo_construir += max(t[0] for t in tiempos)
    tiempo_recorrer = max(t[1] for t in tiempos)

    return aceleraciones, tiempo_construir, tiempo_recorrer


def main():
    if len(sys.argv) < 8:
        raise SystemExit(
            "uso: sim.py <archivo> <N> <pasos> <dt> <theta> <eps> "
            "<hilos> [salida_acc]"
        )

    ruta = sys.argv[1]
    n_pedidas = int(sys.argv[2])
    pasos = int(sys.argv[3])
    dt = float(sys.argv[4])
    theta = float(sys.argv[5])
    epsilon = float(sys.argv[6])
    hilos = int(sys.argv[7])
    salida_acc = sys.argv[8] if len(sys.argv) > 8 else "aceleraciones.txt"

    G = 1.0

    posiciones, velocidades, masas = leer_particulas(ruta)

    if len(posiciones) != n_pedidas:
        raise SystemExit(
            "N=%d pero el archivo tiene %d particulas"
            % (n_pedidas, len(posiciones))
        )

    arbol = None
    paralelo = None

    if hilos > 1:
        from paralelo import Paralelo

        paralelo = Paralelo(
            posiciones, masas, hilos,
            theta=theta, G=G, epsilon=epsilon, lote=LOTE,
        )
    else:
        arbol = ArbolPlano()

    # ---------------------------------------------------------
    # se arranca dentro del try para que cerrar() corra siempre
    # ---------------------------------------------------------
    momento = time.perf_counter()
    tiempo_integrate = 0.0

    try:
        # -----------------------------------------------------
        # aceleraciones iniciales: se escriben en el archivo pedido
        # -----------------------------------------------------
        aceleraciones, tiempo_build, tiempo_walk = calcular_fuerzas(
            arbol, paralelo, posiciones, masas, theta, G, epsilon
        )

        np.savetxt(salida_acc, aceleraciones, fmt="%.17g")

        if len(posiciones) <= 5000:
            energia_inicial = calcular_energia(
                posiciones, velocidades, masas, epsilon, G
            )
        else:
            energia_inicial = None

        # -----------------------------------------------------
        # leapfrog con una sola evaluacion de fuerzas por paso:
        # las aceleraciones del final de un paso son las del
        # inicio del siguiente
        # -----------------------------------------------------
        for _ in range(pasos):
            inicio = time.perf_counter()
            velocidades += 0.5 * aceleraciones * dt
            posiciones += velocidades * dt
            tiempo_integrate += time.perf_counter() - inicio

            aceleraciones, build, walk = calcular_fuerzas(
                arbol, paralelo, posiciones, masas, theta, G, epsilon
            )
            tiempo_build += build
            tiempo_walk += walk

            inicio = time.perf_counter()
            velocidades += 0.5 * aceleraciones * dt
            tiempo_integrate += time.perf_counter() - inicio

        tiempo_total = time.perf_counter() - momento

    finally:
        if paralelo is not None:
            paralelo.cerrar()

    if energia_inicial is not None:
        energia_final = calcular_energia(
            posiciones, velocidades, masas, epsilon, G
        )
        print("energia_inicial=%.17g" % energia_inicial)
        print("energia_final=%.17g" % energia_final)

    print("hilos=%d" % hilos)
    print("tiempo_build=%.6f" % tiempo_build)
    print("tiempo_sort=0.000000")
    print("tiempo_walk=%.6f" % tiempo_walk)
    print("tiempo_integrate=%.6f" % tiempo_integrate)
    print("tiempo_total=%.6f" % tiempo_total)


if __name__ == "__main__":
    main()
