#!/usr/bin/env python3
"""
Genera las distribuciones de particulas del proyecto.

Formato de salida (igual para todas las distribuciones):
    una linea por particula con siete columnas:
    x y z vx vy vz masa
    masa = 1/N en todos los casos, semilla fija para reproducibilidad.

Uso:
    python3 data/generate_particles.py --distribucion uniforme --N 1000
    python3 data/generate_particles.py --distribucion cumulos --N 20000 --out data/cumulos.txt
    python3 data/generate_particles.py --distribucion concentradas --N 200000
"""

import argparse
import os
import random

DISTRIBUCIONES = ("uniforme", "cumulos", "concentradas")


def generar_uniforme(n, rng):
    """Posiciones uniformes en el cubo unitario, velocidades en cero."""
    particulas = []
    for _ in range(n):
        particulas.append([
            rng.uniform(0.0, 1.0),
            rng.uniform(0.0, 1.0),
            rng.uniform(0.0, 1.0),
            0.0, 0.0, 0.0,
        ])
    return particulas


def generar_cumulos(n, rng):
    """
    Dos cumulos gaussianos que se mueven el uno hacia el otro.

    Castiga un reparto ingenuo del trabajo: la mitad de las
    particulas esta en un lado y la otra mitad en el otro, con
    velocidades opuestas, de modo que colisionan en el centro.
    """
    sigma = 0.12
    velocidad = 0.5
    separacion = 0.45

    mitad = n // 2
    conteo = [mitad, n - mitad]
    centros = [
        (0.5 - separacion, 0.5, 0.5),
        (0.5 + separacion, 0.5, 0.5),
    ]
    velocidades = [
        (+velocidad, 0.0, 0.0),
        (-velocidad, 0.0, 0.0),
    ]

    particulas = []

    for cumulo in range(2):
        cx, cy, cz = centros[cumulo]
        vx, vy, vz = velocidades[cumulo]

        for _ in range(conteo[cumulo]):
            particulas.append([
                rng.gauss(cx, sigma),
                rng.gauss(cy, sigma),
                rng.gauss(cz, sigma),
                vx, vy, vz,
            ])

    return particulas


def generar_concentradas(n, rng):
    """
    Casi todas las particulas dentro de un volumen muy pequeno.

    Castiga la profundidad del arbol y los criterios de apertura:
    muchas particulas muy cerca unas de otras.
    """
    sigma = 0.01
    centro = 0.5

    particulas = []

    for _ in range(n):
        particulas.append([
            rng.gauss(centro, sigma),
            rng.gauss(centro, sigma),
            rng.gauss(centro, sigma),
            0.0, 0.0, 0.0,
        ])

    return particulas


GENERADORES = {
    "uniforme": generar_uniforme,
    "cumulos": generar_cumulos,
    "concentradas": generar_concentradas,
}


def guardar(particulas, ruta):
    """Escribe las particulas en el formato de siete columnas."""

    directorio = os.path.dirname(ruta)
    if directorio:
        os.makedirs(directorio, exist_ok=True)

    with open(ruta, "w") as archivo:
        for particula in particulas:
            masa = 1.0 / len(particulas)
            valores = particula + [masa]
            archivo.write(
                " ".join("%.8f" % valor for valor in valores) + "\n"
            )


def main():
    parser = argparse.ArgumentParser(
        description="Genera distribuciones de particulas."
    )
    parser.add_argument(
        "--distribucion",
        choices=DISTRIBUCIONES,
        default="uniforme",
    )
    parser.add_argument("--N", type=int, default=100)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out", default="data/particles.txt")

    args = parser.parse_args()

    if args.N < 1:
        parser.error("--N debe ser al menos 1")

    rng = random.Random(args.seed)
    particulas = GENERADORES[args.distribucion](args.N, rng)

    guardar(particulas, args.out)

    print("Distribucion:", args.distribucion)
    print("Particulas:", args.N)
    print("Semilla:", args.seed)
    print("Masa por particula:", 1.0 / args.N)
    print("Archivo:", args.out)


if __name__ == "__main__":
    main()
