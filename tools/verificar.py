#!/usr/bin/env python3
"""
Verificacion de correctitud de todas las versiones del proyecto.

Comando de puerta antes de cada PR:

    python3 tools/verificar.py --N 2000

Verificar una version externa (C, MPI, CUDA, multiprocessing):

    python3 tools/verificar.py \
        --particulas data/ejemplo_100.txt \
        --programa "python3 python_mp/sim.py" \
        --hilos 4

Comparar dos archivos de aceleraciones ya generados:

    python3 tools/verificar.py --comparar acc_a.txt acc_b.txt

Contrato exigido a los programas externos:

    programa <archivo> <N> <pasos> <dt> <theta> <eps> <hilos> [salida_acc]

  * escribe las aceleraciones iniciales en <salida_acc> (por defecto
    aceleraciones.txt): una linea por particula, en el mismo orden que
    el archivo de entrada, con tres valores "ax ay az";
  * imprime en stdout lineas clave=valor, por ejemplo:

        tiempo_build=0.0123
        tiempo_sort=0.0004
        tiempo_walk=0.4500
        tiempo_integrate=0.0100
        tiempo_total=0.4727
        energia_inicial=-1.234567
        energia_final=-1.234500

  * termina con codigo de salida 0.

Criterios de aprobacion:

  * error relativo promedio contra la suma directa <= --tol-error (1.0 %)
    con N <= 2000, donde la suma directa es viable;
  * derivacion de la energia mecanica <= --tol-energia (1.0 %).
"""

import argparse
import math
import os
import random
import shlex
import subprocess
import sys
import tempfile

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, RAIZ)

from referencia.motor.fuerzas import (
    cargar_particulas,
    calcular_energia_total,
    calcular_todas_las_aceleraciones,
)
from referencia.motor.fuerzas_barnes_hut import (
    calcular_todas_las_aceleraciones_barnes_hut,
)

# La suma directa es O(N^2): mas alla de este N no se evalua.
N_MAXIMO_DIRECTA = 2000
# La energia mecanica tambien es O(N^2).
N_MAXIMO_ENERGIA = 5000


def generar_particulas(n, seed):
    """Distribucion uniforme reproducible, igual que data/generate_particles.py."""
    rng = random.Random(seed)
    masa = 1.0 / n
    particulas = []
    for _ in range(n):
        particulas.append({
            "posicion": [rng.uniform(0.0, 1.0) for _ in range(3)],
            "velocidad": [0.0, 0.0, 0.0],
            "masa": masa,
        })
    return particulas


def escribir_particulas(particulas, ruta):
    with open(ruta, "w") as archivo:
        for p in particulas:
            valores = p["posicion"] + p["velocidad"] + [p["masa"]]
            archivo.write(
                " ".join("%.8f" % v for v in valores) + "\n"
            )


def leer_aceleraciones(ruta):
    """Lee un archivo de aceleraciones: tres valores por linea."""
    aceleraciones = []
    with open(ruta) as archivo:
        for numero, linea in enumerate(archivo, 1):
            linea = linea.strip()
            if not linea or linea.startswith("#"):
                continue
            partes = linea.split()
            # Cabecera opcional con el numero de particulas.
            if numero == 1 and len(partes) == 1:
                continue
            if len(partes) < 3:
                raise ValueError(
                    "%s linea %d: se esperaban 3 valores" % (ruta, numero)
                )
            aceleraciones.append(
                [float(partes[0]), float(partes[1]), float(partes[2])]
            )
    return aceleraciones


def error_relativo(referencia, prueba):
    """Error relativo porcentual por particula: (maximo, promedio, n)."""
    if len(referencia) != len(prueba):
        raise ValueError(
            "cantidad de aceleraciones distinta: %d vs %d"
            % (len(referencia), len(prueba))
        )

    maximo = 0.0
    suma = 0.0
    usadas = 0

    for ref, tes in zip(referencia, prueba):
        diferencia = math.sqrt(
            sum((tes[eje] - ref[eje]) ** 2 for eje in range(3))
        )
        magnitud = math.sqrt(sum(valor ** 2 for valor in ref))
        if magnitud < 1e-300:
            continue
        error = 100.0 * diferencia / magnitud
        maximo = max(maximo, error)
        suma += error
        usadas += 1

    promedio = suma / usadas if usadas else 0.0
    return maximo, promedio, usadas


def ejecutar(programa, ruta_particulas, n, args, salida_acc):
    # Se acepta tanto --programa python3 sim.py como
    # --programa "python3 sim.py".
    if len(programa) == 1 and " " in programa[0]:
        comando = shlex.split(programa[0])
    else:
        comando = list(programa)

    comando += [
        ruta_particulas,
        str(n),
        str(args.pasos),
        str(args.dt),
        str(args.theta),
        str(args.eps),
        str(args.hilos),
        salida_acc,
    ]

    print("Ejecutando:", " ".join(comando))

    resultado = subprocess.run(
        comando,
        capture_output=True,
        text=True,
        cwd=RAIZ,
    )

    if resultado.returncode != 0:
        print(resultado.stdout)
        print(resultado.stderr, file=sys.stderr)
        raise SystemExit(
            "el programa termino con codigo %d" % resultado.returncode
        )

    claves = {}
    for linea in resultado.stdout.splitlines():
        linea = linea.strip()
        if "=" in linea:
            clave, valor = linea.split("=", 1)
            try:
                claves[clave.strip()] = float(valor)
            except ValueError:
                pass

    if args.verbose:
        print(resultado.stdout.rstrip())

    return claves


def resumen_energia(claves):
    inicial = claves.get("energia_inicial")
    final = claves.get("energia_final")
    if inicial is None or final is None or abs(inicial) < 1e-300:
        return None
    return 100.0 * abs(final - inicial) / abs(inicial)


def modo_comparar(ruta_a, ruta_b):
    a = leer_aceleraciones(ruta_a)
    b = leer_aceleraciones(ruta_b)
    maximo, promedio, n = error_relativo(a, b)

    print("COMPARACION DE ACELERACIONES")
    print("  archivo A:", ruta_a, "(%d)" % len(a))
    print("  archivo B:", ruta_b, "(%d)" % len(b))
    print("  error maximo:   %.6f %%" % maximo)
    print("  error promedio: %.6f %%" % promedio)
    print("  particulas:", n)

    return 0 if len(a) == len(b) else 1


def main():
    parser = argparse.ArgumentParser(
        description="Verifica una version contra la suma directa."
    )
    parser.add_argument("--particulas", help="archivo de N filas y 7 columnas")
    parser.add_argument("--N", type=int, default=2000,
                        help="tamano si se generan particulas (defecto 2000)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--programa", nargs="+",
                        help="comando de la version a verificar")
    parser.add_argument("--pasos", type=int, default=5)
    parser.add_argument("--dt", type=float, default=0.001)
    parser.add_argument("--theta", type=float, default=0.3)
    parser.add_argument("--eps", type=float, default=0.01)
    parser.add_argument("--hilos", type=int, default=1)
    parser.add_argument("--tol-error", type=float, default=1.0,
                        help="error relativo promedio maximo en porcentaje")
    parser.add_argument("--tol-energia", type=float, default=1.0,
                        help="derivacion de energia maxima en porcentaje")
    parser.add_argument("--sin-energia", action="store_true",
                        help="no exige el control de energia")
    parser.add_argument("--comparar", nargs=2, metavar=("A", "B"),
                        help="compara dos archivos de aceleraciones")
    parser.add_argument("--verbose", action="store_true")

    args = parser.parse_args()

    if args.comparar:
        raise SystemExit(modo_comparar(*args.comparar))

    temporal = None

    if args.particulas:
        ruta = args.particulas
        particulas = cargar_particulas(ruta)
    else:
        particulas = generar_particulas(args.N, args.seed)
        descriptor, ruta = tempfile.mkstemp(
            prefix="particulas_", suffix=".txt"
        )
        os.close(descriptor)
        escribir_particulas(particulas, ruta)
        temporal = ruta

    n = len(particulas)
    print("VERIFICACION")
    print("  particulas:", n)
    print("  theta:", args.theta, " dt:", args.dt,
          " pasos:", args.pasos, " hilos:", args.hilos)

    try:
        # 1. Aceleraciones de la version bajo prueba.
        if args.programa:
            salida_acc = os.path.join(
                tempfile.gettempdir(), "verificar_aceleraciones.txt"
            )
            if os.path.exists(salida_acc):
                os.remove(salida_acc)
            claves = ejecutar(
                args.programa, ruta, n, args, salida_acc
            )
            if not os.path.exists(salida_acc):
                raise SystemExit(
                    "el programa no escribio %s" % salida_acc
                )
            aceleraciones = leer_aceleraciones(salida_acc)
        else:
            claves = {}
            aceleraciones = calcular_todas_las_aceleraciones_barnes_hut(
                particulas, theta=args.theta, epsilon=args.eps
            )

        # 2. Error contra la suma directa.
        if n > N_MAXIMO_DIRECTA:
            print("  [omitido] suma directa con N > %d" % N_MAXIMO_DIRECTA)
            error_max = error_prom = 0.0
        else:
            directas = calcular_todas_las_aceleraciones(particulas)
            error_max, error_prom, usadas = error_relativo(
                directas, aceleraciones
            )
            print("  error maximo:   %.6f %%" % error_max)
            print("  error promedio: %.6f %%  (%d particulas)"
                  % (error_prom, usadas))

        # 3. Derivacion de la energia.
        if args.sin_energia:
            drift = None
        elif args.programa:
            drift = resumen_energia(claves)
            if drift is None:
                print("  [aviso] el programa no reporto energia_inicial/"
                      "energia_final")
        elif n <= N_MAXIMO_ENERGIA:
            from referencia.motor.fuerzas_barnes_hut import (
                paso_leapfrog_barnes_hut,
            )

            inicial = calcular_energia_total(particulas)
            for _ in range(args.pasos):
                paso_leapfrog_barnes_hut(
                    particulas, args.dt, args.theta,
                    epsilon=args.eps,
                )
            final = calcular_energia_total(particulas)
            drift = (100.0 * abs(final - inicial) / abs(inicial)
                     if abs(inicial) > 1e-300 else 0.0)
        else:
            drift = None
            print("  [omitido] energia con N > %d" % N_MAXIMO_ENERGIA)

        if drift is not None:
            print("  derivacion de energia: %.6f %%" % drift)

        # 4. Veredicto.
        fallos = []
        if n <= N_MAXIMO_DIRECTA and error_prom > args.tol_error:
            fallos.append(
                "error promedio %.4f %% > %.4f %%"
                % (error_prom, args.tol_error)
            )
        if drift is not None and drift > args.tol_energia:
            fallos.append(
                "deriva de energia %.4f %% > %.4f %%"
                % (drift, args.tol_energia)
            )

        print()
        if fallos:
            print("FALLO:")
            for fallo in fallos:
                print("  -", fallo)
            raise SystemExit(1)

        print("OK: la version pasa la verificacion.")
        raise SystemExit(0)

    finally:
        if temporal and os.path.exists(temporal):
            os.remove(temporal)


if __name__ == "__main__":
    main()
