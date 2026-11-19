#!/usr/bin/env python3
"""
Runner de cierre de la version Python (multiprocessing).

Compara las dos estrategias ---compartido y local--- midiendo el escalado
con 1, 2, 4, 8 (y opcionalmente 16) procesos, con varias repeticiones (se
toma la MEDIANA), y calcula speedup y eficiencia respecto al tiempo
secuencial (1 proceso).

El baseline (1 proceso) usa la ruta secuencial de sim.py y es COMUN a las
dos estrategias: es lo correcto para un analisis de escalado honesto, y
resuelve que la estrategia 'local' no tenga su propia fila de 1 proceso.

Deja:
  * un CSV consolidado en resultados/ (una fila por estrategia/hilos), y
  * una tabla comparativa por pantalla, lista para el informe.

Uso:
  python3 tools/medir_cierre.py --archivo data/p20k.txt --N 20000 \
      --pasos 3 --reps 3 --hilos 1 2 4 8
"""

import argparse
import os
import statistics
import subprocess
import sys
from datetime import datetime

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SIM = [sys.executable, "python_mp/sim.py"]
CLAVES = ("tiempo_build", "tiempo_walk", "tiempo_integrate", "tiempo_total")


def correr_una(archivo, n, pasos, dt, theta, eps, hilos, estrategia):
    """Corre sim.py una vez y devuelve los tiempos por fase (segundos)."""

    env = dict(os.environ)
    env["BH_ESTRATEGIA"] = estrategia
    salida_acc = "/tmp/cierre_acc_%d.txt" % os.getpid()

    cmd = SIM + [
        archivo, str(n), str(pasos), str(dt),
        str(theta), str(eps), str(hilos), salida_acc,
    ]
    res = subprocess.run(cmd, cwd=RAIZ, env=env,
                         capture_output=True, text=True)
    if res.returncode != 0:
        sys.stderr.write(res.stderr)
        raise SystemExit(
            "fallo la corrida: hilos=%d estrategia=%s" % (hilos, estrategia)
        )

    claves = {}
    for linea in res.stdout.splitlines():
        if "=" in linea:
            clave, valor = linea.split("=", 1)
            claves[clave.strip()] = valor.strip()

    try:
        os.remove(salida_acc)
    except OSError:
        pass

    return {c: float(claves[c]) for c in CLAVES if c in claves}


def mediana(reps, clave):
    return statistics.median(r[clave] for r in reps)


def fila(estrategia, hilos, reps, base_total):
    total = mediana(reps, "tiempo_total")
    speedup = base_total / total if total > 0 else 0.0
    return {
        "estrategia": estrategia,
        "hilos": hilos,
        "build": mediana(reps, "tiempo_build"),
        "walk": mediana(reps, "tiempo_walk"),
        "integrate": mediana(reps, "tiempo_integrate"),
        "total": total,
        "speedup": speedup,
        "eficiencia": speedup / hilos,
    }


def main():
    p = argparse.ArgumentParser(
        description="Runner de cierre: compartido vs local"
    )
    p.add_argument("--archivo", required=True)
    p.add_argument("--N", type=int, required=True)
    p.add_argument("--pasos", type=int, default=3)
    p.add_argument("--dt", type=float, default=0.001)
    p.add_argument("--theta", type=float, default=0.3)
    p.add_argument("--eps", type=float, default=0.01)
    p.add_argument("--reps", type=int, default=3)
    p.add_argument("--hilos", type=int, nargs="+", default=[1, 2, 4, 8])
    p.add_argument("--estrategias", nargs="+",
                   default=["compartido", "local"])
    p.add_argument("--out", default=None)
    args = p.parse_args()

    def correr(hilos, estrategia):
        return [
            correr_una(args.archivo, args.N, args.pasos, args.dt,
                       args.theta, args.eps, hilos, estrategia)
            for _ in range(args.reps)
        ]

    # 1) Baseline secuencial (1 proceso), comun a las dos estrategias.
    sys.stderr.write(
        "# baseline secuencial (1 proceso) x%d...\n" % args.reps
    )
    base = correr(1, "compartido")
    base_total = mediana(base, "tiempo_total")

    filas = [{
        "estrategia": "secuencial", "hilos": 1,
        "build": mediana(base, "tiempo_build"),
        "walk": mediana(base, "tiempo_walk"),
        "integrate": mediana(base, "tiempo_integrate"),
        "total": base_total, "speedup": 1.0, "eficiencia": 1.0,
    }]

    # 2) Cada estrategia con los conteos > 1.
    for estrategia in args.estrategias:
        for hilos in [h for h in args.hilos if h > 1]:
            sys.stderr.write(
                "# %s / %d proceso(s) x%d...\n"
                % (estrategia, hilos, args.reps)
            )
            filas.append(
                fila(estrategia, hilos, correr(hilos, estrategia),
                     base_total)
            )

    # 3) CSV consolidado.
    marca = datetime.now().strftime("%Y%m%d_%H%M%S")
    out = args.out or "resultados/cierre_python_mp_%s.csv" % marca
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    with open(out, "w") as f:
        f.write("estrategia,hilos,tiempo_build,tiempo_walk,"
                "tiempo_integrate,tiempo_total,speedup,eficiencia\n")
        for r in filas:
            f.write("%s,%d,%.6f,%.6f,%.6f,%.6f,%.4f,%.4f\n" % (
                r["estrategia"], r["hilos"], r["build"], r["walk"],
                r["integrate"], r["total"], r["speedup"], r["eficiencia"],
            ))

    # 4) Tabla por pantalla.
    print()
    print("N=%d  pasos=%d  theta=%s  eps=%s  reps=%d (mediana)"
          % (args.N, args.pasos, args.theta, args.eps, args.reps))
    cab = ("estrategia", "hilos", "build", "walk", "integr",
           "total", "speedup", "eficien")
    print("%-11s %5s %8s %8s %8s %8s %8s %8s" % cab)
    for r in filas:
        print("%-11s %5d %8.3f %8.3f %8.3f %8.3f %8.2f %8.2f" % (
            r["estrategia"], r["hilos"], r["build"], r["walk"],
            r["integrate"], r["total"], r["speedup"], r["eficiencia"],
        ))
    print()
    print("CSV: %s" % out)


if __name__ == "__main__":
    main()
