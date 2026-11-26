#!/usr/bin/env python3
"""
Runner de medicion para las versiones en C (secuencial y MPI).

Corre cada configuracion varias veces, toma la MEDIANA de cada tiempo y
calcula dos speedups:
  speedup                contra MPI con 1 proceso: mismo codigo, mide el
                         escalado (es el que se usa en el analisis);
  speedup_vs_secuencial  contra el C secuencial (c/common/sim_barnes_hut),
                         que recorre en el orden del archivo: incluye la
                         ganancia de localidad del orden Morton de MPI.

Colocacion de procesos MPI en un i5-8250U (4 nucleos fisicos, 8 hilos):
  P <= 4   opciones por defecto: un proceso por NUCLEO fisico.
  P == 8   --use-hwthread-cpus: un proceso por hilo logico (2 por nucleo).
  P > 8    --oversubscribe --bind-to none: mas procesos que hilos; el
           sistema operativo los reparte.
(Con --use-hwthread-cpus y P=4 Open MPI usa los hilos 0-3 = solo 2
nucleos fisicos, por eso no se usa para P <= 4.)

Uso:
  python3 tools/medir.py --archivo /tmp/u200k.txt --N 200000 \
      --pasos 3 --reps 3 --procesos 1 2 4 8 16 --etiqueta uniforme_200k

Deja un CSV en resultados/ y una tabla por pantalla.
"""

import argparse
import csv
import os
import statistics
import subprocess
import sys
import tempfile
from datetime import datetime

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SECUENCIAL = "c/common/sim_barnes_hut"
MPI = "c/mpi/sim_mpi"
NUCLEOS_FISICOS = 4
HILOS_LOGICOS = 8

COLUMNAS = ("tiempo_total", "tiempo_walk", "tiempo_build", "tiempo_sort",
            "tiempo_reparto", "tiempo_let", "tiempo_integrate",
            "desbalance_walk", "nodos_arbol_local_max",
            "nodos_let_recibidos_max", "interacciones_por_particula")


def opciones_mpi(p):
    if p <= NUCLEOS_FISICOS:
        return []
    if p <= HILOS_LOGICOS:
        return ["--use-hwthread-cpus"]
    return ["--oversubscribe", "--bind-to", "none"]


def correr(comando, env):
    r = subprocess.run(comando, capture_output=True, text=True, cwd=RAIZ,
                       env=env)
    if r.returncode != 0:
        sys.stderr.write(r.stdout + r.stderr)
        raise SystemExit("fallo: %s" % " ".join(comando))
    valores = {}
    for linea in r.stdout.splitlines():
        if "=" in linea:
            k, v = linea.split("=", 1)
            try:
                valores[k.strip()] = float(v)
            except ValueError:
                pass
    return valores


def mediana(corridas):
    claves = set().union(*corridas)
    return {k: statistics.median(c[k] for c in corridas if k in c)
            for k in claves}


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--archivo", required=True)
    p.add_argument("--N", type=int, required=True)
    p.add_argument("--pasos", type=int, default=3)
    p.add_argument("--dt", type=float, default=0.001)
    p.add_argument("--theta", type=float, default=0.5)
    p.add_argument("--eps", type=float, default=0.01)
    p.add_argument("--reps", type=int, default=3)
    p.add_argument("--procesos", type=int, nargs="+", default=[1, 2, 4, 8, 16])
    p.add_argument("--etiqueta", default="")
    p.add_argument("--out", default=None)
    a = p.parse_args()

    salida_acc = os.path.join(tempfile.gettempdir(), "medir_acc.txt")
    comun = [a.archivo, str(a.N), str(a.pasos), str(a.dt), str(a.theta),
             str(a.eps)]
    env = dict(os.environ)

    filas = []
    print("# secuencial x%d..." % a.reps, file=sys.stderr)
    base = mediana([correr([SECUENCIAL] + comun + ["1", salida_acc], env)
                    for _ in range(a.reps)])
    base.update(version="secuencial", procesos=1)
    filas.append(base)

    for np_ in a.procesos:
        print("# mpi P=%d x%d..." % (np_, a.reps), file=sys.stderr)
        cmd = (["mpirun"] + opciones_mpi(np_) + ["-np", str(np_), MPI]
               + comun + [str(np_), salida_acc])
        m = mediana([correr(cmd, env) for _ in range(a.reps)])
        m.update(version="mpi", procesos=np_)
        filas.append(m)

    # speedup: contra MPI con 1 proceso (mismo codigo; mide el escalado).
    # speedup_vs_secuencial: contra el C secuencial; incluye ademas la
    # ganancia de recorrer en orden Morton (localidad de cache), que el
    # secuencial no hace.
    t_sec = base["tiempo_total"]
    t_p1 = next((f["tiempo_total"] for f in filas
                 if f["version"] == "mpi" and f["procesos"] == 1), t_sec)
    for f in filas:
        ref = t_sec if f["version"] == "secuencial" else t_p1
        f["speedup"] = ref / f["tiempo_total"]
        f["eficiencia"] = f["speedup"] / f["procesos"]
        f["speedup_vs_secuencial"] = t_sec / f["tiempo_total"]

    sello = datetime.now().strftime("%Y%m%d_%H%M%S")
    nombre = a.out or os.path.join(
        RAIZ, "resultados",
        "medir_%s_%s.csv" % (a.etiqueta or "mpi", sello))
    campos = (["version", "procesos", "speedup", "eficiencia",
               "speedup_vs_secuencial"] + list(COLUMNAS))
    with open(nombre, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=campos, extrasaction="ignore")
        w.writeheader()
        for f in filas:
            w.writerow({k: (round(v, 6) if isinstance(v, float) else v)
                        for k, v in f.items()})

    print("\nN=%d pasos=%d theta=%g reps=%d (mediana)  %s"
          % (a.N, a.pasos, a.theta, a.reps, a.etiqueta))
    print("%-10s %4s %8s %8s %8s %8s %8s %7s %7s %9s %9s" % (
        "version", "P", "total", "walk", "build", "reparto", "let",
        "speedup", "efic", "desbal", "let_nodos"))
    for f in filas:
        print("%-10s %4d %8.3f %8.3f %8.3f %8.3f %8.3f %7.2f %7.2f %9.3f %9d" % (
            f["version"], f["procesos"], f["tiempo_total"], f["tiempo_walk"],
            f["tiempo_build"], f.get("tiempo_reparto", 0.0),
            f.get("tiempo_let", 0.0), f["speedup"], f["eficiencia"],
            f.get("desbalance_walk", 1.0),
            int(f.get("nodos_let_recibidos_max", 0))))
    print("\nCSV: %s" % os.path.relpath(nombre, RAIZ))


if __name__ == "__main__":
    main()
