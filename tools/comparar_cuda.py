#!/usr/bin/env python3
"""
Columna CUDA de la comparacion (correr en Colab, desde la raiz del repo,
despues de `cd cuda && make`).

Mide sim_cuda con los MISMOS datos y parametros que tools/comparar.py
(1 paso, dt 0.001, theta 0.5, eps 0.01, mediana de 3 repeticiones) y
deja resultados/cuda.csv con el mismo formato. Despues, en el portatil:

  python3 tools/comparar.py --cuda resultados/cuda.csv \
      --destino resultados/comparar_<fecha>.csv

Uso en una celda de Colab:
  !python3 tools/comparar_cuda.py               # 128 hilos por bloque
  !python3 tools/comparar_cuda.py --hilos 32 64 128 256   # barrido
"""

import argparse
import csv
import os
import statistics
import subprocess

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PASOS, DT, THETA, EPS = 1, 0.001, 0.5, 0.01
CAMPOS = ("version", "distribucion", "N", "procesos", "reps",
          "t_eval", "t_build", "t_sort", "t_walk", "t_comunicacion")


def correr(archivo, n, hilos):
    out = subprocess.run(
        ["./cuda/sim_cuda", archivo, str(n), str(PASOS), str(DT), str(THETA),
         str(EPS), str(hilos), "/tmp/acc_cuda.txt"],
        cwd=RAIZ, capture_output=True, text=True, check=True)
    v = {}
    for linea in out.stdout.splitlines():
        k, _, x = linea.partition("=")
        try:
            v[k] = float(x)
        except ValueError:
            pass
    e = PASOS + 1
    return {"t_eval": v["tiempo_total"] / e, "t_build": v["tiempo_build"] / e,
            "t_sort": v["tiempo_sort"] / e, "t_walk": v["tiempo_walk"] / e,
            "t_comunicacion": 0.0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, nargs="+", default=[20000, 200000])
    ap.add_argument("--hilos", type=int, nargs="+", default=[128])
    ap.add_argument("--reps", type=int, default=3)
    a = ap.parse_args()

    os.makedirs(os.path.join(RAIZ, "resultados"), exist_ok=True)
    salida = os.path.join(RAIZ, "resultados", "cuda.csv")
    with open(salida, "w", newline="") as fh:
        w = csv.DictWriter(fh, CAMPOS)
        w.writeheader()
        for n in a.N:
            for d in ("uniforme", "cumulos", "concentradas"):
                archivo = f"/tmp/{d}_{n}.txt"
                if not os.path.exists(archivo):
                    subprocess.run(["python3", "data/generate_particles.py",
                                    "--distribucion", d, "--N", str(n),
                                    "--out", archivo], cwd=RAIZ, check=True,
                                   capture_output=True)
                correr(archivo, n, a.hilos[0])          # calentamiento
                for h in a.hilos:
                    m = [correr(archivo, n, h) for _ in range(a.reps)]
                    fila = {"version": "cuda", "distribucion": d, "N": n,
                            "procesos": h, "reps": a.reps}
                    for k in m[0]:
                        fila[k] = round(statistics.median(x[k] for x in m), 6)
                    w.writerow(fila)
                    print(f"cuda {d:13s} N={n:<7d} hilos/bloque={h:<4d} "
                          f"t_eval={fila['t_eval']:.4f} s  "
                          f"(build {fila['t_build']:.4f}, sort {fila['t_sort']:.4f}, "
                          f"walk {fila['t_walk']:.4f})", flush=True)
    print("Guardado en", salida)


if __name__ == "__main__":
    main()
