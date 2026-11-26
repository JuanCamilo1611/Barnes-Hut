#!/usr/bin/env python3
"""
Comparacion de las versiones del proyecto con los mismos datos.

Corre Python-mp, C secuencial, OpenMP y MPI sobre las mismas
distribuciones y tamanos, con 1, 2, 4 y 8 procesos o hilos, y guarda la
MEDIANA de varias repeticiones en un CSV. CUDA se mide en Colab (ver
cuda/barnes_hut_cuda.ipynb) y se agrega al mismo CSV con --cuda.

Metrica: tiempo por evaluacion de fuerzas = tiempo / (pasos + 1), porque
todas las versiones calculan las fuerzas iniciales y una vez por paso.
El tiempo es la suma de las fases (sort + build + walk + integrate)
en todas las versiones salvo MPI, cuyo total es el de pared del proceso
mas lento sin E/S (incluye reparto y LET). En Python se usa la suma de
fases porque su tiempo_total incluye escribir el archivo de aceleraciones.

Uso:
  python3 tools/comparar.py --datos ~/.bh_datos --reps 3       # medir
  python3 tools/comparar.py --cuda resultados/cuda.csv         # agregar CUDA
  python3 tools/comparar.py --graficar resultados/comparar_X.csv

Los archivos de datos se esperan como <datos>/<distribucion>_<N>.txt y se
generan con data/generate_particles.py si no existen.
"""

import argparse
import csv
import os
import statistics
import subprocess
import sys
from datetime import datetime

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DISTRIBUCIONES = ("uniforme", "cumulos", "concentradas")
PASOS, DT, THETA, EPS = 1, 0.001, 0.5, 0.01
CAMPOS = ("version", "distribucion", "N", "procesos", "reps",
          "t_eval", "t_build", "t_sort", "t_walk", "t_comunicacion")


def comando(version, archivo, n, p):
    """Linea de comando y variables de entorno de cada version."""
    args = [archivo, str(n), str(PASOS), str(DT), str(THETA), str(EPS), str(p),
            "/tmp/comparar_acc.txt"]
    env = dict(os.environ)
    if version == "python_mp":
        return ["python3", "python_mp/sim.py"] + args, env
    if version == "c_secuencial":
        return ["c/common/sim_barnes_hut"] + args, env
    if version == "openmp":
        env["OMP_SCHEDULE"] = "dynamic,64"
        return ["c/openmp/sim_openmp"] + args, env
    if version == "mpi":
        # Misma colocacion que tools/medir.py (i5-8250U: 4 nucleos, 8 hilos).
        if p <= 4:
            extra = []
        elif p == 8:
            extra = ["--use-hwthread-cpus"]
        else:
            extra = ["--oversubscribe", "--bind-to", "none"]
        return ["mpirun", "-np", str(p)] + extra + ["c/mpi/sim_mpi"] + args, env
    raise ValueError(version)


def correr(version, archivo, n, p):
    cmd, env = comando(version, archivo, n, p)
    out = subprocess.run(cmd, cwd=RAIZ, env=env, capture_output=True,
                         text=True, timeout=3600)
    if out.returncode != 0:
        raise RuntimeError(f"{' '.join(cmd)}\n{out.stderr[-2000:]}")
    v = {}
    for linea in out.stdout.splitlines():
        if "=" in linea:
            k, _, x = linea.partition("=")
            try:
                v[k.strip()] = float(x)
            except ValueError:
                pass
    fases = v["tiempo_build"] + v["tiempo_sort"] + v["tiempo_walk"] + v["tiempo_integrate"]
    total = v["tiempo_total"] if version == "mpi" else fases
    com = v.get("tiempo_reparto", 0.0) + v.get("tiempo_let", 0.0)
    e = PASOS + 1
    return {"t_eval": total / e, "t_build": v["tiempo_build"] / e,
            "t_sort": v["tiempo_sort"] / e, "t_walk": v["tiempo_walk"] / e,
            "t_comunicacion": com / e}


def plan(args):
    """(version, N, procesos, reps) a medir."""
    filas = []
    for n in args.N:
        filas.append(("c_secuencial", n, 1, args.reps))
        for p in args.procesos:
            filas.append(("openmp", n, p, args.reps))
            filas.append(("mpi", n, p, args.reps))
        # Python es ~100x mas lento: con N grande solo 1 y 8 procesos, 1 rep.
        if n <= 20000:
            for p in args.procesos:
                filas.append(("python_mp", n, p, args.reps))
        elif not args.sin_python_grande:
            for p in (1, max(args.procesos)):
                filas.append(("python_mp", n, p, 1))
    return filas


def medir(args):
    os.makedirs(os.path.join(RAIZ, "resultados"), exist_ok=True)
    salida = os.path.join(RAIZ, "resultados",
                          "comparar_" + datetime.now().strftime("%Y%m%d_%H%M") + ".csv")
    with open(salida, "w", newline="") as fh:
        w = csv.DictWriter(fh, CAMPOS)
        w.writeheader()
        for n in args.N:
            for d in DISTRIBUCIONES:
                archivo = os.path.join(args.datos, f"{d}_{n}.txt")
                if not os.path.exists(archivo):
                    subprocess.run(["python3", "data/generate_particles.py",
                                    "--distribucion", d, "--N", str(n),
                                    "--out", archivo], cwd=RAIZ, check=True,
                                   capture_output=True)
                for version, nn, p, reps in plan(args):
                    if nn != n:
                        continue
                    medidas = [correr(version, archivo, n, p) for _ in range(reps)]
                    fila = {"version": version, "distribucion": d, "N": n,
                            "procesos": p, "reps": reps}
                    for k in medidas[0]:
                        fila[k] = round(statistics.median(m[k] for m in medidas), 6)
                    w.writerow(fila)
                    fh.flush()
                    print(f"{version:13s} {d:13s} N={n:<7d} P={p:<2d} "
                          f"t_eval={fila['t_eval']:.4f} s", flush=True)
    print("Guardado en", salida)


def graficar(ruta):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with open(ruta) as fh:
        filas = list(csv.DictReader(fh))
    for f in filas:
        f["N"] = int(f["N"]); f["procesos"] = int(f["procesos"])
        f["t_eval"] = float(f["t_eval"])

    estilo = {  # Okabe-Ito + marcadores distintos (legible sin color)
        "python_mp":    ("Python-mp",    "#E69F00", "o"),
        "c_secuencial": ("C secuencial", "#000000", "s"),
        "openmp":       ("OpenMP",       "#0072B2", "^"),
        "mpi":          ("MPI",          "#009E73", "D"),
        "cuda":         ("CUDA (T4)",    "#CC79A7", "*"),
    }
    salida = os.path.join(RAIZ, "resultados", "plots")
    os.makedirs(salida, exist_ok=True)

    for n in sorted({f["N"] for f in filas}):
        fig, ejes = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True)
        for ax, d in zip(ejes, DISTRIBUCIONES):
            for v, (nombre, color, m) in estilo.items():
                pts = sorted((f["procesos"], f["t_eval"]) for f in filas
                             if f["version"] == v and f["N"] == n
                             and f["distribucion"] == d)
                if not pts:
                    continue
                if v in ("c_secuencial", "cuda"):   # una sola medida: linea
                    # CUDA: 'procesos' son hilos por bloque; se usa el mejor.
                    ax.axhline(min(t for _, t in pts), color=color, ls="--",
                               lw=1.4, label=nombre)
                    continue
                ax.plot([p for p, _ in pts], [t for _, t in pts], marker=m,
                        color=color, lw=2, label=nombre)
            ax.set_xscale("log", base=2)
            ax.set_yscale("log")
            ax.set_xticks([1, 2, 4, 8]); ax.set_xticklabels(["1", "2", "4", "8"])
            ax.set_title(d)
            ax.set_xlabel("procesos / hilos")
            ax.grid(True, which="major", alpha=0.3)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
        ejes[0].set_ylabel("segundos por evaluacion de fuerzas")
        ejes[0].legend(frameon=False, fontsize=9)
        fig.suptitle(f"Comparacion de las versiones, N = {n:,}".replace(",", " "))
        fig.tight_layout()
        destino = os.path.join(salida, f"comparar_N{n}.png")
        fig.savefig(destino, dpi=150)
        print("Figura:", destino)


def agregar_cuda(ruta_cuda, ruta_csv):
    """Anade filas de CUDA (mismo formato, version=cuda) al CSV indicado."""
    with open(ruta_cuda) as fh:
        nuevas = list(csv.DictReader(fh))
    with open(ruta_csv, "a", newline="") as fh:
        w = csv.DictWriter(fh, CAMPOS)
        for f in nuevas:
            w.writerow({k: f.get(k, "") for k in CAMPOS})
    print(f"{len(nuevas)} filas de CUDA agregadas a {ruta_csv}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--datos", default=os.path.expanduser("~/.bh_datos"))
    ap.add_argument("--N", type=int, nargs="+", default=[20000, 200000])
    ap.add_argument("--procesos", type=int, nargs="+", default=[1, 2, 4, 8])
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--sin-python-grande", action="store_true",
                    help="no medir Python con N > 20000")
    ap.add_argument("--graficar", metavar="CSV")
    ap.add_argument("--cuda", metavar="CSV_CUDA",
                    help="agregar filas de CUDA a --destino")
    ap.add_argument("--destino", metavar="CSV")
    a = ap.parse_args()
    if a.graficar:
        graficar(a.graficar)
    elif a.cuda:
        if not a.destino:
            sys.exit("--cuda necesita --destino <csv de comparar>")
        agregar_cuda(a.cuda, a.destino)
    else:
        medir(a)
