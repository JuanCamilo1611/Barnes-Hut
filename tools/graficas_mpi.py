#!/usr/bin/env python3
"""
Figuras de la version MPI a partir de los CSV de tools/medir.py.

  1. escalado_mpi.png   speedup vs procesos (respecto a MPI con 1 proceso) para
                        cada distribucion, con la linea ideal y las zonas de
                        nucleos fisicos / hilos logicos / sobresuscripcion.
  2. let_mpi.png        (a) fraccion del arbol global que tiene cada proceso
                        (arbol local + LET recibidos) y (b) fraccion del
                        tiempo total en comunicacion (reparto + LET).

Uso:
  python3 tools/graficas_mpi.py            (toma el CSV mas reciente de
                                            cada etiqueta en resultados/)
"""

import csv
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SALIDA = os.path.join(RAIZ, "resultados", "plots")

# Okabe-Ito (validado con el validador de dataviz: CVD en la banda 6-8
# entre verde y rosa -> se compensa con marcadores distintos y etiquetas
# directas; el informe trae ademas la tabla).
SERIES = [
    ("uniforme_20k",      "uniforme 20k",      "#0072B2", "o"),
    ("uniforme_200k",     "uniforme 200k",     "#D55E00", "s"),
    ("cumulos_200k",      "cumulos 200k",      "#009E73", "^"),
    ("concentradas_200k", "concentradas 200k", "#CC79A7", "D"),
]
TINTA = "#222222"
GRIS = "#8a8a8a"
NUCLEOS, HILOS = 4, 8


def ultimo_csv(etiqueta):
    archivos = sorted(glob.glob(os.path.join(
        RAIZ, "resultados", "medir_%s_*.csv" % etiqueta)))
    return archivos[-1] if archivos else None


def leer(ruta):
    with open(ruta) as fh:
        filas = list(csv.DictReader(fh))
    for f in filas:
        for k, v in f.items():
            if k != "version":
                f[k] = float(v) if v not in ("", None) else 0.0
    # El escalado se mide contra MPI con 1 proceso (mismo codigo); se
    # recalcula aqui para no depender de la version del CSV.
    t1 = next(f["tiempo_total"] for f in filas
              if f["version"] == "mpi" and f["procesos"] == 1)
    for f in filas:
        if f["version"] == "mpi":
            f["speedup"] = t1 / f["tiempo_total"]
    return filas


def estilo(ax):
    for lado in ("top", "right"):
        ax.spines[lado].set_visible(False)
    for lado in ("left", "bottom"):
        ax.spines[lado].set_color(GRIS)
    ax.tick_params(colors=TINTA)
    ax.grid(True, color="#e6e6e6", linewidth=0.8)
    ax.set_axisbelow(True)


def zonas(ax, ymax):
    """Sombrea hilos logicos (4-8) y sobresuscripcion (>8)."""
    ax.axvspan(NUCLEOS, HILOS, color="#f2f2f2", zorder=0)
    ax.axvspan(HILOS, 16 * 1.15, color="#e4e4e4", zorder=0)
    ax.text(5.65, ymax * 0.97, "hilos\nlogicos", ha="center", va="top",
            fontsize=8, color=TINTA)
    ax.text(11.6, ymax * 0.97, "mas procesos\nque hilos", ha="center",
            va="top", fontsize=8, color=TINTA)


def main():
    datos = []
    for etiqueta, nombre, color, marca in SERIES:
        ruta = ultimo_csv(etiqueta)
        if ruta:
            datos.append((nombre, color, marca,
                          [f for f in leer(ruta) if f["version"] == "mpi"]))
            print("leido:", os.path.relpath(ruta, RAIZ))
    if not datos:
        raise SystemExit("no hay CSV de tools/medir.py en resultados/")
    os.makedirs(SALIDA, exist_ok=True)

    # ---------- Figura 1: speedup ----------
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    procesos = sorted({f["procesos"] for _, _, _, fs in datos for f in fs})
    ymax = max(5.0, max(f["speedup"] for _, _, _, fs in datos for f in fs) * 1.25)
    zonas(ax, ymax)
    ax.plot([1, NUCLEOS], [1, NUCLEOS], "--", color=GRIS, linewidth=1.5,
            label="ideal hasta 4 nucleos")
    for nombre, color, marca, fs in datos:
        xs = [f["procesos"] for f in fs]
        ys = [f["speedup"] for f in fs]
        ax.plot(xs, ys, "-", marker=marca, color=color, linewidth=2,
                markersize=7, label=nombre, zorder=3)
    ax.set_xscale("log", base=2)
    ax.set_xticks(procesos)
    ax.set_xticklabels([str(int(p)) for p in procesos])
    ax.set_xlim(0.85, 16 * 1.15)
    ax.set_ylim(0, ymax)
    ax.set_xlabel("procesos MPI", color=TINTA)
    ax.set_ylabel("speedup vs MPI con 1 proceso", color=TINTA)
    ax.set_title("Escalado de la version MPI", color=TINTA, fontsize=12)
    estilo(ax)
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(SALIDA, "escalado_mpi.png"), dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    # ---------- Figura 2: LET y comunicacion ----------
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4.2))
    for nombre, color, marca, fs in datos:
        global_ = next(f for f in fs if f["procesos"] == 1)["nodos_arbol_local_max"]
        varios = [f for f in fs if f["procesos"] > 1]
        xs = [f["procesos"] for f in varios]
        frac = [100.0 * (f["nodos_arbol_local_max"] + f["nodos_let_recibidos_max"])
                / global_ for f in varios]
        com = [100.0 * (f["tiempo_reparto"] + f["tiempo_let"]) / f["tiempo_total"]
               for f in fs]
        a1.plot(xs, frac, "-", marker=marca, color=color, linewidth=2,
                markersize=7, label=nombre)
        a2.plot([f["procesos"] for f in fs], com, "-", marker=marca,
                color=color, linewidth=2, markersize=7, label=nombre)
    a1.axhline(100, linestyle="--", color=GRIS, linewidth=1.5)
    a1.text(2.05, 102, "arbol global completo", fontsize=8, color=TINTA,
            va="bottom")
    a1.set_ylim(0, 115)
    a1.set_ylabel("% del arbol global (local + LET)", color=TINTA)
    a1.set_title("Arbol que tiene cada proceso", color=TINTA, fontsize=12)
    a2.set_ylabel("% del tiempo total", color=TINTA)
    a2.set_title("Comunicacion (reparto + LET)", color=TINTA, fontsize=12)
    for ax in (a1, a2):
        ax.set_xscale("log", base=2)
        ax.set_xticks(procesos)
        ax.set_xticklabels([str(int(p)) for p in procesos])
        ax.set_xlabel("procesos MPI", color=TINTA)
        estilo(ax)
    a2.set_ylim(bottom=0)
    con_let = [p for p in procesos if p > 1]     # con 1 proceso no hay LET
    a1.set_xticks(con_let)
    a1.set_xticklabels([str(int(p)) for p in con_let])
    a1.set_xlim(con_let[0] / 1.15, con_let[-1] * 1.15)
    a1.legend(frameon=False, fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(os.path.join(SALIDA, "let_mpi.png"), dpi=150,
                bbox_inches="tight")
    plt.close(fig)
    print("figuras en", os.path.relpath(SALIDA, RAIZ))


if __name__ == "__main__":
    main()
