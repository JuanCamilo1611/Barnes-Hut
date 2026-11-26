#!/usr/bin/env python3
"""
Figuras de escalado de la version Python (compartido vs local).

Lee un CSV consolidado de tools/medir_cierre.py (columnas: estrategia,
hilos, tiempo_build, tiempo_walk, tiempo_integrate, tiempo_total, speedup,
eficiencia) y genera dos figuras en resultados/plots/:

  1. escalado_python_mp.png   speedup y eficiencia vs nro de procesos,
                              con la linea ideal y la marca de nucleos
                              fisicos.
  2. fases_python_mp.png      tiempo de build y total por proceso, que
                              muestra el costo del rebuild por proceso en
                              la estrategia local.

Uso:
  python3 tools/graficas.py --csv resultados/cierre_python_mp_XXXX.csv \
      --nucleos-fisicos 4
"""

import argparse
import csv
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Paleta categorica CVD-safe (Okabe-Ito): identidad por estrategia.
COLOR = {"compartido": "#0072B2", "local": "#D55E00"}
GRIS = "#9A9A9A"
TINTA = "#1A1A1A"
REJILLA = "#D9D9D9"


def leer(csv_path):
    filas = []
    with open(csv_path) as f:
        for r in csv.DictReader(f):
            filas.append({
                "estrategia": r["estrategia"],
                "hilos": int(r["hilos"]),
                "build": float(r["tiempo_build"]),
                "walk": float(r["tiempo_walk"]),
                "total": float(r["tiempo_total"]),
                "speedup": float(r["speedup"]),
                "eficiencia": float(r["eficiencia"]),
            })
    return filas


def serie(filas, estrategia):
    """Puntos de una estrategia incluyendo el baseline secuencial (1)."""
    base = [r for r in filas if r["estrategia"] == "secuencial"]
    propios = [r for r in filas if r["estrategia"] == estrategia]
    return sorted(base + propios, key=lambda r: r["hilos"])


def estilo_ejes(ax):
    ax.grid(True, color=REJILLA, linewidth=0.8, alpha=0.7)
    ax.set_axisbelow(True)
    for lado in ("top", "right"):
        ax.spines[lado].set_visible(False)
    ax.tick_params(colors=TINTA)


def fig_escalado(filas, nucleos, salida):
    hilos = sorted({r["hilos"] for r in filas})
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.2))

    # --- speedup ---
    ax1.plot(hilos, hilos, "--", color=GRIS, linewidth=1.5,
             label="ideal (lineal)", zorder=1)
    # etiqueta directa en el ultimo punto; la de 'local' va debajo para
    # que no se encime cuando las dos curvas terminan casi en el mismo valor
    desfase = {"compartido": (6, 6), "local": (6, -12)}
    for estr in ("compartido", "local"):
        pts = serie(filas, estr)
        xs = [p["hilos"] for p in pts]
        ys = [p["speedup"] for p in pts]
        ax1.plot(xs, ys, "-o", color=COLOR[estr], linewidth=2,
                 markersize=7, label=estr, zorder=3)
        ax1.annotate("%.2f x" % ys[-1], (xs[-1], ys[-1]),
                     textcoords="offset points", xytext=desfase[estr],
                     color=COLOR[estr], fontsize=9)
    ax1.axvline(nucleos, color=GRIS, linewidth=1, alpha=0.6)
    ax1.annotate("%d nucleos\nfisicos" % nucleos, (nucleos, 2.8),
                 color=TINTA, fontsize=8, ha="center",
                 bbox=dict(boxstyle="round,pad=0.2", fc="white",
                           ec="none", alpha=0.85))
    ax1.set_title("Speedup vs procesos", color=TINTA, fontsize=12)
    ax1.set_xlabel("procesos", color=TINTA)
    ax1.set_ylabel("speedup", color=TINTA)
    ax1.set_xticks(hilos)
    estilo_ejes(ax1)
    ax1.legend(frameon=False, fontsize=9)

    # --- eficiencia ---
    ax2.axhline(1.0, linestyle="--", color=GRIS, linewidth=1.5,
                label="ideal (1.0)", zorder=1)
    for estr in ("compartido", "local"):
        pts = serie(filas, estr)
        xs = [p["hilos"] for p in pts]
        ys = [p["eficiencia"] for p in pts]
        ax2.plot(xs, ys, "-o", color=COLOR[estr], linewidth=2,
                 markersize=7, label=estr, zorder=3)
    ax2.axvline(nucleos, color=GRIS, linewidth=1, alpha=0.6)
    ax2.set_title("Eficiencia vs procesos", color=TINTA, fontsize=12)
    ax2.set_xlabel("procesos", color=TINTA)
    ax2.set_ylabel("eficiencia", color=TINTA)
    ax2.set_xticks(hilos)
    ax2.set_ylim(0, 1.1)
    estilo_ejes(ax2)
    ax2.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    fig.savefig(salida, dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_fases(filas, salida):
    paralelos = sorted({r["hilos"] for r in filas if r["hilos"] > 1})
    x = range(len(paralelos))
    ancho = 0.38
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.2))

    def barras(ax, clave, titulo):
        for i, estr in enumerate(("compartido", "local")):
            vals = []
            for h in paralelos:
                r = next((r for r in filas
                          if r["estrategia"] == estr and r["hilos"] == h),
                         None)
                vals.append(r[clave] if r else 0.0)
            pos = [xi + (i - 0.5) * ancho for xi in x]
            ax.bar(pos, vals, ancho, color=COLOR[estr], label=estr)
            for p, v in zip(pos, vals):
                ax.annotate(("%.2f" if v < 1 else "%.1f") % v, (p, v),
                            textcoords="offset points", xytext=(0, 3),
                            ha="center", color=TINTA, fontsize=8)
        ax.set_title(titulo, color=TINTA, fontsize=12)
        ax.set_xlabel("procesos", color=TINTA)
        ax.set_ylabel("segundos", color=TINTA)
        ax.set_xticks(list(x))
        ax.set_xticklabels([str(h) for h in paralelos])
        estilo_ejes(ax1)
        estilo_ejes(ax2)
        ax.legend(frameon=False, fontsize=9)

    barras(ax1, "build", "Tiempo de construccion del arbol")
    barras(ax2, "total", "Tiempo total")
    fig.tight_layout()
    fig.savefig(salida, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description="Figuras de escalado Python")
    p.add_argument("--csv", default=None,
                   help="CSV consolidado (por defecto el mas reciente)")
    p.add_argument("--nucleos-fisicos", type=int, default=4)
    p.add_argument("--dir-salida", default="resultados/plots")
    args = p.parse_args()

    csv_path = args.csv
    if csv_path is None:
        candidatos = sorted(glob.glob("resultados/cierre_python_mp_*.csv"))
        if not candidatos:
            raise SystemExit("no hay CSV de cierre en resultados/")
        csv_path = candidatos[-1]

    filas = leer(csv_path)
    os.makedirs(args.dir_salida, exist_ok=True)

    f1 = os.path.join(args.dir_salida, "escalado_python_mp.png")
    f2 = os.path.join(args.dir_salida, "fases_python_mp.png")
    fig_escalado(filas, args.nucleos_fisicos, f1)
    fig_fases(filas, f2)

    print("CSV leido: %s" % csv_path)
    print("Figura 1 (speedup/eficiencia): %s" % f1)
    print("Figura 2 (fases build/total):  %s" % f2)


if __name__ == "__main__":
    main()
