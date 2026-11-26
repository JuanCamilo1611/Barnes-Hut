#!/usr/bin/env python3
"""
Diagnostico rapido de una distribucion de particulas (para la demo).

Responde, antes de correr nada, a: que tan concentrada esta, si hay
flujos que van a chocar y cuanto va a sufrir el reparto Morton de MPI.

Uso:
  python3 tools/diagnostico.py <archivo> [--procesos 4 8]
  python3 tools/diagnostico.py <archivo> --muestra 2000 --salida /tmp/m.txt
      (ademas escribe una muestra aleatoria para verificar contra la suma
       directa, que solo es viable con N pequeno)
"""

import argparse
import sys

import numpy as np

G = 32   # grilla de 32^3 celdas sobre el cubo raiz


def separar(v):
    x = v.astype(np.uint64) & np.uint64(0x1FFFFF)
    for s, m in [(32, 0x1F00000000FFFF), (16, 0x1F0000FF0000FF),
                 (8, 0x100F00F00F00F00F), (4, 0x10C30C30C30C30C3),
                 (2, 0x1249249249249249)]:
        x = (x | (x << np.uint64(s))) & np.uint64(m)
    return x


def cubo(pos):
    """Mismo cubo raiz que caja_envolvente() del nucleo C."""
    mn = pos.min(0)
    arista = float((pos.max(0) - mn).max()) or 1.0
    return mn, arista


def orden_morton(pos, mn, arista):
    q = np.clip(((pos - mn) / arista * (1 << 21)).astype(np.int64), 0, (1 << 21) - 1)
    c = (separar(q[:, 0]) << np.uint64(2)) | (separar(q[:, 1]) << np.uint64(1)) \
        | separar(q[:, 2])
    return np.argsort(c, kind="stable")


def vecinos_ajenos(celdas, P):
    """Particulas de otros procesos en celdas vecinas (aprox. del LET)."""
    n = len(celdas)
    tramos = np.array_split(np.arange(n), P)
    dueno = np.empty(n, int)
    for r, ix in enumerate(tramos):
        dueno[ix] = r
    res = []
    for r in range(P):
        occ = np.zeros((G + 2,) * 3, bool)
        cc = celdas[tramos[r]] + 1
        occ[cc[:, 0], cc[:, 1], cc[:, 2]] = True
        dil = np.zeros_like(occ)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    dil |= np.roll(np.roll(np.roll(occ, dx, 0), dy, 1), dz, 2)
        co = celdas[dueno != r] + 1
        res.append(int(dil[co[:, 0], co[:, 1], co[:, 2]].sum()))
    return max(res), int(np.mean(res))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("archivo")
    ap.add_argument("--procesos", type=int, nargs="+", default=[4, 8])
    ap.add_argument("--muestra", type=int, default=0)
    ap.add_argument("--salida", default="/tmp/demo_muestra.txt")
    a = ap.parse_args()

    datos = np.loadtxt(a.archivo, dtype=np.float64, ndmin=2)
    if datos.shape[1] != 7:
        sys.exit(f"se esperaban 7 columnas y hay {datos.shape[1]}")
    pos, vel, masa = datos[:, :3], datos[:, 3:6], datos[:, 6]
    n = len(datos)
    mn, arista = cubo(pos)
    print(f"N = {n}   masa total = {masa.sum():.4g}   "
          f"masas {'iguales' if np.ptp(masa) == 0 else 'distintas'}")
    print(f"cubo raiz: min = {np.round(mn, 3)}  arista = {arista:.4g}")

    # Concentracion: que fraccion del volumen contiene la mitad de la masa.
    celdas = np.clip(((pos - mn) / arista * G).astype(int), 0, G - 1)
    lin = (celdas[:, 0] * G + celdas[:, 1]) * G + celdas[:, 2]
    conteo = np.sort(np.bincount(lin, minlength=G ** 3))[::-1]
    acum = np.cumsum(conteo)
    vol50 = (np.searchsorted(acum, n * 0.5) + 1) / G ** 3 * 100
    vol90 = (np.searchsorted(acum, n * 0.9) + 1) / G ** 3 * 100
    ocupadas = (conteo > 0).sum() / G ** 3 * 100
    print(f"concentracion: 50 % de las particulas en el {vol50:.2f} % del volumen,"
          f" 90 % en el {vol90:.2f} %  (celdas ocupadas {ocupadas:.1f} %;"
          f" uniforme ~ 50 / 90 %)")
    print(f"celda mas densa: {conteo[0]} particulas"
          f" ({conteo[0] / max(n / G ** 3, 1e-9):.0f} veces el promedio)")

    # Velocidades: flujos que convergen (choque) o reposo.
    rapidez = np.linalg.norm(vel, axis=1)
    if rapidez.max() == 0:
        print("velocidades: todas en cero (colapso por gravedad)")
    else:
        rel = pos - pos.mean(0)
        conv = np.mean(np.sum(rel * vel, axis=1) < 0) * 100
        print(f"velocidades: |v| media {rapidez.mean():.3g}, max {rapidez.max():.3g};"
              f" {conv:.0f} % se acercan al centro"
              f"{'  -> flujos que CHOCAN' if conv > 70 else ''}")

    # Riesgo para MPI: particulas ajenas pegadas a cada tramo Morton.
    orden = orden_morton(pos, mn, arista)
    for P in a.procesos:
        mx, prom = vecinos_ajenos(celdas[orden], P)
        nivel = "ALTO" if mx > 0.15 * n / P else ("medio" if mx > 0.05 * n / P else "bajo")
        print(f"MPI con P={P}: particulas ajenas vecinas max {mx} (prom {prom})"
              f" = {mx / (n / P):.2f} veces un tramo -> LET {nivel}")

    if a.muestra:
        k = min(a.muestra, n)
        idx = np.sort(np.random.default_rng(7).choice(n, k, replace=False))
        m = datos[idx].copy()
        m[:, 6] *= masa.sum() / m[:, 6].sum()     # misma masa total
        np.savetxt(a.salida, m, fmt="%.17g")
        print(f"muestra de {k} particulas en {a.salida}")


if __name__ == "__main__":
    main()
