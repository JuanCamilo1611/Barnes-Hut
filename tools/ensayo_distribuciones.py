#!/usr/bin/env python3
"""
Distribuciones NUEVAS para ensayar la demo (no son las del desarrollo).

  plummer    esfera de Plummer: nucleo muy denso y halo extendido
  colision   dos cumulos densos (sigma 0.03) que chocan descentrados
  filamento  particulas a lo largo de una curva fina (casi 1D)
  disco      disco delgado en rotacion con un bulbo central

Uso:
  python3 tools/ensayo_distribuciones.py plummer --N 100000 --out /tmp/x.txt
Formato: 7 columnas (x y z vx vy vz masa), masas iguales a 1/N.
"""

import argparse

import numpy as np


def plummer(n, rng):
    a = 0.05
    r = a / np.sqrt(rng.uniform(1e-6, 0.99, n) ** (-2 / 3) - 1)
    r = np.minimum(r, 1.0)
    d = rng.normal(size=(n, 3)); d /= np.linalg.norm(d, axis=1, keepdims=True)
    return 0.5 + r[:, None] * d * 0.5, np.zeros((n, 3))


def colision(n, rng):
    m = n // 2
    p1 = rng.normal([0.3, 0.45, 0.5], 0.03, (m, 3))
    p2 = rng.normal([0.7, 0.55, 0.5], 0.03, (n - m, 3))
    v = np.vstack([np.tile([0.8, 0, 0], (m, 1)), np.tile([-0.8, 0, 0], (n - m, 1))])
    return np.vstack([p1, p2]), v


def filamento(n, rng):
    t = rng.uniform(0, 1, n)
    centro = np.stack([t, 0.5 + 0.2 * np.sin(6 * t), 0.5 + 0.1 * np.cos(4 * t)], 1)
    return centro + rng.normal(0, 0.004, (n, 3)), np.zeros((n, 3))


def disco(n, rng):
    nb = n // 5
    r = np.sqrt(rng.uniform(0.0004, 0.16, n - nb))
    ang = rng.uniform(0, 2 * np.pi, n - nb)
    pd = np.stack([0.5 + r * np.cos(ang), 0.5 + r * np.sin(ang),
                   0.5 + rng.normal(0, 0.005, n - nb)], 1)
    vd = np.stack([-np.sin(ang), np.cos(ang), np.zeros_like(ang)], 1) * np.sqrt(1 / r)[:, None] * 0.3
    pb = rng.normal(0.5, 0.02, (nb, 3))
    return np.vstack([pd, pb]), np.vstack([vd, np.zeros((nb, 3))])


GEN = {"plummer": plummer, "colision": colision, "filamento": filamento, "disco": disco}

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("tipo", choices=GEN)
    ap.add_argument("--N", type=int, default=100000)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    pos, vel = GEN[a.tipo](a.N, rng)
    masa = np.full((a.N, 1), 1.0 / a.N)
    np.savetxt(a.out, np.hstack([pos, vel, masa]), fmt="%.10g")
    print(f"{a.tipo}: {a.N} particulas en {a.out}")
