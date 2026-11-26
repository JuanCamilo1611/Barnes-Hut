#!/usr/bin/env python3
"""
Vecindad entre tramos Morton: explica el tamano del LET de MPI.

Reparte las particulas en P tramos de la curva Morton (como reparto.c con
pesos iguales) y cuenta, para cada proceso, cuantas particulas de OTROS
procesos caen en las celdas vecinas a las suyas (grilla 32^3). Ese conteo
sigue al LET recibido: estable en uniforme, crece con P en cumulos.

Uso: python3 tools/vecinos_corte.py <archivo> [<archivo> ...]
"""

import sys

import numpy as np

def sep(v):
    x = v.astype(np.uint64) & np.uint64(0x1fffff)
    for s,m in [(32,0x1f00000000ffff),(16,0x1f0000ff0000ff),(8,0x100f00f00f00f00f),(4,0x10c30c30c30c30c3),(2,0x1249249249249249)]:
        x = (x | (x << np.uint64(s))) & np.uint64(m)
    return x
G = 32  # grilla 32^3: celda = arista/32
for f in sys.argv[1:]:
    p = np.loadtxt(f, usecols=(0,1,2)); mn=p.min(0); ar=(p.max(0)-mn).max()
    q = np.clip(((p-mn)/ar*(1<<21)).astype(np.int64),0,(1<<21)-1)
    c = (sep(q[:,0])<<np.uint64(2))|(sep(q[:,1])<<np.uint64(1))|sep(q[:,2])
    ps = p[np.argsort(c,kind="stable")]
    cel = np.clip(((ps-mn)/ar*G).astype(int),0,G-1)
    for P in (2,4,8):
        idx = np.array_split(np.arange(len(ps)),P)
        dueno = np.empty(len(ps),int)
        for r,ix in enumerate(idx): dueno[ix]=r
        cerca=[]
        for r in range(P):
            occ = np.zeros((G+2,)*3,bool); cc=cel[idx[r]]+1; occ[cc[:,0],cc[:,1],cc[:,2]]=True
            dil = np.zeros_like(occ)
            for dx in (-1,0,1):
                for dy in (-1,0,1):
                    for dz in (-1,0,1):
                        dil |= np.roll(np.roll(np.roll(occ,dx,0),dy,1),dz,2)
            otros = dueno!=r; co = cel[otros]+1
            cerca.append(int(dil[co[:,0],co[:,1],co[:,2]].sum()))
        print(f"{f.split('/')[-1]:24s} P={P}: particulas ajenas vecinas  max={max(cerca):7d}  prom={int(np.mean(cerca)):7d}")
