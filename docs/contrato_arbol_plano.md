# Contrato del árbol plano (Barnes-Hut)

Estructura de datos **compartida por las 4 versiones** (Python +
multiprocessing, C + OpenMP, MPI, CUDA). Fijarla antes de construir el
octree evita que cada versión invente un árbol distinto e incompatible.

Formato **SoA** (struct-of-arrays): un arreglo por campo, indexado por
nodo. **El nodo 0 es la raíz.**

> Estado: **borrador para ratificar (rev. 2).** Incorpora dos ajustes
> propuestos por Manuel en el PR y validados con medición: numeración en
> **preorden** (para que las ramas queden contiguas, §3) y
> `max_profundidad` = **20** con *bucketing* en el tope (§5). Pendiente la
> aprobación formal de Manuel para empezar el octree.

## 1. Arreglos del árbol (estos SÍ se comparten / serializan)

| Arreglo  | Forma      | Tipo    | Significado                                      |
|----------|------------|---------|--------------------------------------------------|
| `centro` | `[cap, 3]` | float64 | centro geométrico del cubo (x, y, z)             |
| `lado`   | `[cap]`    | float64 | arista del cubo                                  |
| `masa`   | `[cap]`    | float64 | masa total de la celda                           |
| `com`    | `[cap, 3]` | float64 | centro de masa de la celda                       |
| `hijos`  | `[cap, 8]` | int32   | índice de cada octante hijo; **-1 si no existe** |
| `part`   | `[cap]`    | int32   | índice de partícula si hoja ocupada; **-1 si interno/vacío** |

`cap` = capacidad reservada de nodos; `n_nodos` = nodos realmente usados.

**Auxiliares SOLO de construcción** (no se comparten, no se serializan, no
van al LET de MPI ni a la GPU):

| Arreglo         | Forma      | Tipo    | Significado                          |
|-----------------|------------|---------|--------------------------------------|
| `hojas_masa`    | `[cap]`    | float64 | masa acumulada en una hoja           |
| `hojas_momento` | `[cap, 3]` | float64 | masa·posición acumulada en una hoja  |

## 2. Convención de octantes (idéntica en referencia, Python, C y CUDA)

```
indice = 4·(x ≥ cx) + 2·(y ≥ cy) + 1·(z ≥ cz)
```

El bit de **x es el más significativo**. Esto **fija también el orden de
bits Morton** `(x, y, z)`, de modo que el ordenamiento Morton de C/CUDA
produzca el mismo octante que la inserción recursiva de Python.

## 3. Invariantes (obligatorias — habilitan el recorrido sin recursión)

1. **Nodo 0 = raíz.**
2. **Numeración en preorden (DFS):** cada nodo se numera antes que todos
   sus descendientes, y el subárbol de un nodo ocupa un **rango contiguo**
   `[nodo, último_descendiente]`. Esto garantiza a la vez dos cosas:
   - **todo hijo tiene índice mayor que su padre** → masa y centro de masa
     en **una sola pasada en orden inverso** (`n_nodos-1 → 0`), sin pila;
   - **una rama es un bloque contiguo** → el LET de MPI se empaqueta como
     un rango (§6).

   *Nota: la numeración por niveles (BFS) cumpliría "hijo > padre" pero NO
   daría subárboles contiguos; por eso se exige preorden. El orden Morton
   de las hojas coincide con su orden DFS, así que una construcción basada
   en Morton es compatible con preorden.*
3. `part[n] ≥ 0`  ⟺  hoja ocupada; entonces `hijos[n, *]` son todos `-1`.
4. Nodo interno  ⟺  `part[n] = -1`  y  algún `hijos[n, k] ≥ 0`.
5. `masa[n]` y `com[n]` solo son válidos **después** de la pasada de
   agregación.
6. Celda vacía (`masa[n] = 0`) **no aporta aceleración** → se salta en el
   recorrido.

## 4. Recorrido (criterio de apertura θ) — idéntico en todas las versiones

```
s = lado[n];   d = |com[n] − pos_objetivo|
si  s/d < θ   o  el nodo es hoja:   tratar la celda como una masa única (com[n], masa[n])
si no:                              visitar hijos[n, 0..7] con índice ≥ 0
```

Aceleración con suavizado (G = 1):

```
a += masa · (com − pos) / (d² + ε²)^(3/2)
```

Auto-interacción: saltar la hoja cuya `part == objetivo`.

## 5. Parámetros fijos del contrato

- `G = 1.0`.
- Tipos **float64 en CPU**. CUDA puede usar float32, documentando el error
  relativo vs la referencia.
- Centinela `-1` en `int32` para `hijos` y `part`.
- `max_profundidad` = **20**. *(Medido: las distribuciones uniforme,
  cúmulos y concentradas alcanzan profundidad ≤ 14 con N hasta 100k, y
  cap 20 y 40 producen árboles idénticos. Un cap bajo reduce además la
  pila de recorrido en CUDA.)*
- **Comportamiento en la profundidad máxima: *bucketing*, nunca error.**
  Si al llegar a `max_profundidad` la hoja ya está ocupada, se **acumulan**
  masa y momento en esa hoja (varias partículas por hoja) en vez de seguir
  subdividiendo o abortar. En el recorrido, una hoja-bucket se trata como
  una sola masa, y la auto-interacción solo excluye la partícula cuyo
  índice está en `part[n]`. *(Ya implementado así en `python_mp/arbol.py`;
  el núcleo C y CUDA deben replicarlo para no fallar con la distribución
  concentrada de la demo.)*
- La raíz cubre un cubo que contiene **todas** las partículas
  (bounding box + un pequeño margen).

## 6. Notas por versión (cómo cada quien consume el contrato)

### Python + multiprocessing (Camilo)
Los 6 arreglos del contrato viven en `shared_memory`; los auxiliares no.
Los workers crean un `ArbolPlano.con_arreglos(...)` que apunta a los
buffers compartidos, **sin copiar**. *(Ya implementado en
`python_mp/arbol.py`.)*

### Núcleo C (Manuel)
Los mismos 6 arreglos como punteros SoA dentro de un `struct ArbolPlano`:
`double *centro` (`3·cap`), `double *lado`, `double *masa`,
`double *com` (`3·cap`), `int32_t *hijos` (`8·cap`), `int32_t *part`.
La construcción desde orden Morton **debe respetar la invariante 2**
(padre < hijo) para que la agregación de COM sea una pasada lineal.

### MPI / LET (Camilo)
Gracias a la invariante 2, **una rama es un rango contiguo**
`[raíz_subárbol, último_descendiente]`. Para enviar una rama remota se
empaqueta ese rango con **reindexado local de `hijos`** (restar el offset
del subárbol). Fijar esto ahora evita rediseñar el árbol para cumplir
"cada proceso recibe solo las partes remotas del árbol que necesita".

### CUDA (Manuel)
La SoA ya es **coalescente**. Recorrido con **pila explícita** (sin
recursión); `hijos[n, 8]` con índices crecientes permite el walk
iterativo. La pila se dimensiona con `max_profundidad` (20): un cap bajo =
menos presión de registros/shared por hilo. Ordenar las partículas por
Morton antes de construir reduce la **divergencia de warp**.

## 7. Verificación

Toda versión que produzca aceleraciones se valida con `tools/verificar.py`
contra la suma directa O(N²):

- error relativo promedio ≤ 1 % (con N ≤ 2000), y
- deriva de la energía mecánica ≤ 1 %.

Una versión que no pase esta verificación se considera incorrecta.
