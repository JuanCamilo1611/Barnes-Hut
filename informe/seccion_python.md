# Versión Python con multiprocessing

## Diseño

El octree es el **árbol plano** del contrato común
(`docs/contrato_arbol_plano.md`): arreglos NumPy paralelos (`centro`,
`lado`, `masa`, `com`, `hijos[8]`, `part`) numerados en **preorden**, con
profundidad máxima 20 y *bucketing* en el tope. Se construye igual que en
el núcleo C: códigos Morton de 63 bits, ordenamiento estable y creación de
los nodos nivel por nivel como grupos de partículas con el mismo prefijo
Morton; el preorden se obtiene ordenando los nodos por (código de inicio,
profundidad), y masa y centro de masa salen de sumas por rangos contiguos
(`np.add.reduceat`). No hay bucles de Python por partícula ni por nodo.

El recorrido se vectoriza por lotes de partículas y no usa recursión: cada
iteración evalúa toda la frontera (partícula, nodo) con NumPy, de modo que
el número de iteraciones es la profundidad del árbol. La integración es
*leapfrog* KDK con una sola evaluación de fuerzas por paso.

La dificultad de esta versión es **cómo compartir el árbol entre
procesos**. Se implementaron y midieron dos estrategias sobre
`multiprocessing.shared_memory`:

- **compartido**: el proceso principal construye el árbol una vez en
  bloques de memoria compartida y los trabajadores solo lo recorren;
- **local**: cada trabajador reconstruye su propio árbol a partir de las
  posiciones (también en memoria compartida).

## Correctitud

Con N = 2000, θ = 0.3, 50 pasos y 4 procesos, `tools/verificar.py` frente a
la suma directa O(N²):

| distribución | error promedio | error máximo | deriva de energía |
|---|---|---|---|
| uniforme | 0.12 % | 1.91 % | 0.00003 % |
| cúmulos | 0.10 % | 1.32 % | 0.0003 % |
| concentradas | 0.36 % | 0.63 % | 0.42 % |

Las tres pasan (error promedio y deriva < 1 %). Además, el árbol es
**idéntico al del núcleo C**: mismo número de nodos (296 965 con 200 000
partículas concentradas) y misma profundidad, y las aceleraciones coinciden
con las de `sim_barnes_hut` con diferencia 0.000000 % (N = 20 000). Las
estrategias secuencial, compartido y local producen resultados idénticos.

## Mediciones y decisión

Escalado con 1, 2, 4 y 8 procesos (mediana de 3 repeticiones, 3 pasos) con
N = 20 000 uniforme, en un Intel i5-8250U de **4 núcleos físicos / 8 hilos
lógicos** y 6 MiB de L3.

| estrategia | procesos | build (s) | walk (s) | total (s) | speedup | eficiencia |
|---|---|---|---|---|---|---|
| secuencial | 1 | 0.07 | 23.74 | 23.86 | 1.00 | 1.00 |
| compartido | 2 | 0.10 | 18.03 | 18.20 | 1.31 | 0.66 |
| compartido | 4 | 0.12 | 15.84 | 16.01 | 1.49 | 0.37 |
| compartido | 8 | 0.14 | 15.42 | 15.66 | 1.52 | 0.19 |
| local | 2 | 0.09 | 15.96 | 16.12 | 1.48 | 0.74 |
| local | 4 | 0.20 | 15.37 | 15.68 | 1.52 | 0.38 |
| local | 8 | 0.41 | 15.00 | 15.45 | 1.54 | 0.19 |

![Speedup y eficiencia](../resultados/plots/escalado_python_mp.png)

![Tiempo de construcción y total](../resultados/plots/fases_python_mp.png)

**La decisión se midió dos veces.** Con la primera versión del árbol
(inserción partícula por partícula, 2.47 s de *build*), `local` perdía:
reconstruir en cada trabajador llevaba el *build* a 11.2 s con 8 procesos y
el *speedup* a 0.94×. Al alinear el árbol con el contrato y construirlo por
Morton, el *build* bajó 34 veces (2.47 → 0.07 s) y la diferencia en el
tiempo total desapareció: hoy ambas estrategias están a menos de un 2 % con
4 y 8 procesos, dentro de la variación entre corridas de este portátil.

**Se mantiene `compartido`** por tres razones que no dependen del ruido:

1. **Trabajo redundante que crece con P.** El *build* de `local` crece
   linealmente con los procesos (0.09 → 0.20 → 0.41 s), porque son P copias
   del mismo trabajo; el de `compartido` se hace una sola vez.
2. **Memoria.** Cada árbol ocupa unos 100 B por nodo (≈ 30 MB con
   N = 200 000); `local` guarda P copias (≈ 240 MB con 8 procesos) y
   compite por la misma caché L3 que ya satura el recorrido (siguiente
   sección).
3. **Coherencia con MPI.** Construir una vez y recorrer en paralelo es el
   esquema que se extiende a MPI, donde cada proceso recibe solo las ramas
   que necesita (LET), que el preorden deja contiguas.

## Por qué el escalado se aplana

El *speedup* se estanca en ≈1.5× desde 4 procesos. **No es la fracción
secuencial:** *build* e integración son el 0.3 % del tiempo secuencial, y
Amdahl permitiría 3.96× con 4 procesos y 7.8× con 8.

La causa aparece al medir cada trabajador por separado (un cálculo de
fuerzas, estrategia `compartido`):

| procesos | partículas por proceso | tiempo por proceso | partículas/s en total |
|---|---|---|---|
| 1 | 20 000 | 6.06 s | 3 300 |
| 2 | 10 000 | 4.1 s | 4 870 |
| 4 | 5 000 | 4.1 s | 4 880 |
| 8 | 2 500 | 4.0 s | 4 960 |

El reparto está balanceado (todos los trabajadores tardan lo mismo) y el
tiempo de CPU es igual al de pared, así que no hay esperas: **el
rendimiento total se satura en ≈4 900 partículas/s desde 2 procesos**. Los
procesos compiten por un recurso compartido:

1. **Ancho de banda de memoria y L3.** Con N = 20 000 el árbol ocupa solo
   ≈ 3 MB, pero el recorrido vectorizado con NumPy crea en cada nivel de la
   frontera arreglos temporales de cientos de miles de pares (partícula,
   nodo) que se escriben y se leen una vez: es tráfico de memoria, no
   cómputo. Los 4 núcleos comparten la L3 de 6 MiB y el controlador de
   memoria; con N = 200 000 el árbol (≈ 30 MB) tampoco cabe en la L3.
   Reducir el lote para que los temporales quepan en caché y ordenar las
   partículas por Morton mejora el tiempo con 4 procesos solo de forma
   marginal (de 3.9 a 3.0 s), lo que confirma que el límite es el sistema
   de memoria y no el cómputo.
2. **Frecuencia del procesador.** El i5-8250U es un procesador de portátil:
   con un solo núcleo activo llega a 3.4 GHz y baja con todos los núcleos
   cargados, así que cada proceso es más lento cuando hay más procesos.
3. **Hilos lógicos.** De 4 a 8 procesos no hay ganancia (1.49 → 1.52×): los
   8 hilos son *hyperthreading* sobre los mismos 4 núcleos y comparten los
   mismos puertos de memoria, ya saturados.
