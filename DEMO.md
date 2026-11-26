# Guía de la demostración en vivo

Sin internet ni asistentes de IA. Todo se hace desde la raíz del
repositorio en WSL.

## Antes de empezar (una sola vez)

```bash
sudo apt install htop          # monitor que pide el enunciado
git pull && bash demo.sh data/ejemplo_100.txt -m -P "1 2"   # comprobar que compila y corre
```

Durante la demo hay que tener **dos terminales**: una con `htop` y otra
para los comandos. En `htop`, la tecla F2 permite mostrar los 8 hilos
lógicos por separado.

## 1. Correr la distribución desconocida (≈ 1 min con 100k)

```bash
bash demo.sh <archivo>                 # compila, diagnostica, verifica y mide
bash demo.sh <archivo> -r 3            # si hay tiempo: el mejor de 3 corridas
bash demo.sh <archivo> -y -P "1 4"     # incluir Python (lento con N grande)
```

CUDA se corre en Colab: abrir `cuda/barnes_hut_cuda.ipynb`, clonar
**main** y ejecutar
`./cuda/sim_cuda <archivo> <N> 1 0.001 0.5 0.01 128`.

## 2. Cómo leer el diagnóstico y qué decir

| Lo que muestra | Qué significa | Qué predecir |
|---|---|---|
| **50 % en < 1 % del volumen** | Muy concentrada: el árbol es profundo en el núcleo | El recorrido abre más nodos; el walk sube frente a la uniforme |
| **flujos que CHOCAN** | Cúmulos que se acercan | La energía puede derivar en pocos pasos por encuentros cercanos: bajar dt |
| **LET ALTO** | Los cortes Morton caen dentro de zonas densas | MPI escala peor con 8 que con 4, y OpenMP lo supera |
| **LET bajo** | Los cortes caen entre zonas poco densas | MPI ≈ OpenMP; MPI incluso un poco mejor |
| todas las velocidades en 0 | Colapso gravitacional | Nada especial con 1 paso |

En `htop`, señales a comentar:

- **MPI con LET alto:** algunas barras al 100 % y otras esperando. Eso
  es desbalance y colectivas esperando al proceso más lento.
- **OpenMP con `dynamic`:** todas las barras parejas.
- **Con 8:** los 8 hilos lógicos ocupados, pero el tiempo casi no baja.
  Son 4 núcleos físicos y el límite es la memoria.

Ejemplos de ensayo, con `tools/ensayo_distribuciones.py`:

| Distribución | Diagnóstico | Resultado |
|---|---|---|
| **plummer**, 100k | 50 % en el 0.08 % del volumen, LET ALTO | OpenMP 8: 3.93×; MPI 8: 1.45× (peor que con 4: 2.33×) |
| **colisión**, 50k | flujos que CHOCAN, LET ALTO | MPI con bisección y 4 procesos: 25 % del tiempo en comunicación, desbalance 1.32 |

## 3. Modificaciones típicas: dónde se cambia cada cosa

| Pedido | Cómo | Dónde |
|---|---|---|
| **Cambiar θ** | Es un argumento, no hace falta recompilar: `bash demo.sh <archivo> -t 0.7` | — |
| **Criterio de apertura** | `abrir = (s*s >= theta2*d2)` | `c/common/recorrido.c:37` (C, OpenMP, MPI local) |
| | El mismo criterio en el LET | `c/mpi/let.c:57` |
| | El mismo criterio en la GPU | `cuda/sim_cuda.cu:127` |
| **Particionado de MPI** | Sin recompilar: `-b` (bisección en x) o `BH_PARTICION=biseccion` | — |
| | Estrategia nueva: agregar un `case` en `clave_particula()` (por ejemplo, bisección en y: `pos[3*i+1]`) | `c/mpi/reparto.c:48-65` |
| **Balance de MPI** | `-c` (por número de partículas) o `BH_BALANCE=particulas` | `c/mpi/sim_mpi.c:284` |
| **Reparto de OpenMP** | Sin recompilar: `OMP_SCHEDULE=static` (o `guided`, `dynamic,16`) | La directiva `schedule(runtime)` está en `c/openmp/sim_openmp.c:129` |
| **Hilos por bloque (CUDA)** | Último argumento de `sim_cuda` (32, 64, 128, 256) | Valor por defecto en `cuda/sim_cuda.cu:41` |
| **Asignación de partículas a warps** | El hilo `j` toma la partícula `orden[j]`, en orden Morton. Con `int i = j;` se usa el orden del archivo: más divergencia y peor caché | `cuda/sim_cuda.cu:102` |
| **Profundidad máxima del árbol** | `MAX_PROFUNDIDAD` | `c/common/arbol.h:18` (recompilar: `make -C c/common`, etc.) |
| **Estrategia de Python** | `BH_ESTRATEGIA=local` o `compartido` | `python_mp/paralelo.py:12` |

Después de cambiar C, `demo.sh` recompila solo. Para comprobar que el
cambio no rompió nada sin medir, se usa `bash demo.sh <archivo> -v`.

## 4. Si algo falla

| Síntoma | Causa probable y salida |
|---|---|
| Solo falla la **energía** | dt grande con encuentros cercanos: `-d 0.00025`. La suma directa deriva igual |
| **Error > 1 %** con θ > 0.5 | Es lo esperado. Lo que importa es que todas las versiones den el mismo error que el C secuencial |
| MPI con P = 8 dice "not enough slots" | Falta `--use-hwthread-cpus`; `demo.sh` ya lo pone |
| MPI con P > 8 | Usa `--oversubscribe`; esperar que empeore. Es el efecto de usar más procesos que núcleos |
| Tiempos raros entre corridas | Ruido térmico de hasta 25 %: repetir con `-r 3` |
