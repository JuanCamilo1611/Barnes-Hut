# Simulación de N cuerpos con Barnes-Hut

Proyecto final de la asignatura de Concurrente: simulación gravitacional
de N partículas con el algoritmo de Barnes-Hut, paralelizada en cuatro
versiones:

| Versión | Carpeta | Estado |
|---|---|---|
| Referencia secuencial (Python) | `referencia/` | ✅ correcta |
| Python + multiprocessing | `python_mp/` | ⏳ pendiente |
| C con OpenMP | `c/openmp/` | ⏳ pendiente |
| MPI en C | `c/mpi/` | ⏳ pendiente |
| CUDA (Google Colab) | `cuda/` | ⏳ pendiente |

## Estructura

```
Barnes-Hut/
  README.md            este archivo
  data/                generadores y distribuciones
  referencia/          suma directa y octree de referencia (Python puro)
  tools/               verificar.py, medir.sh, graficas.py
  python_mp/           versión multiprocessing
  c/common/            núcleo C compartido: io, morton, arbol, recorrido
  c/openmp/            versión OpenMP
  c/mpi/               versión MPI
  cuda/                notebook de Colab
  informe/             borrador y figuras
  resultados/          CSV de mediciones (los grandes no se suben)
```

## Contrato común de entrada y salida

Todos los programas reciben los mismos argumentos posicionales:

```
programa <archivo> <N> <pasos> <dt> <theta> <eps> <hilos> [salida_acc]
```

| Argumento | Significado |
|---|---|
| `archivo` | partículas: N filas con 7 columnas `x y z vx vy vz masa` |
| `N` | número de partículas (debe coincidir con el archivo) |
| `pasos` | número de pasos de tiempo (100 para las mediciones) |
| `dt` | intervalo de tiempo, p. ej. `0.001` |
| `theta` | criterio de apertura, p. ej. `0.3` |
| `eps` | suavizado gravitacional, p. ej. `0.01` |
| `hilos` | procesos o hilos a usar (1 si no paraleliza) |
| `salida_acc` | archivo de aceleraciones, defecto `aceleraciones.txt` |

**Salida estándar** en líneas `clave=valor`:

```
tiempo_build=0.0123
tiempo_sort=0.0004
tiempo_walk=0.4500
tiempo_integrate=0.0100
tiempo_total=0.4727
energia_inicial=-1.234567
energia_final=-1.234500
```

`energia_*` solo es obligatorio con `N <= 5000` (la energía es O(N²)).

**Archivo de aceleraciones**: una línea por partícula, en el mismo orden
que el archivo de entrada, con tres valores `ax ay az`.

## Distribuciones

```bash
# uniforme en el cubo unitario, v = 0, masa = 1/N (desarrollo)
python3 data/generate_particles.py --distribucion uniforme --N 1000

# dos cúmulos en colisión (castiga el reparto de trabajo)
python3 data/generate_particles.py --distribucion cumulos --N 20000 --out data/cumulos.txt

# partículas concentradas en un volumen muy pequeño (castiga la profundidad)
python3 data/generate_particles.py --distribucion concentradas --N 200000
```

Siempre con `--seed` fijo para que las corridas sean reproducibles.
`data/ejemplo_100.txt` es la distribución pequeña versionada en el repo.

## Verificación

Puerta obligatoria antes de abrir un PR:

```bash
python3 tools/verificar.py --N 2000
```

Verificar una versión concreta:

```bash
python3 tools/verificar.py --particulas data/ejemplo_100.txt \
    --programa "python3 python_mp/sim.py" --hilos 4
```

Comparar dos archivos de aceleraciones:

```bash
python3 tools/verificar.py --comparar acc_a.txt acc_b.txt
```

Criterios: error relativo promedio ≤ 1 % contra la suma directa
(con N ≤ 2000) y derivación de la energía ≤ 1 %.

## Flujo de trabajo

- Una rama por tarea: `camilo/python-shm`, `manuel/nucleo-c`, `manuel/cuda-recorrido`.
- Commits pequeños y frecuentes; el historial de git puntúa en el informe.
- Antes de cada PR: `python3 tools/verificar.py --N 2000` en verde.
- Cada sesión empieza con `git pull origin main` y termina con `git push`.
- Sin copias `_backup`: git guarda el historial.
- Los PR se revisan en menos de 24 horas.
- Chat corto de 15 minutos cada noche: qué se hizo, qué sigue y bloqueos.
