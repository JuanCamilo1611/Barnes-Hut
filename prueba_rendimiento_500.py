
import random
import time
import math

from motor.fuerzas import calcular_todas_las_aceleraciones
from motor.fuerzas_barnes_hut import construir_arbol

# Configuracion de la prueba
N = 4000
SEMILLA = 42
THETA = 0.3

# Generar particulas reproducibles
generador = random.Random(SEMILLA)

particulas = []

for _ in range(N):
    particulas.append({
        "posicion": [
            generador.random(),
            generador.random(),
            generador.random()
        ],
        "velocidad": [0.0, 0.0, 0.0],
        "masa": 1.0 / N
    })

print(f"Numero de particulas: {N}")

# 1. Medir motor secuencial
inicio = time.perf_counter()

aceleraciones_directas = calcular_todas_las_aceleraciones(
    particulas
)

tiempo_directo = time.perf_counter() - inicio

# 2. Medir construccion del Octree
inicio = time.perf_counter()

arbol = construir_arbol(particulas)

tiempo_arbol = time.perf_counter() - inicio

# 3. Medir recorrido del mismo Octree
inicio = time.perf_counter()

aceleraciones_bh = []

for particula in particulas:
    aceleracion = arbol.calcular_aceleracion(
        particula,
        theta=THETA
    )

    aceleraciones_bh.append(aceleracion)

tiempo_recorrido = time.perf_counter() - inicio

# 4. Calcular tiempo total Barnes-Hut
tiempo_bh_total = tiempo_arbol + tiempo_recorrido

# 5. Calcular errores
errores = []

for directa, aproximada in zip(
    aceleraciones_directas,
    aceleraciones_bh
):
    diferencia = math.sqrt(sum(
        (aproximada[eje] - directa[eje]) ** 2
        for eje in range(3)
    ))

    magnitud = math.sqrt(sum(
        valor ** 2 for valor in directa
    ))

    error = diferencia / magnitud if magnitud > 0 else 0.0
    errores.append(error)

error_promedio = sum(errores) / len(errores)
error_maximo = max(errores)

# 6. Mostrar resultados
print("\nRESULTADOS")
print("------------------------------")
print(f"Masa total: {arbol.masa_total:.6f}")
print(f"Theta: {THETA}")
print(f"Tiempo secuencial: {tiempo_directo:.6f} segundos")
print(f"Tiempo construccion Octree: {tiempo_arbol:.6f} segundos")
print(f"Tiempo recorrido Octree: {tiempo_recorrido:.6f} segundos")
print(f"Tiempo total Barnes-Hut: {tiempo_bh_total:.6f} segundos")
print(f"Error promedio: {error_promedio * 100:.6f}%")
print(f"Error maximo: {error_maximo * 100:.6f}%")

if tiempo_bh_total > 0:
    relacion = tiempo_directo / tiempo_bh_total
    print(f"Relacion de tiempos: {relacion:.4f}x")

    if relacion > 1:
        print("Barnes-Hut fue mas rapido.")
    else:
        print("El motor secuencial fue mas rapido.")
