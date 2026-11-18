
from motor.fuerzas import (
    cargar_particulas,
    calcular_energia_total,
    ejecutar_simulacion
)

from motor.fuerzas_barnes_hut import (
    ejecutar_simulacion_barnes_hut
)

import math

# Parametros comunes
pasos = 100
dt = 0.001
theta = 0.3

# Cargar dos copias identicas de las particulas
particulas_directas = cargar_particulas("data/particles.txt")
particulas_bh = cargar_particulas("data/particles.txt")

# Energias iniciales
energia_inicial = calcular_energia_total(particulas_directas)

# Ejecutar motor secuencial
tiempo_directo = ejecutar_simulacion(
    particulas_directas, pasos, dt
)

# Ejecutar motor Barnes-Hut
tiempo_bh = ejecutar_simulacion_barnes_hut(
    particulas_bh, pasos, dt, theta
)

# Energias finales
energia_directa = calcular_energia_total(particulas_directas)
energia_bh = calcular_energia_total(particulas_bh)

# Diferencia maxima entre posiciones finales
diferencias = []

for i in range(len(particulas_directas)):
    diferencia = math.sqrt(
        sum(
            (
                particulas_directas[i]["posicion"][eje]
                - particulas_bh[i]["posicion"][eje]
            ) ** 2
            for eje in range(3)
        )
    )

    diferencias.append(diferencia)

print("COMPARACION DE MOTORES")
print("-----------------------------")

print("\nTIEMPOS")
print("Motor secuencial:", tiempo_directo)
print("Motor Barnes-Hut:", tiempo_bh)

print("\nENERGIA FINAL")
print("Secuencial:", energia_directa)
print("Barnes-Hut:", energia_bh)

print("\nDIFERENCIAS")
print("Diferencia maxima de posicion:", max(diferencias))
print("Diferencia promedio de posicion:",
      sum(diferencias) / len(diferencias))
