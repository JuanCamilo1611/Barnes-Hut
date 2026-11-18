import math

from motor.fuerzas import calcular_todas_las_aceleraciones
from motor.fuerzas_barnes_hut import calcular_todas_las_aceleraciones_barnes_hut


# Las mismas cuatro particulas de nuestra prueba anterior
particulas = [
    {
        "posicion": [-2.0, 1.0, 0.5],
        "velocidad": [0.0, 0.0, 0.0],
        "masa": 0.25
    },
    {
        "posicion": [3.0, 2.0, 1.5],
        "velocidad": [0.0, 0.0, 0.0],
        "masa": 0.25
    },
    {
        "posicion": [1.0, -3.0, 2.0],
        "velocidad": [0.0, 0.0, 0.0],
        "masa": 0.25
    },
    {
        "posicion": [2.0, 1.0, -4.0],
        "velocidad": [0.0, 0.0, 0.0],
        "masa": 0.25
    }
]

# Calcular aceleraciones con el motor secuencial
aceleraciones_directas = calcular_todas_las_aceleraciones(particulas)

print("VALIDACION DE ACELERACIONES")
print("---------------------------")

for theta in [0.1, 0.3, 0.5]:

    # Calcular aceleraciones con Barnes-Hut
    aceleraciones_bh = calcular_todas_las_aceleraciones_barnes_hut(
        particulas,
        theta=theta
    )

    errores = []

    for directa, aproximada in zip(
        aceleraciones_directas,
        aceleraciones_bh
    ):

        diferencia = [
            aproximada[eje] - directa[eje]
            for eje in range(3)
        ]

        magnitud_diferencia = math.sqrt(
            sum(valor ** 2 for valor in diferencia)
        )

        magnitud_directa = math.sqrt(
            sum(valor ** 2 for valor in directa)
        )

        error = magnitud_diferencia / magnitud_directa
        errores.append(error)

    error_promedio = sum(errores) / len(errores)
    error_maximo = max(errores)

    print(f"\nTheta: {theta}")
    print(f"Error promedio: {error_promedio * 100:.6f}%")
    print(f"Error maximo: {error_maximo * 100:.6f}%")
