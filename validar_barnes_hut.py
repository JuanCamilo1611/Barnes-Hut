
from motor.fuerzas import cargar_particulas, calcular_aceleracion
from motor.octree import NodoOctree
import math

particulas = cargar_particulas("data/particles.txt")

arbol = NodoOctree([0.5, 0.5, 0.5], 1.0)

for particula in particulas:
    arbol.insertar(particula)

arbol.actualizar_masa_centro_masa()

thetas = [0.1, 0.3, 0.5, 0.8, 1.0]

print("VALIDACION GLOBAL BARNES-HUT")
print("Cantidad de particulas:", len(particulas))

for theta in thetas:

    errores = []

    for particula in particulas:

        directa = calcular_aceleracion(
            particula, particulas
        )

        aproximada = arbol.calcular_aceleracion(
            particula, theta=theta
        )

        diferencia = math.sqrt(
            sum(
                (aproximada[eje] - directa[eje]) ** 2
                for eje in range(3)
            )
        )

        magnitud_directa = math.sqrt(
            sum(componente ** 2 for componente in directa)
        )

        if magnitud_directa > 1e-12:
            error = diferencia / magnitud_directa * 100
            errores.append(error)

    promedio = sum(errores) / len(errores)
    maximo = max(errores)

    print(f"\nTheta: {theta}")
    print(f"Error promedio: {promedio:.6f}%")
    print(f"Error maximo: {maximo:.6f}%")
