
import time
from referencia.motor.octree import NodoOctree
from referencia.motor.fuerzas import calcular_energia_total

def construir_arbol(particulas):

    # Obtener los limites del espacio en cada eje
    minimo = [
        min(p["posicion"][eje] for p in particulas)
        for eje in range(3)
    ]

    maximo = [
        max(p["posicion"][eje] for p in particulas)
        for eje in range(3)
    ]

    # Calcular el centro del espacio ocupado
    centro = [
        (minimo[eje] + maximo[eje]) / 2
        for eje in range(3)
    ]

    # Calcular la mayor dimension para formar un cubo
    lado = max(
        maximo[eje] - minimo[eje]
        for eje in range(3)
    )

    # Agregar un pequeno margen de seguridad
    lado = max(lado, 1e-6) * 1.01

    # Crear la raiz adaptada al espacio ocupado
    arbol = NodoOctree(centro, lado)

    # Insertar todas las particulas
    for particula in particulas:
        arbol.insertar(particula)

    # Calcular masa total y centro de masa
    arbol.actualizar_masa_centro_masa()

    return arbol

def calcular_todas_las_aceleraciones_barnes_hut(
    particulas, theta=0.3, G=1.0, epsilon=0.01
):
    """
    Construye el Octree y calcula las aceleraciones
    de todas las particulas.
    """
    arbol = construir_arbol(particulas)

    aceleraciones = []

    for particula in particulas:
        aceleracion = arbol.calcular_aceleracion(
            particula,
            theta=theta,
            G=G,
            epsilon=epsilon
        )

        aceleraciones.append(aceleracion)

    return aceleraciones


def paso_leapfrog_barnes_hut(
    particulas, dt, theta=0.3, G=1.0, epsilon=0.01
):
    """
    Ejecuta un paso Leapfrog usando Barnes-Hut.
    """

    # Primera evaluacion de aceleraciones.
    aceleraciones_iniciales = (
        calcular_todas_las_aceleraciones_barnes_hut(
            particulas, theta, G, epsilon
        )
    )

    # Primera mitad de la actualizacion de velocidades.
    for i, particula in enumerate(particulas):
        for eje in range(3):
            particula["velocidad"][eje] += (
                0.5 * aceleraciones_iniciales[i][eje] * dt
            )

    # Actualizacion de posiciones.
    for particula in particulas:
        for eje in range(3):
            particula["posicion"][eje] += (
                particula["velocidad"][eje] * dt
            )

    # Segunda evaluacion: reconstruye el arbol con posiciones nuevas.
    aceleraciones_finales = (
        calcular_todas_las_aceleraciones_barnes_hut(
            particulas, theta, G, epsilon
        )
    )

    # Segunda mitad de la actualizacion de velocidades.
    for i, particula in enumerate(particulas):
        for eje in range(3):
            particula["velocidad"][eje] += (
                0.5 * aceleraciones_finales[i][eje] * dt
            )


def ejecutar_simulacion_barnes_hut(
    particulas, pasos, dt, theta=0.3
):
    """
    Ejecuta la simulacion y mide su tiempo.
    """
    inicio = time.perf_counter()

    for _ in range(pasos):
        paso_leapfrog_barnes_hut(
            particulas, dt, theta
        )

    fin = time.perf_counter()

    return fin - inicio
