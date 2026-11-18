def cargar_particulas(ruta):
    """
    Carga las partículas desde un archivo de texto.

    Cada fila contiene:
    x, y, z, vx, vy, vz, masa
    """

    particulas = []

    with open(ruta, "r") as archivo:
        for linea in archivo:
            valores = list(map(float, linea.split()))

            particula = {
                "posicion": valores[0:3],
                "velocidad": valores[3:6],
                "masa": valores[6]
            }

            particulas.append(particula)

    return particulas


"""Calcular la gravedad"""

def calcular_aceleracion(particula, particulas, G=1.0, epsilon=0.01):
    """
    Calcula la aceleración gravitacional sobre una partícula
    mediante suma directa.
    """

    ax = 0.0
    ay = 0.0
    az = 0.0

    xi, yi, zi = particula["posicion"]

    for otra in particulas:

        if otra is particula:
            continue

        xj, yj, zj = otra["posicion"]
        masa = otra["masa"]

        dx = xj - xi
        dy = yj - yi
        dz = zj - zi

        distancia_cuadrada = dx**2 + dy**2 + dz**2

        distancia_suavizada = (
            distancia_cuadrada + epsilon**2
        )

        factor = G * masa / (distancia_suavizada**1.5)

        ax += factor * dx
        ay += factor * dy
        az += factor * dz

    return [ax, ay, az]


def calcular_todas_las_aceleraciones(particulas):

    aceleraciones = []

    for particula in particulas:

        aceleracion = calcular_aceleracion(
            particula,
            particulas
        )

        aceleraciones.append(aceleracion)

    return aceleraciones


def paso_leapfrog(particulas, dt):

    # 1. Calcular aceleraciones actuales
    aceleraciones = calcular_todas_las_aceleraciones(particulas)

    velocidades_medias = []

    # 2. Calcular velocidades a medio paso
    # 3. Actualizar posiciones
    for i, particula in enumerate(particulas):

        velocidad = particula["velocidad"]
        aceleracion = aceleraciones[i]

        velocidad_media = [
            velocidad[e] + 0.5 * aceleracion[e] * dt
            for e in range(3)
        ]

        velocidades_medias.append(velocidad_media)

        particula["posicion"] = [
            particula["posicion"][e] + velocidad_media[e] * dt
            for e in range(3)
        ]

    # 4. Calcular nuevas aceleraciones
    nuevas_aceleraciones = calcular_todas_las_aceleraciones(
        particulas
    )

    # 5. Completar actualización de velocidades
    for i, particula in enumerate(particulas):

        particula["velocidad"] = [
            velocidades_medias[i][e]
            + 0.5 * nuevas_aceleraciones[i][e] * dt
            for e in range(3)
        ]



if __name__ == "__main__":

    particulas = cargar_particulas("data/particles.txt")

    print("Particulas cargadas:", len(particulas))

    # Guardar posicion inicial de la primera particula
    posicion_inicial = particulas[0]["posicion"].copy()

    # Ejecutar un paso de simulacion
    dt = 0.001

    paso_leapfrog(particulas, dt)

    print("\nSIMULACION LEAPFROG")
    print("Intervalo de tiempo:", dt)

    print("\nPrimera particula:")
    print("Posicion inicial:", posicion_inicial)
    print("Posicion final:", particulas[0]["posicion"])
    print("Velocidad final:", particulas[0]["velocidad"])
