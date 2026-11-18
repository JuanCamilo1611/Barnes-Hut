
import time


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



def calcular_energia_cinetica(particulas):
    """
    Calcula la energia cinetica total del sistema.
    """
    energia = 0.0

    for particula in particulas:
        masa = particula["masa"]
        vx, vy, vz = particula["velocidad"]

        rapidez_cuadrada = vx**2 + vy**2 + vz**2

        energia += 0.5 * masa * rapidez_cuadrada

    return energia


def calcular_energia_potencial(particulas, G=1.0, epsilon=0.01):
    """
    Calcula la energia potencial gravitacional total.
    Cada pareja se cuenta una sola vez.
    """
    energia = 0.0

    for i in range(len(particulas)):
        pi = particulas[i]
        xi, yi, zi = pi["posicion"]

        for j in range(i + 1, len(particulas)):
            pj = particulas[j]
            xj, yj, zj = pj["posicion"]

            dx = xj - xi
            dy = yj - yi
            dz = zj - zi

            distancia_cuadrada = dx**2 + dy**2 + dz**2
            distancia_suavizada = distancia_cuadrada + epsilon**2

            energia -= (
                G * pi["masa"] * pj["masa"]
                / (distancia_suavizada**0.5)
            )

    return energia


def calcular_energia_total(particulas):
    """
    Calcula la energia mecanica total del sistema.
    """
    cinetica = calcular_energia_cinetica(particulas)
    potencial = calcular_energia_potencial(particulas)

    return cinetica + potencial



def ejecutar_simulacion(particulas, pasos, dt):

    tiempo_inicio = time.perf_counter()

    for paso in range(pasos):
        paso_leapfrog(particulas, dt)

    tiempo_fin = time.perf_counter()

    tiempo_total = tiempo_fin - tiempo_inicio

    return tiempo_total



if __name__ == "__main__":
    particulas = cargar_particulas("data/particles.txt")

    print("Particulas cargadas:", len(particulas))

    posicion_inicial = particulas[0]["posicion"].copy()

    # Configuracion de la simulacion
    pasos = 100
    dt = 0.001

    # Energia inicial
    energia_inicial = calcular_energia_total(particulas)

    # Ejecutar simulacion
    tiempo_total = ejecutar_simulacion(particulas, pasos, dt)

    # Energia final
    energia_final = calcular_energia_total(particulas)

    # Calcular variacion porcentual
    variacion = (
        abs(energia_final - energia_inicial)
        / abs(energia_inicial)
    ) * 100

    print("\nRESULTADOS DE LA SIMULACION")
    print("Pasos ejecutados:", pasos)
    print("Intervalo de tiempo:", dt)
    print("Tiempo simulado:", pasos * dt)

    print("\nPRIMERA PARTICULA")
    print("Posicion inicial:", posicion_inicial)
    print("Posicion final:", particulas[0]["posicion"])
    print("Velocidad final:", particulas[0]["velocidad"])

    print("\nENERGIA DEL SISTEMA")
    print("Energia inicial:", energia_inicial)
    print("Energia final:", energia_final)
    print("Variacion porcentual:", variacion, "%")

    print("\nRENDIMIENTO")
    print("Tiempo de ejecucion:", tiempo_total, "segundos")





