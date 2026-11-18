
import random
import os

# ==========================================
# CONFIGURACION DE LA SIMULACION
# ==========================================

N = 100                  # Cantidad de particulas
SEED = 42                # Semilla para reproducibilidad
OUTPUT = "data/particles.txt"

# ==========================================
# GENERACION DE PARTICULAS
# ==========================================

def generate_particles(n, seed):

    # Inicializar generador aleatorio
    rng = random.Random(seed)

    particles = []

    # Masa igual para todas las particulas
    mass = 1.0 / n

    for i in range(n):

        # Posiciones aleatorias dentro del cubo unitario
        x = rng.uniform(0.0, 1.0)
        y = rng.uniform(0.0, 1.0)
        z = rng.uniform(0.0, 1.0)

        # Velocidades iniciales iguales a cero
        vx = 0.0
        vy = 0.0
        vz = 0.0

        # Guardar los datos de la particula
        particle = [x, y, z, vx, vy, vz, mass]

        particles.append(particle)

    return particles


# ==========================================
# GUARDAR PARTICULAS EN UN ARCHIVO
# ==========================================

def save_particles(particles, filename):

    # Crear directorio si no existe
    os.makedirs(os.path.dirname(filename), exist_ok=True)

    with open(filename, "w") as file:

        for particle in particles:

            line = " ".join(f"{value:.8f}" for value in particle)

            file.write(line + "\n")


# ==========================================
# PROGRAMA PRINCIPAL
# ==========================================

def main():

    print("================================")
    print(" GENERADOR DE PARTICULAS")
    print("================================")

    particles = generate_particles(N, SEED)

    save_particles(particles, OUTPUT)

    print(f"Particulas generadas: {N}")
    print(f"Semilla utilizada: {SEED}")
    print(f"Masa por particula: {1.0/N}")
    print(f"Archivo generado: {OUTPUT}")

    print("================================")
    print("Generacion completada.")
    print("================================")


if __name__ == "__main__":
    main()
