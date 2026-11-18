
from motor.fuerzas import cargar_particulas, calcular_aceleracion
from motor.octree import NodoOctree
import math

particulas = cargar_particulas("data/particles.txt")

arbol = NodoOctree([0.5, 0.5, 0.5], 1.0)

for p in particulas:
    arbol.insertar(p)

arbol.actualizar_masa_centro_masa()

objetivo = particulas[0]

directa = calcular_aceleracion(objetivo, particulas)

magnitud_directa = math.sqrt(
    sum(x**2 for x in directa)
)

print("Aceleracion directa:", directa)
print("\nComparacion de theta")

for theta in [0.1, 0.3, 0.5, 0.8, 1.0]:
    aproximada = arbol.calcular_aceleracion(
        objetivo, theta=theta
    )

    error = math.sqrt(
        sum(
            (aproximada[i] - directa[i])**2
            for i in range(3)
        )
    )

    error_relativo = error / magnitud_directa * 100

    print(
        f"Theta: {theta:.1f} | "
        f"Error relativo: {error_relativo:.6f}%"
    )
