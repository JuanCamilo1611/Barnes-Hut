import math
from motor.fuerzas import cargar_particulas, calcular_todas_las_aceleraciones
from motor.fuerzas_barnes_hut import calcular_todas_las_aceleraciones_barnes_hut, construir_arbol

# Cargar los datos originales
particulas = cargar_particulas("data/particles.txt")

# Crear un desplazamiento espacial sin modificar el archivo original
desplazamiento = [2.0, -3.0, 4.0]

for p in particulas:
    for eje in range(3):
        p["posicion"][eje] += desplazamiento[eje]

# Calcular referencia con suma directa
directas = calcular_todas_las_aceleraciones(particulas)

# Construir árbol y comprobar masa
arbol = construir_arbol(particulas)

print("PARTICULAS:", len(particulas))
print("CENTRO DEL OCTREE:", arbol.centro)
print("LADO DEL OCTREE:", arbol.lado)
print("MASA TOTAL:", arbol.masa_total)

print("\nCOMPARACION DE ACELERACIONES")

for theta in [0.1, 0.3, 0.5]:
    aproximadas = calcular_todas_las_aceleraciones_barnes_hut(
        particulas, theta=theta
    )

    errores = []

    for directa, aproximada in zip(directas, aproximadas):
        diferencia = math.sqrt(sum(
            (aproximada[e] - directa[e]) ** 2 for e in range(3)
        ))

        magnitud = math.sqrt(sum(valor ** 2 for valor in directa))

        errores.append(diferencia / magnitud if magnitud > 0 else 0.0)

    print(f"\nTheta: {theta}")
    print(f"Error promedio: {sum(errores) / len(errores) * 100:.6f}%")
    print(f"Error maximo: {max(errores) * 100:.6f}%")
