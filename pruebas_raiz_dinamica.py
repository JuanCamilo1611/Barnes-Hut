from motor.fuerzas_barnes_hut import construir_arbol

# Particulas de prueba con posiciones fuera del cubo original
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

# Construir el Octree
arbol = construir_arbol(particulas)

print("CENTRO DEL OCTREE:", arbol.centro)
print("LADO DEL CUBO:", arbol.lado)
print("MASA TOTAL:", arbol.masa_total)

print("\nVERIFICACION DE PARTICULAS:")

for i, particula in enumerate(particulas):
    posicion = particula["posicion"]

    dentro = all(
        arbol.centro[eje] - arbol.lado / 2 <= posicion[eje]
        <= arbol.centro[eje] + arbol.lado / 2
        for eje in range(3)
    )

    print(f"Particula {i + 1}: {posicion} -> {'DENTRO' if dentro else 'FUERA'}")
