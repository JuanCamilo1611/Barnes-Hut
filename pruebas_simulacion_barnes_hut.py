
from motor.fuerzas import cargar_particulas, calcular_energia_total
from motor.fuerzas_barnes_hut import ejecutar_simulacion_barnes_hut

particulas = cargar_particulas("data/particles.txt")

dt = 0.001
pasos = 100
theta = 0.3

energia_inicial = calcular_energia_total(particulas)

posicion_inicial = particulas[0]["posicion"].copy()

tiempo = ejecutar_simulacion_barnes_hut(
    particulas,
    pasos,
    dt,
    theta
)

energia_final = calcular_energia_total(particulas)

variacion = (
    abs(energia_final - energia_inicial)
    / abs(energia_inicial)
) * 100

print("SIMULACION BARNES-HUT")
print("-----------------------------")
print("Particulas:", len(particulas))
print("Pasos:", pasos)
print("Theta:", theta)
print("Tiempo simulado:", pasos * dt)
print("Tiempo de ejecucion:", tiempo, "segundos")

print("\nPRIMERA PARTICULA")
print("Posicion inicial:", posicion_inicial)
print("Posicion final:", particulas[0]["posicion"])
print("Velocidad final:", particulas[0]["velocidad"])

print("\nENERGIA")
print("Energia inicial:", energia_inicial)
print("Energia final:", energia_final)
print("Variacion porcentual:", variacion, "%")
