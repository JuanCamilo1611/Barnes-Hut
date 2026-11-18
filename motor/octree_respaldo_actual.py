

class NodoOctree:

    def __init__(self, centro, lado):
        """
        Representa una region cubica del espacio 3D.
        """

        # Centro geometrico del cubo
        self.centro = centro

        # Longitud del lado del cubo
        self.lado = lado

        # Masa total contenida en el nodo
        self.masa_total = 0.0

        # Centro de masa
        self.centro_masa = [0.0, 0.0, 0.0]

        # Ocho posibles nodos hijos
        self.hijos = [None] * 8

    # Indica si el nodo es una hoja
        self.es_hoja = True

        # Particula almacenada en este nodo
        self.particula = None

    def subdividir(self):
        """
        Divide el nodo actual en ocho cubos hijos.
        """

        self.es_hoja = False

        mitad = self.lado / 4
        nuevo_lado = self.lado / 2

        indice = 0

        for x in [-1, 1]:
            for y in [-1, 1]:
                for z in [-1, 1]:

                    centro_hijo = [
                        self.centro[0] + x * mitad,
                        self.centro[1] + y * mitad,
                        self.centro[2] + z * mitad
                    ]

                    self.hijos[indice] = NodoOctree(
                        centro_hijo,
                        nuevo_lado
                    )

                    indice += 1



    def obtener_indice_hijo(self, particula):
        """
        Determina en cual de los ocho octantes
        debe ubicarse una particula.
        """

        x, y, z = particula["posicion"]
        cx, cy, cz = self.centro

        indice = 0

        if x >= cx:
            indice += 4

        if y >= cy:
            indice += 2

        if z >= cz:
            indice += 1

        return indice

    def insertar(self, particula, profundidad=0, max_profundidad=20):
        """
        Inserta una particula dentro del Octree.
        """


        # Caso 1: nodo vacio y sin subdivisiones
        if self.es_hoja and self.particula is None:
            self.particula = particula
            return

        # Caso 2: nodo hoja que ya contiene una particula
        if self.es_hoja:

            if profundidad >= max_profundidad:
                raise ValueError(
                    "Se alcanzo la profundidad maxima del Octree"
                )

            # Recuperar la particula anterior
            particula_anterior = self.particula

            # Convertir el nodo en un nodo interno
            self.particula = None
            self.subdividir()

            # Reubicar la particula anterior
            indice_anterior = self.obtener_indice_hijo(
                particula_anterior
            )

            self.hijos[indice_anterior].insertar(
                particula_anterior,
                profundidad + 1,
                max_profundidad
            )

        # Caso 3: insertar la nueva particula en su hijo
        indice = self.obtener_indice_hijo(particula)

        self.hijos[indice].insertar(
            particula,
            profundidad + 1,
            max_profundidad
        )



    def calcular_aceleracion(
        self,
        objetivo,
        theta=0.5,
        G=1.0,
        epsilon=0.01
    ):
        """
        Calcula la aceleracion gravitacional sobre una particula
        utilizando el criterio de apertura Barnes-Hut.
        """

        # Nodo vacio: no aporta aceleracion.
        if self.masa_total == 0:
            return [0.0, 0.0, 0.0]

        posicion = objetivo["posicion"]

        # Vector desde la particula objetivo hasta el centro de masa.
        desplazamiento = [
            self.centro_masa[eje] - posicion[eje]
            for eje in range(3)
        ]

        distancia_cuadrada = sum(
            componente ** 2 for componente in desplazamiento
        )

        # Identificamos si el nodo contiene a la particula objetivo.
        contiene_objetivo = (
            all(
                abs(posicion[eje] - self.centro[eje])
                <= self.lado / 2
                for eje in range(3)
            )
        )

        # Caso 1: es una hoja.
        if self.es_hoja:

            if self.particula is None or self.particula is objetivo:
                return [0.0, 0.0, 0.0]

            masa = self.particula["masa"]

        else:
            # Caso 2: nodo interno.
            distancia = distancia_cuadrada ** 0.5

            # Si contiene al objetivo, debemos abrirlo.
            # Si no lo contiene, aplicamos el criterio theta.
            if contiene_objetivo or distancia == 0 or (
                self.lado / distancia >= theta
            ):
                aceleracion = [0.0, 0.0, 0.0]

                for hijo in self.hijos:
                    if hijo is not None:
                        contribucion = hijo.calcular_aceleracion(
                            objetivo, theta, G, epsilon
                        )

                        for eje in range(3):
                            aceleracion[eje] += contribucion[eje]

                return aceleracion

            masa = self.masa_total

        # Calculo de la aceleracion gravitacional.
        denominador = (
            distancia_cuadrada + epsilon ** 2
        ) ** 1.5

        factor = G * masa / denominador

        return [
            factor * desplazamiento[eje]
            for eje in range(3)
        ]



    def actualizar_masa_centro_masa(self):
        """
        Calcula la masa total y el centro de masa
        del nodo a partir de sus particulas o hijos.
        """

        # Caso 1: el nodo es una hoja con una particula
        if self.particula is not None:
            self.masa_total = self.particula["masa"]
            self.centro_masa = self.particula["posicion"].copy()
            return

        # Caso 2: el nodo es interno
        masa_total = 0.0
        suma_ponderada = [0.0, 0.0, 0.0]

        for hijo in self.hijos:
            if hijo is not None:

                hijo.actualizar_masa_centro_masa()

                masa = hijo.masa_total
                masa_total += masa

                for eje in range(3):
                    suma_ponderada[eje] += (
                        masa * hijo.centro_masa[eje]
                    )

        self.masa_total = masa_total

        # Evitar division entre cero
        if masa_total > 0:
            self.centro_masa = [
                suma_ponderada[eje] / masa_total
                for eje in range(3)
            ]
        else:
            self.centro_masa = [0.0, 0.0, 0.0]
