"""
Octree plano (arreglos) para Barnes-Hut.

Sustituye al arbol de objetos de referencia/motor/octree.py por arreglos
planos, que es lo que exige el enunciado para la version multiprocessing
(memoria compartida) y la misma estructura que usaran las versiones C y
CUDA.

Cada nodo n guarda:

    centro[n, 3]        centro geometrico del cubo
    lado[n]             arista del cubo
    masa[n]             masa total de la celda
    com[n, 3]           centro de masa de la celda
    hojas_masa[n]       masa acumulada en una hoja
    hojas_momento[n,3]  masa * posicion acumulada en una hoja
    hijos[n, 8]         indice de los ocho hijos, -1 si no existen
    part[n]             particula de la hoja, -1 si el nodo es interno

Los hijos siempre se crean despues que el padre, por lo que la masa y el
centro de masa se calculan recorriendo los nodos en orden inverso de
creacion, sin recursion.
"""

import numpy as np


class ArbolPlano:

    def __init__(self, capacidad=1024, max_profundidad=40):
        self.max_profundidad = max_profundidad
        self.capacidad = 0
        self.n_nodos = 0
        self.compartido = False
        self.centro = None
        self.lado = None
        self.masa = None
        self.com = None
        self.hojas_masa = None
        self.hojas_momento = None
        self.hijos = None
        self.part = None
        self._reservar(capacidad)

    @classmethod
    def con_arreglos(cls, centro, lado, masa, com, hijos, part,
                     n_nodos, capacidad, compartido=False,
                     construible=False, max_profundidad=40):
        """
        Arbol que apunta a arreglos ya existentes.

        Es lo que usan los trabajadores para recorrer el arbol que
        construyo el proceso principal en memoria compartida, sin
        copiar nada.  `construible` reserva ademas los arreglos
        auxiliares de la construccion, que no se comparten.
        """

        arbol = cls.__new__(cls)
        arbol.max_profundidad = max_profundidad
        arbol.capacidad = capacidad
        arbol.n_nodos = n_nodos
        arbol.compartido = compartido
        arbol.centro = centro
        arbol.lado = lado
        arbol.masa = masa
        arbol.com = com
        arbol.hijos = hijos
        arbol.part = part

        if construible:
            arbol.hojas_masa = np.zeros(capacidad, dtype=np.float64)
            arbol.hojas_momento = np.zeros(
                (capacidad, 3), dtype=np.float64
            )
        else:
            # Solo se necesitan durante la construccion.
            arbol.hojas_masa = None
            arbol.hojas_momento = None

        return arbol

    # ------------------------------------------------------------
    # almacenamiento
    # ------------------------------------------------------------

    def _reservar(self, necesarios):
        """Amplia los arreglos si no cabe el numero pedido de nodos."""

        if necesarios <= self.capacidad:
            return

        if self.compartido:
            # Los arreglos viven en memoria compartida: no se pueden
            # redimensionar sin invalidar a los demas procesos.
            raise RuntimeError(
                "capacidad compartida agotada: %d nodos disponibles, "
                "se necesitaban %d" % (self.capacidad, necesarios)
            )

        nueva = max(self.capacidad * 2, necesarios)

        if self.capacidad == 0:
            self.centro = np.zeros((nueva, 3), dtype=np.float64)
            self.lado = np.zeros(nueva, dtype=np.float64)
            self.masa = np.zeros(nueva, dtype=np.float64)
            self.com = np.zeros((nueva, 3), dtype=np.float64)
            self.hojas_masa = np.zeros(nueva, dtype=np.float64)
            self.hojas_momento = np.zeros((nueva, 3), dtype=np.float64)
            self.hijos = np.full((nueva, 8), -1, dtype=np.int32)
            self.part = np.full(nueva, -1, dtype=np.int32)
        else:
            self.centro = np.resize(self.centro, (nueva, 3))
            self.lado = np.resize(self.lado, nueva)
            self.masa = np.resize(self.masa, nueva)
            self.com = np.resize(self.com, (nueva, 3))
            self.hojas_masa = np.resize(self.hojas_masa, nueva)
            self.hojas_momento = np.resize(self.hojas_momento, (nueva, 3))
            self.hijos = np.resize(self.hijos, (nueva, 8))
            self.part = np.resize(self.part, nueva)

        self.capacidad = nueva

    def _nuevo_nodo(self, centro, lado):
        indice = self.n_nodos
        self._reservar(indice + 1)
        self.n_nodos = indice + 1

        self.centro[indice] = centro
        self.lado[indice] = lado
        self.masa[indice] = 0.0
        self.com[indice] = 0.0
        self.hojas_masa[indice] = 0.0
        self.hojas_momento[indice] = 0.0
        self.hijos[indice] = -1
        self.part[indice] = -1

        return indice

    # ------------------------------------------------------------
    # construccion
    # ------------------------------------------------------------

    @staticmethod
    def _octante(centro, posicion):
        indice = 0
        if posicion[0] >= centro[0]:
            indice += 4
        if posicion[1] >= centro[1]:
            indice += 2
        if posicion[2] >= centro[2]:
            indice += 1
        return indice

    def _subdividir(self, nodo):
        centro = self.centro[nodo]
        mitad = self.lado[nodo] / 4.0
        nuevo_lado = self.lado[nodo] / 2.0

        creados = np.empty(8, dtype=np.int32)
        posicion = 0

        for sx in (-1.0, 1.0):
            for sy in (-1.0, 1.0):
                for sz in (-1.0, 1.0):
                    hijo_centro = (
                        centro[0] + sx * mitad,
                        centro[1] + sy * mitad,
                        centro[2] + sz * mitad,
                    )
                    creados[posicion] = self._nuevo_nodo(
                        hijo_centro, nuevo_lado
                    )
                    posicion += 1

        self.hijos[nodo] = creados

    def _insertar(self, indice, posiciones, masas):
        nodo = 0
        profundidad = 0

        while True:
            if self.hijos[nodo, 0] == -1:
                # Nodo hoja.
                anterior = self.part[nodo]

                if anterior == -1:
                    self.part[nodo] = indice
                    self.hojas_masa[nodo] = masas[indice]
                    self.hojas_momento[nodo] = (
                        masas[indice] * posiciones[indice]
                    )
                    return

                if profundidad >= self.max_profundidad:
                    # Se alcanzo la profundidad maxima: se acumulan.
                    self.hojas_masa[nodo] += masas[indice]
                    self.hojas_momento[nodo] += (
                        masas[indice] * posiciones[indice]
                    )
                    return

                # Convertir la hoja en nodo interno.
                self._subdividir(nodo)

                hijo_anterior = self.hijos[
                    nodo, self._octante(self.centro[nodo],
                                        posiciones[anterior])
                ]
                self.part[hijo_anterior] = anterior
                self.hojas_masa[hijo_anterior] = masas[anterior]
                self.hojas_momento[hijo_anterior] = (
                    masas[anterior] * posiciones[anterior]
                )

                self.part[nodo] = -1
                self.hojas_masa[nodo] = 0.0
                self.hojas_momento[nodo] = 0.0

                nodo = self.hijos[
                    nodo, self._octante(self.centro[nodo],
                                        posiciones[indice])
                ]
                profundidad += 1
                continue

            nodo = self.hijos[
                nodo, self._octante(self.centro[nodo],
                                    posiciones[indice])
            ]
            profundidad += 1

    def _calcular_masa_com(self):
        """
        Masa y centro de masa de cada nodo.

        Los hijos se crean despues que el padre, así que recorrer los
        nodos al reves garantiza que los hijos ya esten resueltos.
        """
        for nodo in range(self.n_nodos - 1, -1, -1):

            if self.hijos[nodo, 0] == -1:
                masa = self.hojas_masa[nodo]
                self.masa[nodo] = masa

                if masa > 0.0:
                    self.com[nodo] = self.hojas_momento[nodo] / masa
                else:
                    self.com[nodo] = 0.0
                continue

            masa_total = 0.0
            momento = np.zeros(3, dtype=np.float64)
            hijos = self.hijos[nodo]

            for k in range(8):
                hijo = hijos[k]
                if hijo == -1:
                    continue

                masa_hijo = self.masa[hijo]
                if masa_hijo == 0.0:
                    continue

                masa_total += masa_hijo
                momento += masa_hijo * self.com[hijo]

            self.masa[nodo] = masa_total

            if masa_total > 0.0:
                self.com[nodo] = momento / masa_total
            else:
                self.com[nodo] = 0.0

    def construir(self, posiciones, masas):
        """Construye el arbol completo a partir de las posiciones."""

        posiciones = np.asarray(posiciones, dtype=np.float64)
        masas = np.asarray(masas, dtype=np.float64)

        minimo = posiciones.min(axis=0)
        maximo = posiciones.max(axis=0)

        centro = (minimo + maximo) * 0.5
        lado = float((maximo - minimo).max())
        lado = max(lado, 1e-6) * 1.01

        self.n_nodos = 0
        self._nuevo_nodo(centro, lado)

        for indice in range(len(posiciones)):
            self._insertar(indice, posiciones, masas)

        self._calcular_masa_com()

        return self

    # ------------------------------------------------------------
    # recorrido
    # ------------------------------------------------------------

    def aceleraciones(self, posiciones, theta=0.3, G=1.0, epsilon=0.01,
                      lote=2048, desplazamiento=0):
        """
        Aceleraciones sobre todas las particulas, sin recursion y sin
        bucles de Python por nodo.

        Se procesan las particulas en lotes.  Cada iteracion evalua
        de golpe toda la frontera de pares (particula, nodo) con
        operaciones de numpy: se calcula el vector distancia, se
        deciden los nodos que se abren y se expanden sus ocho hijos.
        El numero de iteraciones es la profundidad del arbol.

        `desplazamiento` es el indice global de la primera particula
        de `posiciones` cuando el arreglo viene recortado (lo usan los
        trabajadores para excluir bien la particula propia).
        """

        posiciones = np.asarray(posiciones, dtype=np.float64)
        salida = np.zeros_like(posiciones)

        for inicio in range(0, len(posiciones), lote):
            fin = min(inicio + lote, len(posiciones))
            self._aceleraciones_lote(
                posiciones[inicio:fin],
                salida[inicio:fin],
                theta, G, epsilon,
                desplazamiento=desplazamiento + inicio,
            )

        return salida

    def _aceleraciones_lote(self, posiciones, salida, theta, G, epsilon,
                            desplazamiento=0):
        masa = self.masa
        com = self.com
        hijos = self.hijos
        centro = self.centro
        lado = self.lado
        part = self.part

        theta_cuadrado = theta * theta
        epsilon_cuadrado = epsilon * epsilon

        n_particulas = len(posiciones)

        # Frontera inicial: todas las particulas en la raiz.
        # Los indices de particula son locales al lote para recorrer
        # posiciones y salida; globales para comparar con part[].
        pares_particula = np.arange(n_particulas, dtype=np.int64)
        pares_globales = np.arange(
            desplazamiento, desplazamiento + n_particulas, dtype=np.int64
        )
        pares_nodo = np.zeros(n_particulas, dtype=np.int64)

        while pares_particula.size:

            distancia = com[pares_nodo] - posiciones[pares_particula]
            distancia_cuadrada = (
                distancia[:, 0] ** 2
                + distancia[:, 1] ** 2
                + distancia[:, 2] ** 2
            )
            masas = masa[pares_nodo]
            es_hoja = hijos[pares_nodo, 0] == -1

            # --------------------------------------------------
            # hojas: aportan su masa salvo que sean la propia
            # particula
            # --------------------------------------------------
            hojas = np.nonzero(es_hoja)[0]

            if hojas.size:
                mismas = (
                    part[pares_nodo[hojas]]
                    == pares_globales[hojas]
                )
                hojas = hojas[~mismas]

            if hojas.size:
                self._aportar(
                    salida,
                    pares_particula[hojas],
                    masas[hojas],
                    distancia[hojas],
                    distancia_cuadrada[hojas],
                    G, epsilon_cuadrado,
                )

            # --------------------------------------------------
            # nodos internos: abrir o aproximar
            # --------------------------------------------------
            internos = np.nonzero(~es_hoja)[0]

            if not internos.size:
                break

            nodos = pares_nodo[internos]
            lados = lado[nodos]

            distancia_interna = distancia_cuadrada[internos]
            abre = distancia_interna == 0.0

            # Si la particula esta dentro de la celda, la distancia al
            # centro de masa es como maximo la diagonal del cubo, y
            # para theta <= 1/sqrt(3) eso ya obliga a abrir el nodo.
            # El test de contencion solo hace falta con theta grandes.
            if theta > 0.5773502691896258:
                centros = centro[nodos]
                contiene = np.all(
                    np.abs(
                        posiciones[pares_particula[internos]] - centros
                    )
                    <= (lados * 0.5)[:, None],
                    axis=1,
                )
                abre = abre | contiene

            abre = (
                abre
                | (lados * lados >= theta_cuadrado * distancia_interna)
            )

            # Los que no se abren se aproximan por su centro de masa.
            aproximados = internos[~abre]

            if aproximados.size:
                self._aportar(
                    salida,
                    pares_particula[aproximados],
                    masas[aproximados],
                    distancia[aproximados],
                    distancia_cuadrada[aproximados],
                    G, epsilon_cuadrado,
                )

            # Los que se abren expanden sus ocho hijos.
            abiertos = internos[abre]

            if not abiertos.size:
                break

            hijos_abiertos = hijos[pares_nodo[abiertos]].reshape(-1)
            nuevas_particulas = np.repeat(pares_particula[abiertos], 8)
            nuevos_globales = np.repeat(pares_globales[abiertos], 8)

            # Solo siguen los hijos con masa.
            vivos = masa[hijos_abiertos] != 0.0

            pares_nodo = hijos_abiertos[vivos]
            pares_particula = nuevas_particulas[vivos]
            pares_globales = nuevos_globales[vivos]

    @staticmethod
    def _aportar(salida, indice, masas, vector, distancia_cuadrada,
                 G, epsilon_cuadrado):
        """Suma las aceleraciones de un conjunto de pares."""

        factor = G * masas / (
            distancia_cuadrada + epsilon_cuadrado
        ) ** 1.5

        for eje in range(3):
            salida[:, eje] += np.bincount(
                indice,
                weights=factor * vector[:, eje],
                minlength=salida.shape[0],
            )

    def aceleraciones_python(self, posiciones, theta=0.3, G=1.0,
                             epsilon=0.01):
        """
        Mismo recorrido escrito como bucle de Python.

        Existe solo para comprobar que la version vectorizada da lo
        mismo; no se usa en la simulacion.  Los arreglos se pasan a
        listas porque el acceso a elementos de numpy es mas lento
        dentro de un bucle.
        """

        posiciones = np.asarray(posiciones, dtype=np.float64)

        n = self.n_nodos
        masa_l = self.masa[:n].tolist()
        com_l = self.com[:n].tolist()
        centro_l = self.centro[:n].tolist()
        lado_l = self.lado[:n].tolist()
        hijos_l = self.hijos[:n].tolist()
        part_l = self.part[:n].tolist()

        theta_cuadrado = theta * theta
        epsilon_cuadrado = epsilon * epsilon

        salida = [[0.0, 0.0, 0.0] for _ in range(len(posiciones))]
        pila = []

        for indice, posicion in enumerate(posiciones):
            px = posicion[0]
            py = posicion[1]
            pz = posicion[2]

            ax = 0.0
            ay = 0.0
            az = 0.0

            pila.append(0)

            while pila:
                nodo = pila.pop()
                masa = masa_l[nodo]

                if masa == 0.0:
                    continue

                centro_masa = com_l[nodo]
                dx = centro_masa[0] - px
                dy = centro_masa[1] - py
                dz = centro_masa[2] - pz
                distancia_cuadrada = dx * dx + dy * dy + dz * dz

                hijos = hijos_l[nodo]

                if hijos[0] == -1:
                    # Hoja: si es la propia particula no aporta nada.
                    if part_l[nodo] == indice:
                        continue

                    factor = G * masa / (
                        distancia_cuadrada + epsilon_cuadrado
                    ) ** 1.5
                    ax += factor * dx
                    ay += factor * dy
                    az += factor * dz
                    continue

                lado = lado_l[nodo]
                centro_geometrico = centro_l[nodo]

                contiene = (
                    abs(px - centro_geometrico[0]) <= lado * 0.5
                    and abs(py - centro_geometrico[1]) <= lado * 0.5
                    and abs(pz - centro_geometrico[2]) <= lado * 0.5
                )

                if (
                    contiene
                    or distancia_cuadrada == 0.0
                    or lado * lado >= theta_cuadrado * distancia_cuadrada
                ):
                    for k in range(8):
                        hijo = hijos[k]
                        if hijo != -1 and masa_l[hijo] != 0.0:
                            pila.append(hijo)
                    continue

                factor = G * masa / (
                    distancia_cuadrada + epsilon_cuadrado
                ) ** 1.5
                ax += factor * dx
                ay += factor * dy
                az += factor * dz

            salida[indice][0] = ax
            salida[indice][1] = ay
            salida[indice][2] = az

        return np.array(salida, dtype=np.float64)
