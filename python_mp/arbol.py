"""
Octree plano (arreglos) para Barnes-Hut.

Sigue docs/contrato_arbol_plano.md, igual que el nucleo C:

    centro[n, 3]   centro geometrico del cubo
    lado[n]        arista del cubo
    masa[n]        masa total de la celda
    com[n, 3]      centro de masa de la celda
    hijos[n, 8]    hijo por octante, -1 si el octante esta vacio
    part[n]        particula de la hoja (>= 0), -1 si el nodo es interno

Invariantes:
  - nodo 0 = raiz;
  - numeracion en PREORDEN: cada nodo va antes que sus descendientes y el
    subarbol de un nodo es un rango contiguo de indices (habilita el LET
    de MPI);
  - solo se crean octantes ocupados (no hay nodos vacios);
  - MAX_PROFUNDIDAD = 20 con bucketing: en el tope una hoja puede guardar
    varias particulas; part[n] es la primera en orden Morton.

Construccion (sin bucles de Python por particula ni por nodo):
  1. codigos Morton de 63 bits (21 por eje, x el mas significativo) y
     ordenamiento estable -> mismo orden que ordenar_por_morton() en C;
  2. nivel por nivel, los nodos de profundidad d son los grupos de
     particulas con el mismo prefijo de 3*d bits; un nodo es hoja si
     tiene una particula o si d == MAX_PROFUNDIDAD;
  3. ordenar los nodos por (codigo de inicio, profundidad) da el preorden;
  4. como cada nodo es un rango contiguo de particulas en orden Morton,
     masa y centro de masa salen de una suma por rangos (reduceat).
"""

import numpy as np

MAX_PROFUNDIDAD = 20      # contrato, seccion 5
MORTON_BITS = 21          # bits por eje: 3 * 21 = 63 bits
MARGEN = 1e-6             # fraccion de la arista agregada al cubo raiz


# ------------------------------------------------------------
# Morton (mismas operaciones que c/common/morton.c)
# ------------------------------------------------------------

def _separar_bits(v):
    """...b2 b1 b0 -> ...b2 0 0 b1 0 0 b0, sobre uint64."""
    x = v.astype(np.uint64) & np.uint64(0x1fffff)
    x = (x | (x << np.uint64(32))) & np.uint64(0x1f00000000ffff)
    x = (x | (x << np.uint64(16))) & np.uint64(0x1f0000ff0000ff)
    x = (x | (x << np.uint64(8))) & np.uint64(0x100f00f00f00f00f)
    x = (x | (x << np.uint64(4))) & np.uint64(0x10c30c30c30c30c3)
    x = (x | (x << np.uint64(2))) & np.uint64(0x1249249249249249)
    return x


def _juntar_bits(x):
    """Inversa de _separar_bits: extrae un eje de un codigo Morton."""
    x = x & np.uint64(0x1249249249249249)
    x = (x | (x >> np.uint64(2))) & np.uint64(0x10c30c30c30c30c3)
    x = (x | (x >> np.uint64(4))) & np.uint64(0x100f00f00f00f00f)
    x = (x | (x >> np.uint64(8))) & np.uint64(0x1f0000ff0000ff)
    x = (x | (x >> np.uint64(16))) & np.uint64(0x1f00000000ffff)
    x = (x | (x >> np.uint64(32))) & np.uint64(0x1fffff)
    return x


def caja_envolvente(posiciones):
    """Cubo raiz: (minimo[3], lado), igual que caja_envolvente() en C."""
    minimo = posiciones.min(axis=0)
    maximo = posiciones.max(axis=0)
    arista = float((maximo - minimo).max())
    if arista == 0.0:
        arista = 1.0          # todas las particulas en el mismo punto
    margen = MARGEN * arista
    return minimo - margen, arista + 2.0 * margen


def codigos_morton(posiciones, minimo, lado):
    escala = float(1 << MORTON_BITS)
    t = (posiciones - minimo) / lado * escala
    t = np.clip(t, 0.0, escala - 1.0).astype(np.uint64)
    return ((_separar_bits(t[:, 0]) << np.uint64(2))
            | (_separar_bits(t[:, 1]) << np.uint64(1))
            | _separar_bits(t[:, 2]))


class ArbolPlano:

    def __init__(self, capacidad=1024, max_profundidad=MAX_PROFUNDIDAD):
        self.max_profundidad = max_profundidad
        self.capacidad = 0
        self.n_nodos = 0
        self.compartido = False
        self.centro = None
        self.lado = None
        self.masa = None
        self.com = None
        self.hijos = None
        self.part = None
        self._reservar(capacidad)

    @classmethod
    def con_arreglos(cls, centro, lado, masa, com, hijos, part,
                     n_nodos, capacidad, compartido=False,
                     construible=False, max_profundidad=MAX_PROFUNDIDAD):
        """
        Arbol que apunta a arreglos ya existentes.

        Es lo que usan los trabajadores para recorrer el arbol que
        construyo el proceso principal en memoria compartida, sin
        copiar nada.  `construible` se conserva por compatibilidad: la
        construccion por Morton no necesita arreglos auxiliares.
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
        self.centro = np.zeros((nueva, 3), dtype=np.float64)
        self.lado = np.zeros(nueva, dtype=np.float64)
        self.masa = np.zeros(nueva, dtype=np.float64)
        self.com = np.zeros((nueva, 3), dtype=np.float64)
        self.hijos = np.full((nueva, 8), -1, dtype=np.int32)
        self.part = np.full(nueva, -1, dtype=np.int32)
        self.capacidad = nueva

    # ------------------------------------------------------------
    # construccion
    # ------------------------------------------------------------

    def construir(self, posiciones, masas):
        """Construye el arbol completo a partir de las posiciones."""

        posiciones = np.asarray(posiciones, dtype=np.float64)
        masas = np.asarray(masas, dtype=np.float64)
        n = len(posiciones)
        prof_max = min(self.max_profundidad, MORTON_BITS - 1)

        # 1. Morton y orden (estable: empates por indice, como en C).
        minimo, lado_raiz = caja_envolvente(posiciones)
        codigos = codigos_morton(posiciones, minimo, lado_raiz)
        orden = np.argsort(codigos, kind="stable")
        codigos_ord = codigos[orden]

        # 2. Nodos nivel por nivel.  `activas` son posiciones en el
        #    arreglo ordenado de las particulas que aun no estan en una
        #    hoja.  Por nivel se guarda: inicio, cantidad y prefijo.
        niveles = []
        activas = np.arange(n, dtype=np.int64)
        for d in range(prof_max + 1):
            desplaz = np.uint64(3 * (MORTON_BITS - d))
            prefijo = codigos_ord[activas] >> desplaz
            nuevo = np.empty(len(activas), dtype=bool)
            nuevo[0] = True
            nuevo[1:] = prefijo[1:] != prefijo[:-1]
            arranques = np.flatnonzero(nuevo)
            cantidades = np.diff(np.append(arranques, len(activas)))

            es_hoja = cantidades == 1
            if d == prof_max:
                es_hoja[:] = True            # bucketing en el tope
            niveles.append((activas[arranques], cantidades,
                            prefijo[arranques], es_hoja))

            activas = activas[~np.repeat(es_hoja, cantidades)]
            if activas.size == 0:
                break

        inicio = np.concatenate([nv[0] for nv in niveles])
        cantidad = np.concatenate([nv[1] for nv in niveles])
        es_hoja = np.concatenate([nv[3] for nv in niveles])
        prof = np.concatenate([np.full(len(nv[0]), d, dtype=np.int64)
                               for d, nv in enumerate(niveles)])
        n_nodos = len(inicio)

        # 3. Preorden = orden por (codigo de inicio del cubo, profundidad):
        #    un padre comparte codigo de inicio con su primer hijo pero
        #    tiene menor profundidad; los descendientes de un nodo caen
        #    antes que su siguiente hermano.
        cod_inicio = np.concatenate([
            nv[2] << np.uint64(3 * (MORTON_BITS - d))
            for d, nv in enumerate(niveles)
        ])
        preorden = np.lexsort((prof, cod_inicio))   # viejo -> posicion
        nuevo_id = np.empty(n_nodos, dtype=np.int64)
        nuevo_id[preorden] = np.arange(n_nodos)

        self._reservar(n_nodos)
        self.n_nodos = n_nodos

        # Geometria: el cubo de profundidad d con codigo de inicio c.
        lado = lado_raiz / np.exp2(prof)
        esquina = np.stack([_juntar_bits(cod_inicio >> np.uint64(s))
                            for s in (2, 1, 0)], axis=1)
        celda = lado_raiz / float(1 << MORTON_BITS)
        centro = minimo + esquina.astype(np.float64) * celda \
            + 0.5 * lado[:, None]

        # 4. Masa y centro de masa por rangos de particulas.
        masas_ord = masas[orden]
        tabla = np.zeros((n + 1, 4), dtype=np.float64)
        tabla[:n, 0] = masas_ord
        tabla[:n, 1:] = masas_ord[:, None] * posiciones[orden]
        cortes = np.empty(2 * n_nodos, dtype=np.int64)
        cortes[0::2] = inicio
        cortes[1::2] = inicio + cantidad
        sumas = np.add.reduceat(tabla, cortes, axis=0)[0::2]
        masa = sumas[:, 0]
        con_masa = masa > 0.0
        com = centro.copy()          # celda sin masa: com = centro
        com[con_masa] = sumas[con_masa, 1:] / masa[con_masa, None]

        # Escritura en los arreglos del contrato (posiblemente
        # compartidos), ya en preorden.
        k = nuevo_id
        self.centro[k] = centro
        self.lado[k] = lado
        self.masa[k] = masa
        self.com[k] = com
        self.hijos[:n_nodos] = -1
        self.part[:n_nodos] = -1
        self.part[k[es_hoja]] = orden[inicio[es_hoja]]

        # hijos[padre, octante]: el padre de un nodo de profundidad d es
        # el nodo de profundidad d-1 cuyo prefijo es el suyo >> 3.
        base = 0
        for d in range(1, len(niveles)):
            prefijos_padre = niveles[d - 1][2]
            base_padre = base
            base += len(prefijos_padre)
            prefijo = niveles[d][2]
            padre = base_padre + np.searchsorted(
                prefijos_padre, prefijo >> np.uint64(3))
            octante = (prefijo & np.uint64(7)).astype(np.int64)
            hijo = base + np.arange(len(prefijo))
            self.hijos[k[padre], octante] = k[hijo]

        return self

    def profundidad(self):
        """Profundidad maxima alcanzada (para el informe y las pruebas)."""
        lados = self.lado[:self.n_nodos]
        return int(np.round(np.log2(lados[0] / lados.min())))

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
        deciden los nodos que se abren y se expanden sus hijos.
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
            hoja_part = part[pares_nodo]
            es_hoja = hoja_part >= 0

            # --------------------------------------------------
            # hojas: aportan su masa salvo que sean la propia
            # particula
            # --------------------------------------------------
            hojas = np.nonzero(es_hoja & (hoja_part != pares_globales))[0]

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

            # Los que se abren expanden sus hijos existentes (los
            # octantes vacios tienen -1 y se descartan).
            abiertos = internos[abre]

            if not abiertos.size:
                break

            hijos_abiertos = hijos[pares_nodo[abiertos]].reshape(-1)
            nuevas_particulas = np.repeat(pares_particula[abiertos], 8)
            nuevos_globales = np.repeat(pares_globales[abiertos], 8)

            vivos = hijos_abiertos >= 0

            pares_nodo = hijos_abiertos[vivos].astype(np.int64)
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
        Mismo recorrido escrito como bucle de Python con pila explicita.

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

                if part_l[nodo] >= 0:
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
                    for hijo in hijos_l[nodo]:
                        if hijo >= 0:
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
