"""
Paralelismo con multiprocessing y memoria compartida.

Dos estrategias, medibles entre si:

  compartido  el proceso principal construye el arbol una sola vez
              dentro de bloques shared_memory y todos los trabajadores
              solo lo recorren;
  local       cada trabajador reconstruye su propio arbol a partir de
              las posiciones, que tambien estan en memoria compartida.

La estrategia se elige con la variable de entorno BH_ESTRATEGIA
(compartido por defecto) y el numero de trabajadores con el argumento
<hilos> del contrato comun.
"""

import os
import sys
import time

import numpy as np
from multiprocessing import Pool, shared_memory

from arbol import ArbolPlano

ESTRATEGIA_POR_DEFECTO = "compartido"

# Estado de cada trabajador, llenado por el inicializador.
_ESTADO = {}


def _inicializar(estrategia, bloques, theta, G, epsilon, lote, capacidad):
    """
    Adjunta los bloques de memoria compartida en el trabajador.

    bloques: {clave: (nombre, forma, tipo)} con la informacion
    necesaria para reconstruir cada arreglo numpy.
    """

    global _ESTADO

    arreglos = {}
    memorias = []

    for clave, (nombre, forma, tipo) in bloques.items():
        memoria = shared_memory.SharedMemory(name=nombre, create=False)
        arreglo = np.ndarray(
            tuple(forma), dtype=np.dtype(tipo), buffer=memoria.buf
        )
        memorias.append(memoria)
        arreglos[clave] = arreglo

    _ESTADO = {
        "estrategia": estrategia,
        "theta": theta,
        "G": G,
        "epsilon": epsilon,
        "lote": lote,
        "capacidad": capacidad,
        "arreglos": arreglos,
        "memorias": memorias,
    }


def _arbol_compartido():
    """Arbol apuntando directamente a los bloques compartidos."""

    arreglos = _ESTADO["arreglos"]

    return ArbolPlano.con_arreglos(
        centro=arreglos["centro"],
        lado=arreglos["lado"],
        masa=arreglos["masa"],
        com=arreglos["com"],
        hijos=arreglos["hijos"],
        part=arreglos["part"],
        n_nodos=int(arreglos["meta"][0]),
        capacidad=int(arreglos["meta"][1]),
    )


def _trabajador(tarea):
    """
    Calcula las aceleraciones de un tramo de particulas.

    Devuelve los segundos que tardo en construir y en recorrer, para
    que el proceso principal pueda desglosar las fases.
    """

    inicio, fin = tarea
    arreglos = _ESTADO["arreglos"]

    posiciones = arreglos["posiciones"]
    salida = arreglos["salida"]

    momento = time.perf_counter()

    if _ESTADO["estrategia"] == "compartido":
        arbol = _arbol_compartido()
        tiempo_construir = 0.0
    else:
        arbol = ArbolPlano(capacidad=_ESTADO["capacidad"])
        arbol.construir(posiciones, arreglos["masas"])
        tiempo_construir = time.perf_counter() - momento

    momento = time.perf_counter()
    cpu = time.process_time()
    salida[inicio:fin] = arbol.aceleraciones(
        posiciones[inicio:fin],
        theta=_ESTADO["theta"],
        G=_ESTADO["G"],
        epsilon=_ESTADO["epsilon"],
        lote=_ESTADO["lote"],
        desplazamiento=inicio,
    )
    tiempo_recorrer = time.perf_counter() - momento
    tiempo_cpu = time.process_time() - cpu

    return tiempo_construir, tiempo_recorrer, tiempo_cpu


def _dividir(n, partes):
    """Reparte [0, n) en `partes` tramos lo mas equilibrados posible."""

    base = n // partes
    resto = n % partes
    tareas = []
    inicio = 0

    for indice in range(partes):
        tamano = base + (1 if indice < resto else 0)
        if tamano:
            tareas.append((inicio, inicio + tamano))
        inicio += tamano

    return tareas


class Paralelo:
    """Ejecuta los calculos de fuerzas en varios procesos."""

    def __init__(self, posiciones, masas, hilos, theta=0.3, G=1.0,
                 epsilon=0.01, lote=2048, estrategia=None):

        if estrategia is None:
            estrategia = os.environ.get(
                "BH_ESTRATEGIA", ESTRATEGIA_POR_DEFECTO
            )

        if estrategia not in ("compartido", "local"):
            raise SystemExit(
                "BH_ESTRATEGIA debe ser 'compartido' o 'local', no %r"
                % estrategia
            )

        self.hilos = max(1, int(hilos))
        self.estrategia = estrategia
        self.theta = theta
        self.G = G
        self.epsilon = epsilon
        self.lote = lote
        self.n = len(posiciones)

        # Capacidad del arbol: se mide con una construccion previa y
        # se deja holgura para los siguientes pasos.
        prueba = ArbolPlano()
        prueba.construir(posiciones, masas)
        self.capacidad = int(prueba.n_nodos * 1.5) + 64

        formas = {
            "posiciones": (self.n, 3),
            "masas": (self.n,),
            "salida": (self.n, 3),
            "meta": (2,),
        }

        if estrategia == "compartido":
            formas.update({
                "centro": (self.capacidad, 3),
                "lado": (self.capacidad,),
                "masa": (self.capacidad,),
                "com": (self.capacidad, 3),
                "hijos": (self.capacidad, 8),
                "part": (self.capacidad,),
            })

        self._memorias = []
        self._arreglos = {}

        for clave, forma in formas.items():
            tipo = np.int64 if clave == "meta" else np.float64
            if clave in ("hijos", "part"):
                tipo = np.int32
            arreglo, memoria = self._crear_bloque(clave, forma, tipo)
            self._arreglos[clave] = arreglo
            self._memorias.append(memoria)

        self._arreglos["posiciones"][:] = posiciones
        self._arreglos["masas"][:] = masas
        self._arreglos["meta"][0] = 0
        self._arreglos["meta"][1] = self.capacidad

        self.arbol = None
        if estrategia == "compartido":
            self.arbol = ArbolPlano.con_arreglos(
                centro=self._arreglos["centro"],
                lado=self._arreglos["lado"],
                masa=self._arreglos["masa"],
                com=self._arreglos["com"],
                hijos=self._arreglos["hijos"],
                part=self._arreglos["part"],
                n_nodos=0,
                capacidad=self.capacidad,
                compartido=True,
                construible=True,
            )

        bloques = {
            clave: (
                self._nombres[clave],
                forma,
                str(self._arreglos[clave].dtype),
            )
            for clave, forma in formas.items()
        }

        self.tareas = _dividir(self.n, self.hilos)

        self.pool = Pool(
            processes=self.hilos,
            initializer=_inicializar,
            initargs=(
                estrategia,
                bloques,
                theta,
                G,
                epsilon,
                lote,
                self.capacidad,
            ),
        )

    # ------------------------------------------------------------

    def _crear_bloque(self, clave, forma, tipo):
        tamano = int(np.prod(forma)) * np.dtype(tipo).itemsize
        nombre = "bh_%s_%d_%d" % (clave, os.getpid(), id(self))
        memoria = shared_memory.SharedMemory(
            name=nombre, create=True, size=tamano
        )
        arreglo = np.ndarray(forma, dtype=tipo, buffer=memoria.buf)

        if not hasattr(self, "_nombres"):
            self._nombres = {}
        self._nombres[clave] = memoria.name

        return arreglo, memoria

    # ------------------------------------------------------------

    def construir(self, posiciones):
        """
        Construye el arbol en memoria compartida.

        Con la estrategia local no hace nada: cada trabajador
        construye el suyo.  Devuelve los segundos del proceso
        principal.
        """

        if self.estrategia != "compartido":
            return 0.0

        self._arreglos["posiciones"][:] = posiciones

        momento = time.perf_counter()
        self.arbol.construir(posiciones, self._arreglos["masas"])
        tiempo = time.perf_counter() - momento

        if self.arbol.n_nodos > self.capacidad:
            raise SystemExit(
                "el arbol necesito %d nodos y la capacidad era %d"
                % (self.arbol.n_nodos, self.capacidad)
            )

        self._arreglos["meta"][0] = self.arbol.n_nodos

        return tiempo

    def recorrer(self, posiciones=None):
        """
        Reparte el recorrido entre los trabajadores.

        Devuelve (aceleraciones, tiempos) donde tiempos es la lista de
        (construir, recorrer) de cada tarea.
        """

        if posiciones is not None:
            self._arreglos["posiciones"][:] = posiciones

        momento = time.perf_counter()
        tiempos = self.pool.map(_trabajador, self.tareas)
        tiempo_pared = time.perf_counter() - momento

        if os.environ.get("BH_DEBUG") == "1":
            detalle = " ".join(
                "%.2f/%.2f/cpu%.2f" % (construir, recorrer, cpu)
                for construir, recorrer, cpu in tiempos
            )
            print(
                "# trabajadores (construir/recorrer/cpu): %s  pared=%.3f"
                % (detalle, tiempo_pared),
                file=sys.stderr,
            )

        aceleraciones = self._arreglos["salida"].copy()

        return aceleraciones, tiempos, tiempo_pared

    def cerrar(self):
        self.pool.close()
        self.pool.join()

        for memoria in self._memorias:
            try:
                memoria.close()
                memoria.unlink()
            except FileNotFoundError:
                pass
