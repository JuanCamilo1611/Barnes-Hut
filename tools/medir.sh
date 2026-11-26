#!/bin/bash
#
# Mide una version con varias cantidades de hilos o procesos y deja el
# resultado en un CSV dentro de resultados/.
#
# Uso:
#   tools/medir.sh <version> <archivo> <N> <pasos> <dt> <theta> <eps> <hilos...>
#
# Ejemplo:
#   tools/medir.sh python_mp data/p20k.txt 20000 3 0.001 0.3 0.01 1 2 4 8
#
# La variable BH_ESTRATEGIA controla la estrategia del paralelismo
# (compartido o local); si se quiere medir las dos, se llama dos veces.

set -euo pipefail

if [ $# -lt 8 ]; then
    echo "uso: tools/medir.sh <version> <archivo> <N> <pasos> <dt> <theta> <eps> <hilos...>"
    echo "versiones: python_mp"
    exit 1
fi

version=$1
archivo=$2
n=$3
pasos=$4
dt=$5
theta=$6
eps=$7
shift 7

case $version in
    python_mp)
        programa=(python3 python_mp/sim.py)
        ;;
    *)
        echo "version desconocida: $version"
        echo "versiones: python_mp"
        exit 1
        ;;
esac

raiz=$(cd "$(dirname "$0")/.." && pwd)
cd "$raiz"

estrategia=${BH_ESTRATEGIA:-secuencial}
marca=$(date +%Y%m%d_%H%M%S)
salida="resultados/medicion_${version}_${estrategia}_${marca}.csv"

mkdir -p resultados
echo "hilos,estrategia,tiempo_build,tiempo_sort,tiempo_walk,tiempo_integrate,tiempo_total" > "$salida"

for hilos in "$@"; do
    echo "==> $version / $estrategia / $hilos proceso(s)"

    salida_acc=$(mktemp)
    registro=$("${programa[@]}" "$archivo" "$n" "$pasos" "$dt" "$theta" "$eps" "$hilos" "$salida_acc")
    rm -f "$salida_acc"

    valor() {
        echo "$registro" | sed -n "s/^$1=//p"
    }

    echo "$hilos,$estrategia,$(valor tiempo_build),$(valor tiempo_sort),$(valor tiempo_walk),$(valor tiempo_integrate),$(valor tiempo_total)" >> "$salida"
done

echo
echo "Resultados en $salida"
column -s, -t "$salida" || cat "$salida"
