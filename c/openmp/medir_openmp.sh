#!/usr/bin/env bash
# Mide la version OpenMP con N = 200 000 y guarda un CSV.
#
# Uso (desde la raiz del repositorio):
#   bash c/openmp/medir_openmp.sh [N] [repeticiones]
#
# Mide: 1, 2, 4 y 8 hilos con las distribuciones uniforme y concentradas
# (reparto dynamic,64), y ademas static vs dynamic vs guided con el
# maximo de hilos. Cada corrida hace 1 paso (= 2 calculos de fuerzas).
set -e
N=${1:-200000}
REP=${2:-3}
PASOS=1; DT=0.001; THETA=0.5; EPS=0.01
SIM=./c/openmp/sim_openmp
SALIDA=resultados/openmp_N${N}.csv
mkdir -p resultados

for d in uniforme concentradas; do
    [ -f /tmp/${d}_${N}.txt ] || \
        python3 data/generate_particles.py --distribucion $d --N $N \
            --out /tmp/${d}_${N}.txt > /dev/null
done

# Extrae "clave=valor" de la salida del programa.
valor() { grep "^$1=" <<< "$2" | cut -d= -f2; }

correr() {   # correr <distribucion> <hilos> <reparto>
    for r in $(seq 1 $REP); do
        out=$(OMP_SCHEDULE="$3" $SIM /tmp/$1_${N}.txt $N $PASOS $DT $THETA $EPS $2 /tmp/acc_omp.txt)
        linea="$1,$2,\"$3\",$r,$(valor tiempo_build "$out"),$(valor tiempo_sort "$out"),$(valor tiempo_walk "$out"),$(valor tiempo_integrate "$out"),$(valor tiempo_total "$out")"
        echo "$linea" | tee -a $SALIDA
    done
}

echo "distribucion,hilos,reparto,rep,build,sort,walk,integrate,total" | tee $SALIDA
for d in uniforme concentradas; do
    for h in 1 2 4 8 16; do
        correr $d $h "dynamic,64"
    done
done
for rep in static guided; do
    correr concentradas 16 "$rep"
done
echo "Guardado en $SALIDA"