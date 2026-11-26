#!/usr/bin/env bash
# =====================================================================
# demo.sh — todo el proyecto con un solo comando (para la demo en vivo)
#
#   bash demo.sh <archivo> [opciones]
#
# Hace, en orden:
#   1. compila C, OpenMP y MPI (avisa si algo no compila);
#   2. diagnostica la distribucion (concentracion, choques, riesgo MPI);
#   3. verifica cada version contra la suma directa con una muestra de
#      2000 particulas (la suma directa con N grande no es viable);
#   4. mide cada version sobre el archivo completo y muestra una tabla.
# Todo queda tambien en resultados/demo_<fecha>.log.
#
# Opciones:
#   -t THETA        criterio de apertura           (0.5)
#   -d DT           paso de tiempo                 (0.001)
#   -s PASOS        pasos de simulacion            (1)
#   -P "1 4 8"      procesos / hilos a medir       ("1 4 8")
#   -r REPS         repetir cada medicion y tomar la mejor (1)
#   -b              MPI con biseccion en x en vez de Morton
#   -c              MPI con balance por particulas en vez de por costo
#   -y              incluir Python-mp (lento: ~25 s por evaluacion con 200k)
#   -v              solo diagnostico y verificacion (sin medir)
#   -m              solo medir (sin verificar)
#
# Ejemplos:
#   bash demo.sh datos.txt
#   bash demo.sh datos.txt -t 0.7 -P "1 2 4 8" -b
#
# CUDA se corre en Colab: cuda/barnes_hut_cuda.ipynb (misma interfaz:
#   ./cuda/sim_cuda <archivo> <N> <pasos> <dt> <theta> <eps> <hilos/bloque>).
# =====================================================================
set -u
cd "$(dirname "$0")"

THETA=0.5; DT=0.001; PASOS=1; EPS=0.01; PROCS="1 4 8"
PARTICION=morton; BALANCE=costo; PYTHON=0; VERIFICAR=1; MEDIR=1; REPS=1

ARCHIVO=${1:-}
if [ -z "$ARCHIVO" ] || [ ! -f "$ARCHIVO" ]; then
    sed -n '3,30p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
fi
shift
while getopts "t:d:s:P:r:bcyvm" op; do
    case $op in
        t) THETA=$OPTARG ;; d) DT=$OPTARG ;; s) PASOS=$OPTARG ;; r) REPS=$OPTARG ;;
        P) PROCS=$OPTARG ;; b) PARTICION=biseccion ;; c) BALANCE=particulas ;;
        y) PYTHON=1 ;; v) MEDIR=0 ;; m) VERIFICAR=0 ;;
        *) exit 1 ;;
    esac
done
N=$(grep -c '[0-9]' "$ARCHIVO")
mkdir -p resultados
LOG=resultados/demo_$(date +%Y%m%d_%H%M%S).log
exec > >(tee "$LOG") 2>&1

titulo() { printf '\n\033[1m== %s\033[0m\n' "$*"; }

# Colocacion de procesos MPI (4 nucleos fisicos, 8 hilos logicos):
mpi_opts() {
    if [ "$1" -le 4 ]; then echo ""
    elif [ "$1" -le 8 ]; then echo "--use-hwthread-cpus"
    else echo "--oversubscribe --bind-to none"; fi
}
valor() { grep "^$1=" <<< "$2" | head -1 | cut -d= -f2; }

echo "archivo=$ARCHIVO  N=$N  theta=$THETA  dt=$DT  pasos=$PASOS  eps=$EPS"
echo "procesos=[$PROCS]  MPI: particion=$PARTICION balance=$BALANCE"

# ---------------------------------------------------------------- 1
titulo "1. Compilacion"
ok=1
for d in c/common c/openmp c/mpi; do
    if out=$(make -s -C "$d" 2>&1); then
        w=$(grep -ci warning <<< "$out" || true)
        echo "  $d: OK${w:+ ($w warnings)}"
    else
        echo "  $d: ERROR"; echo "$out" | tail -15; ok=0
    fi
done
[ $ok = 1 ] || { echo "Arreglar la compilacion antes de seguir."; exit 1; }

# ---------------------------------------------------------------- 2
titulo "2. Diagnostico de la distribucion"
python3 tools/diagnostico.py "$ARCHIVO" --procesos $(echo $PROCS | tr ' ' '\n' | awk '$1>1') \
    --muestra 2000 --salida /tmp/demo_muestra.txt

# ---------------------------------------------------------------- 3
if [ $VERIFICAR = 1 ]; then
    titulo "3. Verificacion (muestra de 2000 particulas, theta=$THETA, dt=$DT, 5 pasos)"
    P4=$(echo $PROCS | tr ' ' '\n' | awk '$1>1' | head -1); P4=${P4:-4}
    verif() {   # verif <nombre> <programa> <hilos>
        r=$(python3 tools/verificar.py --particulas /tmp/demo_muestra.txt \
              --programa "$2" --hilos "$3" --theta "$THETA" --dt "$DT" 2>&1)
        e=$(grep "^  error promedio" <<< "$r" | awk '{print $3}')
        m=$(grep "^  error maximo" <<< "$r" | awk '{print $3}')
        g=$(grep "^  derivacion" <<< "$r" | awk '{print $4}')
        res=$(grep -E "^(OK|FALLO)" <<< "$r" | sed 's/:.*//')
        printf "  %-16s error prom %8s %%  max %8s %%  energia %10s %%  %s\n" \
               "$1" "$e" "$m" "$g" "$res"
        grep -E "^  - " <<< "$r" | sed 's/^/      /'
    }
    verif "C secuencial"     "./c/common/sim_barnes_hut" 1
    verif "OpenMP ($P4)"     "./c/openmp/sim_openmp" "$P4"
    verif "MPI ($P4)"        "env BH_PARTICION=$PARTICION BH_BALANCE=$BALANCE mpirun -np $P4 $(mpi_opts $P4) ./c/mpi/sim_mpi" 1
    [ $PYTHON = 1 ] && verif "Python-mp ($P4)" "python3 python_mp/sim.py" "$P4"
    echo "  Si falla solo la ENERGIA con particulas muy juntas: bajar dt (p. ej. -d 0.00025);"
    echo "  la suma directa deriva igual, es el integrador y no Barnes-Hut."
    if awk -v t="$THETA" 'BEGIN{exit !(t > 0.5)}'; then
        echo "  Con theta=$THETA > 0.5 un error promedio algo mayor que 1 % es esperado:"
        echo "  theta mas grande abre menos celdas (mas rapido, menos preciso). Lo que"
        echo "  importa es que todas las versiones den el MISMO error que el C secuencial."
    fi
fi

# ---------------------------------------------------------------- 4
if [ $MEDIR = 1 ]; then
    titulo "4. Rendimiento sobre el archivo completo (N=$N)"
    printf "  %-13s %3s %9s %8s %8s %8s %8s  %s\n" \
           version P total_s build_s walk_s speedup com_% "extra"
    fila() {   # fila <version> <P> <salida> <t_ref>
        local t b w c ex
        t=$(valor tiempo_total "$3"); b=$(valor tiempo_build "$3"); w=$(valor tiempo_walk "$3")
        c=$(awk -v r="$(valor tiempo_reparto "$3")" -v l="$(valor tiempo_let "$3")" -v t="$t" \
              'BEGIN{ if (r=="") print "-"; else printf "%.1f", 100*(r+l)/t }')
        ex=""
        if [ -n "$(valor desbalance_walk "$3")" ]; then
            ex="desbalance $(valor desbalance_walk "$3"), LET $(valor nodos_let_recibidos_max "$3") nodos"
        fi
        printf "  %-13s %3s %9.3f %8.3f %8.3f %7.2fx %8s  %s\n" "$1" "$2" "$t" "$b" "$w" \
               "$(awk -v a="$4" -v b="$t" 'BEGIN{print a/b}')" "$c" "$ex"
    }
    # corre "<comando>" <hilos>: contrato comun del README; repite REPS veces
    # y devuelve la salida con el MENOR tiempo_total (filtra el ruido termico).
    corre() {
        local mejor="" tmejor=1e99 o t i
        for i in $(seq "$REPS"); do
            o=$($1 "$ARCHIVO" "$N" "$PASOS" "$DT" "$THETA" "$EPS" "$2" /tmp/demo_acc.txt 2>/dev/null)
            t=$(valor tiempo_total "$o")
            if awk -v a="$t" -v b="$tmejor" 'BEGIN{exit !(a < b)}'; then mejor=$o; tmejor=$t; fi
        done
        echo "$mejor"
    }

    out=$(corre ./c/common/sim_barnes_hut 1)
    fila "C secuencial" 1 "$out" "$(valor tiempo_total "$out")"

    ref=""
    for p in $PROCS; do
        out=$(corre "env OMP_SCHEDULE=dynamic,64 ./c/openmp/sim_openmp" "$p")
        ref=${ref:-$(valor tiempo_total "$out")}; fila OpenMP "$p" "$out" "$ref"
    done
    ref=""
    for p in $PROCS; do
        out=$(corre "env BH_PARTICION=$PARTICION BH_BALANCE=$BALANCE mpirun -np $p $(mpi_opts "$p") ./c/mpi/sim_mpi" "$p")
        ref=${ref:-$(valor tiempo_total "$out")}; fila "MPI" "$p" "$out" "$ref"
    done
    if [ $PYTHON = 1 ]; then
        ref=""
        for p in $PROCS; do
            out=$(corre "python3 python_mp/sim.py" "$p")
            ref=${ref:-$(valor tiempo_total "$out")}; fila "Python-mp" "$p" "$out" "$ref"
        done
    fi
    echo "  tiempos: el mejor de $REPS corrida(s) (-r 3 si hay tiempo: el portatil varia ~25 %)."
    echo "  speedup: contra la misma version con el primer P de la lista."
    echo "  com_%: reparto + LET de MPI.  desbalance: max/promedio del walk."
fi
echo; echo "Registro: $LOG"
