#!/usr/bin/env bash
# Exporta informe/informe.md a DOCX, HTML y PDF y cuenta palabras.
#
#   bash tools/exportar_informe.sh            # docx + html + pdf
#   bash tools/exportar_informe.sh --contar   # solo palabras por seccion
#
# El PDF se imprime con Microsoft Edge de Windows en modo headless (no hace
# falta LaTeX). Salida en informe/salida/ (ignorada por git).
set -e
cd "$(dirname "$0")/.."
MD=${INFORME:-informe/informe.md}   # INFORME=otro.md para probar variantes
OUT=informe/salida
mkdir -p "$OUT"

contar() {
    # Palabras del cuerpo (sin comentarios HTML, encabezado YAML ni referencias).
    python3 - "$MD" <<'EOF'
import re, sys
t = open(sys.argv[1], encoding="utf-8").read()
t = re.sub(r"<!--.*?-->", "", t, flags=re.S)
t = re.sub(r"^---.*?---", "", t, count=1, flags=re.S)
t = t.split("## Referencias")[0]
total = 0
for sec in re.split(r"^## ", t, flags=re.M)[1:]:
    titulo, _, cuerpo = sec.partition("\n")
    cuerpo = re.sub(r"!\[.*?\]\(.*?\)(\{.*?\})?", "", cuerpo)
    n = len(re.findall(r"\w+", cuerpo))
    total += n
    print(f"  {n:5d}  {titulo.strip()}")
print(f"  {total:5d}  TOTAL (cabe hasta ~2000 con 3 figuras y 2 tablas)")
EOF
}

if [ "$1" = "--contar" ]; then contar; exit 0; fi

pandoc "$MD" -o "$OUT/informe.docx" --resource-path=informe
pandoc "$MD" -s -o "$OUT/informe.html" --resource-path=informe \
    --css=informe/apa.css --embed-resources

# Edge no lee bien archivos dentro de \\wsl.localhost: se imprime desde %TEMP%.
EDGE="/mnt/c/Program Files (x86)/Microsoft/Edge/Application/msedge.exe"
rm -f "$OUT/informe.pdf"
if [ -x "$EDGE" ]; then
    WTEMP=$(cmd.exe /c 'echo %TEMP%' 2>/dev/null | tr -d '\r')
    WDIR="$(wslpath "$WTEMP")/bh_informe"
    mkdir -p "$WDIR"
    cp "$OUT/informe.html" "$WDIR/informe.html"
    rm -f "$WDIR/informe.pdf"
    "$EDGE" --headless=new --disable-gpu --no-pdf-header-footer \
        --user-data-dir="$WTEMP\\bh_informe\\perfil" \
        --print-to-pdf="$WTEMP\\bh_informe\\informe.pdf" \
        "file:///$(echo "$WTEMP" | tr '\\' '/')/bh_informe/informe.html" >/dev/null 2>&1 || true
    for _ in 1 2 3 4 5; do [ -f "$WDIR/informe.pdf" ] && break; sleep 1; done
    [ -f "$WDIR/informe.pdf" ] && cp "$WDIR/informe.pdf" "$OUT/informe.pdf"
fi

contar
if [ -f "$OUT/informe.pdf" ]; then
    paginas=$(grep -a -c "/Type /Page$\|/Type/Page[^s]" "$OUT/informe.pdf" || true)
    echo "  PDF: $OUT/informe.pdf  (paginas: $paginas, limite 6)"
fi
echo "  DOCX: $OUT/informe.docx"
