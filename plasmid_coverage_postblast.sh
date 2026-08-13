#!/usr/bin/env bash
set -euo pipefail

# Uso:
# ./plasmid_coverage_postblast.sh \
#   -b resultados_uymma_blast.tsv \
#   -i isla_uymma02A.tsv \
#   -o resultados_mosaico_uymma

# --- helper: convertir a ruta absoluta ---
abs_path() {
  case "${1:-}" in
    /*) echo "$1" ;;
    *)  echo "$(pwd)/$1" ;;
  esac
}

BLAST_TSV=""
ISLAND_TSV=""
OUTDIR=""

while getopts ":b:i:o:" opt; do
  case $opt in
    b) BLAST_TSV=$(abs_path "$OPTARG") ;;
    i) ISLAND_TSV=$(abs_path "$OPTARG") ;;
    o) OUTDIR="$OPTARG" ;;
    \?) echo "Opción inválida: -$OPTARG" >&2; exit 1 ;;
    :)  echo "La opción -$OPTARG requiere un argumento." >&2; exit 1 ;;
  esac
done

if [[ -z "${BLAST_TSV:-}" || -z "${ISLAND_TSV:-}" || -z "${OUTDIR:-}" ]]; then
  echo "Faltan argumentos. Uso: $0 -b <blast.tsv> -i <isla.tsv> -o <outdir>" >&2
  exit 1
fi

[[ -f "$BLAST_TSV" ]] || { echo "No existe -b $BLAST_TSV" >&2; exit 1; }
[[ -f "$ISLAND_TSV" ]] || { echo "No existe -i $ISLAND_TSV" >&2; exit 1; }

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# 1. Leer QID (qseqid) y PLEN (qlen) del archivo BLAST
# En tu outfmt: col 1 es qseqid, col 11 es qlen
QID=$(head -n1 "$BLAST_TSV" | cut -f1)
PLEN=$(head -n1 "$BLAST_TSV" | cut -f11)
echo "[info] qseqid=$QID  qlen=$PLEN"

echo "[1/6] HSPs → BED (0-based) + Agrupación por Genoma (WGS)…"
# Mapeo de tus columnas: q=$1, sid=$2, pid=$3, qs=$5, qe=$6
awk 'BEGIN{OFS="\t"}{
  q=$1; sid=$2; pid=$3; qs=$5; qe=$6;
  if(qs<=qe){s=qs-1; e=qe} else {s=qe-1; e=qs}
  
  clean=sid; 
  # Quitar pipes si vienen del NCBI
  if(match(sid, /\|([^|]+)\|/, m)) clean=m[1];
  
  # AGRUPACIÓN WGS: Extraer prefijo de proyecto WGS (ej. JAUYZV02)
  if (match(clean, /^([A-Za-z]{2}_)?[A-Za-z]{4,6}[0-9]{2}/)) {
      clean = substr(clean, RSTART, RLENGTH);
  }
  
  print q, s, e, clean, pid, "."
}' "$BLAST_TSV" | bedtools sort -i - > hsp_query.clean.bed

echo "[2/6] Isla: TSV 1-based → BED (0-based) + merge…"
awk 'BEGIN{OFS="\t"}{print $1,$2-1,$3,$4,$5,$6}' "$ISLAND_TSV" \
| bedtools sort -i - > island_genes.bed

bedtools merge -i island_genes.bed > island_merged.bed
bedtools sort -i island_merged.bed > island_merged.sorted.bed

ISLEN=$(awk '{s+=($3-$2)} END{print s+0}' island_merged.sorted.bed)
echo "[info] island length (nt) = $ISLEN"

ISID=$(cut -f1 island_merged.sorted.bed | head -n1)
if [[ "$ISID" != "$QID" ]]; then
  echo "[warn] island seqid=$ISID difiere de qseqid=$QID → rescribiendo a $QID"
  awk -v Q="$QID" 'BEGIN{OFS="\t"}{$1=Q; print}' island_merged.sorted.bed > tmp && mv tmp island_merged.sorted.bed
fi

echo "[3/6] Cobertura en la isla (única) por sujeto…"
: > coverage_in_island.tsv
cut -f4 hsp_query.clean.bed | sort -u | while read -r sid; do
  awk -v S="$sid" '$4==S{print $1"\t"$2"\t"$3}' hsp_query.clean.bed \
  | bedtools sort -i - \
  | bedtools intersect -sorted -a - -b island_merged.sorted.bed \
  | bedtools sort -i - \
  | bedtools merge -i - \
  | awk -v S="$sid" -v L="$ISLEN" 'BEGIN{len=0}{len+=($3-$2)} END{printf "%s\t%d\t%.6f\n",S,len,(L>0?100*len/L:0)}' \
  >> coverage_in_island.tsv
done

echo "[4/6] Cobertura total del plásmido (única) por sujeto…"
: > coverage_total.tsv
cut -f4 hsp_query.clean.bed | sort -u | while read -r sid; do
  awk -v S="$sid" '$4==S{print $1"\t"$2"\t"$3}' hsp_query.clean.bed \
  | bedtools sort -i - \
  | bedtools merge -i - \
  | awk -v S="$sid" -v L="$PLEN" 'BEGIN{len=0}{len+=($3-$2)} END{printf "%s\t%d\t%.6f\n",S,len,(L>0?100*len/L:0)}' \
  >> coverage_total.tsv
done

echo "[5/6] Integrar nt y % (denominador: plásmido)…"

# 1. Crear el archivo final escribiendo los headers primero
echo -e "sseqid\tnt_total\tnt_in_island\tnt_out_island\tpct_total\tpct_isla_plasmid\tpct_out_plasmid\tpct_isla_relative" > coverage_by_subject.tsv

# 2. Calcular los datos, ordenarlos y anexarlos (>>) debajo de los headers
awk -v PLEN="$PLEN" '
  NR==FNR { inis[$1]=$2; pir[$1]=$3; next }
  {
    sid=$1; nt_tot=$2; ptot=$3;
    nt_is = (sid in inis ? inis[sid] : 0);
    nt_out = nt_tot - nt_is; if (nt_out < 0) nt_out = 0;
    pisP = (PLEN>0 ? 100*nt_is/PLEN : 0);
    pout = (PLEN>0 ? 100*nt_out/PLEN : 0);
    pirv = (sid in pir ? pir[sid] : 0);
    printf "%s\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\n", sid, nt_tot, nt_is, nt_out, ptot, pisP, pout, pirv;
  }
' coverage_in_island.tsv coverage_total.tsv \
| sort -k7,7nr >> coverage_by_subject.tsv
