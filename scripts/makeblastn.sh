#!/bin/bash
# ==============================================================================
# Script: makeblastn.sh
# Descripción: Ejecuta búsquedas locales con blastn contra bases de datos NCBI (nt)
#              generando salidas tabulares delimitadas por tabulador (outfmt 6 o 7).
# Uso: ./makeblastn.sh <database_path> <query_fasta> <outfmt_code>
# Ejemplo: ./makeblastn.sh /export/storage/data/blast/db/nt ../assembly_pilon.ffn 6
# ==============================================================================

set -euo pipefail

if [ "$#" -lt 3 ]; then
    echo "Error: Argumentos insuficientes." >&2
    echo "Uso: $0 <database_path> <query_fasta> <outfmt_format>" >&2
    echo "  outfmt_format: 6 (tabular sin headers), 7 (tabular con comentarios)" >&2
    exit 1
fi

# Argumentos posicionales
db="$1"
fasta="$2"
format="$3"

# Extraer el nombre base del archivo de secuencias
file=$(basename "$fasta")
name="${file%.*}"

# Ejecución de blastn
blastn \
  -db "$db" \
  -query "$fasta" \
  -outfmt "$format" \
  -max_target_seqs 1 \
  -out "${name}.blastn_m${format}.txt"

echo "Búsqueda BLASTn completada: ${name}.blastn_m${format}.txt"
