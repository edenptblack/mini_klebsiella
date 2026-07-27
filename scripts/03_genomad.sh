# 03_genomad.sh - Find prophage regions in each genome

# Run from project root directory with: bash scripts/03_genomad.sh

# Stop the pipeline if any command fails or if an unset variable is used
set -euo pipefail

# --- Configuration settings ---
GENOME_DIR="data/genomes"
ACCESSIONS="data/accession_selected.txt"
OUT_DIR="results/genomad"
DB="databases/genomad_db"
ENV="kp-genomad"
THREADS=2
SPLITS=8

# --- Run tool ---
mkdir -p "$OUT_DIR"

for sample in $(cat "$ACCESSIONS"); do
    
    genome="$GENOME_DIR/$sample.fna"

    if [[ ! -s "$genome" ]]; then
        echo "[genomad]: $genome not found, skipping"
        continue
    fi

    if [[ -d "$OUT_DIR/$sample/${sample}_summary" ]]; then
        echo "[genomad] $sample already done, skipping"
        continue
    fi

    echo "[genomad] $sample"

    conda run --live-stream -n "$ENV" \
        genomad end-to-end \
            --threads "$THREADS" \
            --splits "$SPLITS" \
            --disable-nn-classification \         # NN computationally expensive and doesn't seem to be used for prophage identification
            "$genome" "$OUT_DIR/$sample" "$DB"
done

echo "[genomad]: All complete -> $OUT_DIR/"