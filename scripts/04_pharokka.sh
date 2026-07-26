# 04_pharokka.sh - Annotate phage regions identified by genomad

# Run from project root directory with: bash scripts/04_pharokka.sh

# Stop the pipeline if any command fails or if an unset variable is used
set -euo pipefail

# --- Configuration settings ---
GENOMAD_DIR="results/genomad"
OUT_DIR="results/pharokka"
DB="databases/pharokka_db"
ENV="kp-pharokka"
THREADS=2

# --- Run tool ---
mkdir -p "$OUT_DIR"

for sample_dir in $GENOMAD_DIR/*/; do

    sample=$(basename "$sample_dir")

    prophages=$sample_dir/${sample}_summary/${sample}_virus.fna

    if [[ ! -s "$prophages" ]]; then
        echo "[pharokka]: No prophages found in $sample, skipping"
        continue
    fi

    echo "[pharokka] $sample"

    # -m: meta mode. Use all available threads
    # -f: Force. Ensure repeat work on rerun
    conda run --live-stream -n "$ENV" \
        pharokka.py \
            -i "$prophages" \
            -o "$OUT_DIR/$sample" \
            -d "$DB" \
            -t "$THREADS" \
            -m \
            -f

done

echo "[pharokka]: All complete -> $OUT_DIR"