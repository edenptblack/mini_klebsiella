# 04_defensefinder - Find defense and anti-defense systems

# NB: This analysis is run on whole-genome proteins, not extracted prophage proteins.
# Defence systems are often multi-gene operons, and defensefinder locates these using genomic co-location. 
# Hits will be inspected for intersection which prophages in downstream analysis
# geNomad already produces whole-genome gene calling, so we can use that

# Run from project root directory with: bash scripts/05_defensefinder.sh

set -euo pipefail

# --- Configuration settings ---
GENOMAD_DIR="results/genomad"
OUT_DIR="results/defensefinder"
ENV="kp-defensefinder"

# --- Run tool ---
mkdir -p "$OUT_DIR"

for sample_dir in "$GENOMAD_DIR"/*/; do

    sample=$(basename "$sample_dir")

    proteins="$sample_dir/${sample}_annotate/${sample}_proteins.faa"

    if [[ ! -s "$proteins" ]]; then
        echo "[defensefinder]: no $sample protein file found, skipping"
        continue
    fi

    echo "[defensefinder] $sample"

    # antidefensefinder locates anti-defense elements.
    # preserve-raw keeps underlying HMMER hit information for later inspection if necessary
    conda run --live-stream -n "$ENV" \
        defense-finder run \
            -o "$OUT_DIR/$sample" \
            --antidefensefinder \
            --preserve-raw \
            "$proteins"

done

echo "[defensefinder]: All complete -> $OUT_DIR"