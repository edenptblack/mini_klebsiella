# 02_mlst.sh - Assign a sequence type to each downloaded genome

# Run from project root directory with: bash scripts/02_mlst.sh

# Stop the pipeline if any command fails or if an unset variable is used
set -euo pipefail

# --- Configuration settings ---
GENOME_DIR="data/genomes"
OUT_DIR="results/mlst"
ENV="kp-typing"
SCHEME="klebsiella"

# --- Run tool ---
mkdir -p "$OUT_DIR"
echo "mlst is typing genomes in $GENOME_DIR with scheme $SCHEME"

conda run -n "$ENV" mlst --scheme $SCHEME "$GENOME_DIR"/*.fna > "$OUT_DIR/mlst.tsv"

echo "MLST complete -> $OUT_DIR/mlst.tsv"