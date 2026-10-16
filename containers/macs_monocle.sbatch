#!/bin/bash
# One allocation. Rustroo finishes, then Panaroo starts, so each tool gets
# all 32 CPUs and the whole 64 GB. Running them together would split both.
# 24 hours is the site task cap; the 800- and full-set Panaroo runs are the long ones.
#SBATCH --account=brics.b36ad
#SBATCH --partition=genoa
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out

set -euo pipefail

: "${SCRATCHDIR:?SCRATCHDIR is unset}"
MACS="${RUSTROO_MACS:-$SCRATCHDIR/rustaroo-macs}"
BIN="${RUSTROO_BIN:-$MACS/rustroo}"
SIF="${PANAROO_SIF:-$MACS/sif/panaroo-1.5.1.sif}"
INPUT="${RUSTROO_INPUT:-$MACS/monocle_gff}"
OUTPUT="${RUSTROO_OUTPUT:-$MACS/runs/curve_${SLURM_JOB_ID}}"
SCRIPT="${RUSTROO_CURVE:-$MACS/rustaroo/benchmarks/macs_curve.py}"
THREADS="${SLURM_CPUS_PER_TASK:?Submit with --cpus-per-task}"

test "$(uname -m)" = x86_64
test -x "$BIN"
test -s "$SIF"
test -d "$INPUT"
test -f "$SCRIPT"
command -v python3 >/dev/null
command -v apptainer >/dev/null

export APPTAINER_TMPDIR="$MACS/tmp"
export APPTAINER_CACHEDIR="$MACS/apptainer-cache"
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"

python3 "$SCRIPT" \
    --input "$INPUT" \
    --output "$OUTPUT" \
    --rustroo "$BIN" \
    --panaroo-sif "$SIF" \
    --threads "$THREADS" \
    --seed 1729
