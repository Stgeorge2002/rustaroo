#!/bin/bash
# One node. Rustroo finishes, then Panaroo starts, so each tool gets
# all 32 CPUs and the whole 64 GB. Running them together would split both.
# Submit with --array=0-19. Sizes grow exponentially from 1 to 1500.
# 24 hours is the site task cap; the large Panaroo runs are the long ones.
#SBATCH --account=brics.b36ad
#SBATCH --partition=genoa
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --job-name=rustroo_curve
#SBATCH --output=logs/%x_%A_%a.out

set -euo pipefail

: "${SCRATCHDIR:?SCRATCHDIR is unset}"
MACS="${RUSTROO_MACS:-$SCRATCHDIR/rustaroo-macs}"
BIN="${RUSTROO_BIN:-$MACS/rustaroo/target/release/rustroo}"
SIF="${PANAROO_SIF:-$MACS/sif/panaroo-1.5.1.sif}"
INPUT="${RUSTROO_INPUT:-$MACS/monocle_gff}"
OUTPUT="${RUSTROO_OUTPUT:-$MACS/runs/curve_${SLURM_JOB_ID}}"
SCRIPT="${RUSTROO_CURVE:-$MACS/rustaroo/benchmarks/macs_curve.py}"
THREADS="${SLURM_CPUS_PER_TASK:?Submit with --cpus-per-task}"
QUALITY="${RUSTROO_QUALITY:-standard}"
# An array task is one sample size. Rustroo and Panaroo still share that node,
# one after the other, so each keeps all of its CPUs and memory.
TASK=""
N=""
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    TASK="$SLURM_ARRAY_TASK_ID"
    N="task"
else
    N="${RUSTROO_N:-}"
fi
if [[ -n "$N" ]]; then
    OUTPUT="${RUSTROO_OUTPUT:-$MACS/runs/curve_${QUALITY}}"
fi

# The job log is written outside the run directory. Copy it in when the
# task ends, including after a failed tool, so it downloads with the cohort.
copy_job_log() {
    mkdir -p "$OUTPUT/logs" || return 0
    if [[ -n "${SLURM_ARRAY_JOB_ID:-}" && -n "${SLURM_ARRAY_TASK_ID:-}" && -n "${SLURM_JOB_NAME:-}" ]]; then
        local src="logs/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out"
        if [[ -f "$src" ]]; then
            cp -f "$src" "$OUTPUT/logs/slurm_${SLURM_ARRAY_TASK_ID}.out" || true
        fi
    fi
}
trap copy_job_log EXIT

test "$(uname -m)" = x86_64
if [[ ! -x "$BIN" ]]; then
    echo "No rustroo binary at $BIN" >&2
    echo "From the rustaroo checkout: cargo build --release --locked" >&2
    exit 1
fi
if [[ "$QUALITY" == "high" ]] && ! "$BIN" --help 2>&1 | grep -q -- '--quality'; then
    echo "$BIN has no --quality flag. Rebuild the checkout: cargo build --release --locked" >&2
    exit 1
fi
test -s "$SIF"
test -d "$INPUT"
test -f "$SCRIPT"
command -v python3 >/dev/null
command -v apptainer >/dev/null

export APPTAINER_TMPDIR="$MACS/tmp"
export APPTAINER_CACHEDIR="$MACS/apptainer-cache"
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"

args=(
    python3 "$SCRIPT"
    --input "$INPUT"
    --output "$OUTPUT"
    --rustroo "$BIN"
    --panaroo-sif "$SIF"
    --threads "$THREADS"
    --quality "$QUALITY"
    --seed 1729
)
if [[ -n "$TASK" ]]; then
    args+=(--task "$TASK")
elif [[ -n "$N" ]]; then
    args+=(--n "$N")
fi
"${args[@]}"
