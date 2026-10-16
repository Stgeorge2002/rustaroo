#!/bin/bash
# One rustroo measurement of the uploaded monocle GFF directory.
# Submit from $SCRATCHDIR/rustaroo-macs so the relative log path lands on scratch:
#   sbatch --job-name=rustroo_monocle macs_monocle.sbatch
#SBATCH --account=brics.b36ad
#SBATCH --partition=genoa
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out

set -euo pipefail

: "${SCRATCHDIR:?SCRATCHDIR is unset}"
MACS="${RUSTROO_MACS:-$SCRATCHDIR/rustaroo-macs}"
INPUT="${RUSTROO_INPUT:-$MACS/monocle_gff}"
OUTPUT="${RUSTROO_OUTPUT:-$MACS/runs/monocle_${SLURM_JOB_ID}}"
THREADS="${SLURM_CPUS_PER_TASK:?Submit with --cpus-per-task}"

if [ -n "${RUSTROO_SIF:-}" ]; then
    IMAGE=$RUSTROO_SIF
else
    mapfile -t images < <(find "$MACS/sif" -maxdepth 1 -type f -name 'rustroo-*-linux-amd64.sif' | sort)
    if [ "${#images[@]}" -ne 1 ]; then
        printf 'rustroo: expected one image in %s/sif, found %s\n' "$MACS" "${#images[@]}" >&2
        exit 1
    fi
    IMAGE=${images[0]}
fi

test "$(uname -m)" = x86_64
test -s "$IMAGE"
test -d "$INPUT"
command -v apptainer >/dev/null

if [ -f "${IMAGE}.sha256" ]; then
    (
        cd "$(dirname "$IMAGE")"
        sha256sum --check "$(basename "$IMAGE").sha256"
    )
fi

mapfile -t gffs < <(find "$INPUT" -maxdepth 1 -type f \( -name '*.gff' -o -name '*.gff3' -o -name '*.gff.gz' -o -name '*.gff3.gz' \) | sort)
if [ "${#gffs[@]}" -eq 0 ]; then
    printf 'rustroo: no GFF files in %s\n' "$INPUT" >&2
    exit 1
fi

mkdir -p "$OUTPUT" "$MACS/tmp" "$MACS/apptainer-cache"
printf '%s\n' "${gffs[@]##*/}" > "$OUTPUT/inputs.txt"
{
    printf 'job_id\t%s\n' "$SLURM_JOB_ID"
    printf 'node\t%s\n' "$(hostname)"
    printf 'arch\t%s\n' "$(uname -m)"
    printf 'cpus\t%s\n' "$THREADS"
    printf 'mem\t%s\n' "${SLURM_MEM_PER_NODE:-}"
    printf 'image\t%s\n' "$IMAGE"
    printf 'input\t%s\n' "$INPUT"
    printf 'gff_count\t%s\n' "${#gffs[@]}"
    date -Is
} > "$OUTPUT/job.txt"

export APPTAINER_TMPDIR="$MACS/tmp"
export APPTAINER_CACHEDIR="$MACS/apptainer-cache"

apptainer exec --cleanenv --no-home \
    --bind "$INPUT:/input:ro" \
    --bind "$OUTPUT:/output" \
    "$IMAGE" \
    env RAYON_NUM_THREADS="$THREADS" \
    rustroo -i /input -o /output \
        --threads "$THREADS" --audit-members
