#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

BIN=./target/release/rustroo
OUTBASE=bench_out
rm -rf "$OUTBASE"
mkdir -p "$OUTBASE"
RESULTS="$OUTBASE/results.tsv"
echo -e "n_samples\tn_files\tseconds\tgenes\tsamples_reported\tbroken" > "$RESULTS"

mapfile -t CHROM_FILES < <(ls gff_test/chromosome/*.gff3 | sort)

for n in 3 4 5 6 7 8 9 10; do
  FILES=()
  for ((i = 0; i < n; i++)); do
    f="${CHROM_FILES[$i]}"
    FILES+=("$f")
    base=$(basename "$f" .gff3)
    for p in gff_test/plasmid/"${base}"__plasmid_*.gff3; do
      [ -e "$p" ] && FILES+=("$p")
    done
  done

  OUTDIR="$OUTBASE/n$n"
  mkdir -p "$OUTDIR"

  START=$(date +%s.%N)
  LOG=$("$BIN" -i "${FILES[@]}" -o "$OUTDIR" 2>&1)
  END=$(date +%s.%N)
  ELAPSED=$(echo "$END - $START" | bc)

  GENES=$(echo "$LOG" | grep -m1 '^genes' | cut -f2)
  SAMP=$(echo "$LOG" | grep -m1 '^samples' | cut -f2)
  BROKEN=$(echo "$LOG" | grep -m1 'had no intact copy' | grep -oE '^[0-9]+' || echo 0)

  echo -e "${n}\t${#FILES[@]}\t${ELAPSED}\t${GENES}\t${SAMP}\t${BROKEN}" >> "$RESULTS"
  echo "n=${n} files=${#FILES[@]} time=${ELAPSED}s genes=${GENES} samples=${SAMP}" >&2
done

echo "---"
column -t "$RESULTS"
