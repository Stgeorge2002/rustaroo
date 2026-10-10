# Chromosome-only rustroo vs Panaroo benchmark

This protocol measures computational scaling and reference-sequence concordance
on the existing chromosome GFF3 collection. It does not establish biological
ground-truth accuracy. Neither tool is a gold standard.

For the original paper's simulations, real-data controls and scaling tests,
see [Panaroo's published testing methods](PANAROO_PUBLISHED_TESTING.md).
That guide distinguishes verified publication details from reproduction gaps
and from this repository's modern paired benchmark.

## Eight-pair timing-only run

For one run per tool at **5, 10, 15, 20, 25, 30, 35, and 40 genomes**,
with no warmups or reference comparisons, use:

```bash
cd ~/rustaroo
RUN="$PWD/benchmark_8pairs_$(date +%Y%m%d_%H%M%S)"
python3 benchmarks/benchmark.py prepare \
  --input gff_test/chromosome \
  --binary target/release/rustroo \
  --output "$RUN" \
  --image rustroo-panaroo-benchmark:1.5.1-time \
  --threads 26 \
  --sizes 5 10 15 20 25 30 35 40 \
  --subsets 1 \
  --repeats 1 \
  --no-warmups \
  --no-concordance \
  --seed 1729 \
  --clean-mode strict && \
python3 "$RUN/protocol/benchmark.py" run --output "$RUN"
```

This performs **16 tool executions total**, with one seeded nested subset series.
Preparation reports `16 measured runs, 0 warmups, 0 comparisons`; BLAST is not
required or invoked. Runtime/resource tables, speedups, counts and SVG figures
are still generated, but concordance artifacts are not.
With one run per point there is **no estimate of timing variability or
statistical uncertainty**; the median/min/max columns all equal the single
measurement. This is a descriptive scaling run, not the replicated publication
protocol below.

If an older, larger experiment is still running, interrupt it in its terminal
with Ctrl+C and verify it has stopped before starting this one. Do not modify
its saved protocol or configuration: use a new output directory.

## Design

- Sample counts: **3, 5, 10, 15, 20, 25, 30, 35, 40**.
- Five random subset series, seeds **1729--1733**, sampled without replacement.
  Each series is nested: its 3 genomes are included in its 5, and so on.
  Both tools receive the same sorted files at every point.
- Five measured repetitions per tool and subset: **450 measured executions**.
- One excluded warmup per tool and subset: **90 additional executions**.
- Subset/size blocks and the tool order within each paired repetition are
  randomized with a saved schedule. Tools run sequentially, never concurrently.
- Reference concordance uses repetition 1 for each of the 45 subsets, after all
  performance runs. BLAST time is excluded from performance measurements.

These are five independently drawn subset series from **one 44-genome dataset**,
not five independent biological cohorts. Subsets overlap, especially at 40
genomes, and adjacent sizes within a series are dependent.

## Annotation-level family audit

Reference-only reciprocal best hits cannot explain excess families: they lose
many-to-one relationships, ties, and the original gene memberships. The separate
`audit_families.py` workflow compares annotations without altering either tool's
clustering. It requires Rustroo's new `--audit-members` export and Panaroo's
`gene_presence_absence.csv`, `gene_data.csv`, and `pan_genome_reference.fa`.
Older Rustroo outputs cannot supply the lost original memberships and must be
rerun with that flag. Use the **same input GFF files**, not just the same sample
count. The matrices' sample sets must match exactly; do not mix 40- and
44-chromosome outputs or chromosome and plasmid cohorts.

### Generate membership evidence

Validation and build, for user execution:

```bash
cd ~/rustaroo
cargo test --lib && \
python3 -m unittest discover -s benchmarks -p 'test_*.py' && \
cargo build --release --locked
```

After successful validation:

```bash
set -o pipefail
THREADS=26
RUN="$PWD/full44_membership_audit_$(date +%Y%m%d_%H%M%S)"
PAN="$PWD/chromosome_compare_20261009_140422/panaroo"
mkdir -p "$RUN"

RAYON_NUM_THREADS="$THREADS" /usr/bin/time \
  -f 'Wall seconds: %e\nUser seconds: %U\nSystem seconds: %S\nMaximum-process RSS (KiB): %M' \
  -o "$RUN/rustroo_time.txt" \
  ./target/release/rustroo \
    -i gff_test/chromosome -o "$RUN/rustroo" \
    --threads "$THREADS" --audit-members \
  2>&1 | tee "$RUN/rustroo.log"
```

The existing full-44 Panaroo output above can be reused **if the GFF inputs have
not changed**. Membership-enabled timing includes extra export/storage work;
do not substitute it for the earlier timing-only benchmark.

### Sequence comparison

Run this only after the Rustroo command succeeds. Use Rustroo references as
queries against **all** Panaroo references. Do not request only one target:
that would hide repeated genes and possible partitions. This command runs BLAST
only, not Panaroo:

```bash
(
  set -euo pipefail
  : "${RUN:?Set RUN to the completed membership-enabled Rustroo run}"
  : "${PAN:?Set PAN to the matching Panaroo output directory}"
  THREADS="${THREADS:-26}"
  IMAGE=staphb/panaroo:1.5.1
  IMAGE_ID=$(docker image inspect --format '{{.Id}}' "$IMAGE")
  for file in "$RUN/rustroo/genes.fa" "$RUN/rustroo/members.tsv" \
              "$PAN/pan_genome_reference.fa"; do
    if [ ! -s "$file" ]; then
      printf 'Missing or empty input: %s\n' "$file" >&2
      exit 1
    fi
  done
  STAMP=$(date +%Y%m%d_%H%M%S)
  SEQ="$RUN/sequence_$STAMP"
  REPORT="$RUN/report_$STAMP"
  mkdir "$SEQ"
  printf '%s\n' "$IMAGE_ID" > "$SEQ/image_id.txt"

  if ! docker run --rm --pull=never --network=none \
    --user "$(id -u):$(id -g)" \
    -v "$RUN/rustroo:/rustroo:ro" \
    -v "$PAN:/panaroo:ro" \
    -v "$SEQ:/work" \
    "$IMAGE_ID" bash -euc '
      blastn -version
      makeblastdb -in /panaroo/pan_genome_reference.fa -dbtype nucl -out /work/panaroo_db
      targets=$(grep -c "^>" /panaroo/pan_genome_reference.fa)
      blastn -task blastn -query /rustroo/genes.fa -db /work/panaroo_db \
        -num_threads "$1" -evalue 1e-10 -max_target_seqs "$targets" -max_hsps 1 \
        -outfmt "6 qseqid sseqid pident length nident qlen slen qstart qend sstart send evalue bitscore" \
        -out /work/rustroo_to_panaroo.tsv
    ' bash "$THREADS" > "$SEQ/blast.log" 2>&1; then
    printf 'BLAST failed; diagnostic log follows:\n' >&2
    cat "$SEQ/blast.log" >&2
    exit 1
  fi
  python3 benchmarks/audit_families.py \
    --rustroo "$RUN/rustroo" --panaroo "$PAN" \
    --hits "$SEQ/rustroo_to_panaroo.tsv" \
    --output "$REPORT"
)
```

This uses the original Panaroo image: GNU time is not needed for BLAST.
Image inspection must succeed before any container or audit starts; the
container uses that inspected image ID. The fail-fast subshell prevents an
image/BLAST failure from falling through to the audit without changing the
interactive shell's error settings. Each attempt has fresh sequence/report
directories, preserving any interrupted attempt instead of trusting partial
hits. This step does not pull or build an image automatically.

If BLAST fails, inspect the printed log or that attempt's
`sequence_<timestamp>/blast.log`; do not audit a partial hit file.
The audit can also run without `--hits` for a membership-first report, but then
only identical references are sequence-confirmed and other pairs are explicitly
`not_aligned`, not failures. Sequence thresholds default to >=70% nucleotide
identity and >=90% coordinate coverage of **both** references. These report
thresholds do not modify homology parameters.

The audit prints flushed stage messages, periodic Panaroo row counts, and
hashing/rechecking progress to stderr. Reading the large `gene_data.csv` and
verifying its hash may take time; a quiet interval is not proof of a hang.
Identical Panaroo DNA strings are pooled rather than retained once per
occurrence, and comparisons still use exact strings, not hash-only equality.
Missing input files and malformed BLAST hits fail before loading the large
Panaroo membership table. Ctrl+C returns status 130 with an explicit
incomplete-audit message. Restart just this comparison step; no pangenome
rerun is required if the membership-enabled output is already complete.

### Outputs and interpretation

| File | Evidence |
|------|----------|
| `summary.json` | Counts, exact group-delta decomposition, exceptions, thresholds, input/report/driver hashes |
| `components.tsv` | Connected membership components: one-to-one, Rustroo partitions, Rustroo unions, many-to-many, or unlinked families |
| `families.tsv` | Every family, component, source coverage, missing/refound/fusion counts, reference length and tied best sequence targets |
| `membership_pairs.tsv` | Every shared-source family pair: annotation overlap, presence Jaccard, exact member DNA checks, reference identity and both coverages |
| `sequence_hits.tsv` | Best single HSP for every reported reference pair, including tied and below-threshold hits |
| `annotation_exceptions.tsv` | Sample-scoped IDs that are unshared, ambiguous, fusion-derived, refound, missing, or sequence/scaffold-discordant |

A Rustroo partition of one Panaroo family is **not automatically an error**.
Co-occurring Rustroo groups are highlighted so they are not merged merely to
reduce the count. A sequence-supported partition with disjoint sample sets is
a useful case for annotation/neighborhood review, not an automatic union.
Components with source sequence/scaffold disagreements are flagged first:
matching IDs do not prove that both tools read identical inputs or interpreted
the CDS in the same way.

The graph uses `(sample, original annotation ID)`, never bare IDs (the Bakta
`KPN_...` identifiers repeat across samples). Refound occurrences and missing
IDs are excluded from membership edges. Fusion or nonunique source assignments
are excluded on **both** sides to avoid bridging unrelated families; they stay
visible as exceptions. Panaroo refound tokens follow its `_refound_` naming
convention. Unknown non-refound matrix IDs fail explicitly rather than being
matched by name similarity. Families with no unambiguous edge are labeled
unlinked even if their references resemble each other.

Presence Jaccard measures sample overlap, not copy-level membership agreement.
Membership fractions use only unambiguous original source annotations.
Excluded source counts are reported separately; even a Jaccard of 1 does not
certify complete agreement when either family has excluded occurrences.
Sequence coverage is from one HSP, not the sum of potentially overlapping HSPs.
All target ties are retained; no reciprocal-best-hit filter discards excess
families. A missing Panaroo reference is distinct from a missing Panaroo
family. No reported hit is not proof that homology is absent. The audit retains
failed/unavailable evidence explicitly and refuses to overwrite a report.
Neither tool is biological ground truth, and this audit cannot assign causes
such as "invalid CDS filtering" solely from absence in an output.

## Running in WSL/Linux

Use a Linux filesystem, not a Windows-mounted drive. Python 3.9+ and its standard
library, Docker, GNU `/usr/bin/time`, `taskset`, `lscpu`, Cargo and rustc are
required. The existing image must contain GNU `/usr/bin/time`, Panaroo, BLASTN
and makeblastdb. Preparation checks this; it does not install or pull anything.
If GNU time is absent from the image, build and document a derived image adding
it, then provide that image with `--image`. The provided Dockerfile adds only
the GNU `time` package to the Debian/Ubuntu-based Panaroo image:

```bash
docker build --pull=false \
  -f benchmarks/panaroo.Dockerfile \
  -t rustroo-panaroo-benchmark:1.5.1-time \
  benchmarks
```

The build needs network access to the distribution package repositories; it
does not run a pangenome analysis. Use
`--image rustroo-panaroo-benchmark:1.5.1-time` during preparation below.
The base image and distribution packages should be archived for publication;
repository package versions can change between builds. The resulting image ID,
Dockerfile, and tool/timer versions are recorded by preparation.
When the original tag has a repository digest, pin it at build time with
`--build-arg PANAROO_BASE=staphb/panaroo@sha256:<repository-digest>`.
Use a repository digest from `docker image inspect`, not a local image ID.
No package installation or fallback timer happens automatically in the harness.

Build separately from the benchmark:

```bash
cd ~/rustaroo
cargo build --release --locked

RUN="$PWD/benchmark_publication_$(date +%Y%m%d_%H%M%S)"
python3 benchmarks/benchmark.py prepare \
  --input gff_test/chromosome \
  --binary target/release/rustroo \
  --output "$RUN" \
  --image staphb/panaroo:1.5.1 \
  --threads 26 \
  --subsets 5 \
  --repeats 5 \
  --seed 1729 \
  --clean-mode strict && \
python3 "$RUN/protocol/benchmark.py" run --output "$RUN"
```

**Preparation does not run the pangenome pipelines.** It checks versions,
copies inputs and the executable, and records the protocol. Only `run` launches
the actual benchmark. Allow several hours and substantial disk space: all
warmup and repetition outputs are retained. No existing directory is deleted.
Do not start `run` if preparation fails: the saved runner is created only after
the dependency preflight succeeds. Preflight failures now include the actual
command output as well as the exit status.

For separate performance and concordance stages:

```bash
python3 "$RUN/protocol/benchmark.py" run --output "$RUN" --stage performance
python3 "$RUN/protocol/benchmark.py" run --output "$RUN" --stage concordance
python3 "$RUN/protocol/benchmark.py" report --output "$RUN"
```

Rerunning skips verified successful performance jobs and successful comparisons.
Failures and interruptions stop execution, retain logs, and are never replaced
by zero timings. Inspect the named failed job directory, ensure its process or
named container has stopped, then move **that exact directory** outside `jobs`
or `concordance` before explicitly retrying. A still-running attempt is not safe
to retry. Reports require every planned execution and comparison to succeed.

The saved executable, input hashes, driver, image ID, CPU set, or protocol must
not be changed mid-experiment. Use a new output directory for another build or
cleaning mode. Preparation snapshots Rust sources and Cargo manifests, but a
binary hash does not prove those sources produced the binary: build immediately
before preparation and archive the build/toolchain information.

## Timing and resource definitions

- Primary runtime: GNU time's wall time around the **tool invocation**, including
  input parsing, analysis, output writes, and tool-internal subprocesses.
  Container startup/teardown, compilation, warmups, and BLAST are excluded.
- `launcher_wall_seconds` separately records host-observed invocation time,
  including Docker startup/teardown for Panaroo.
- User/system CPU times are GNU time measurements inside the Panaroo container
  and on the native host for rustroo, not measurements of the Docker client.
- `average_cpu_cores = (user_seconds + system_seconds) / wall_seconds`.
- `max_process_rss_kib` is GNU time `%M`: the maximum-process resident-set
  measurement, **not the simultaneous sum across the process tree**. Panaroo
  launches subprocesses, so this can substantially underestimate aggregate
  pipeline memory. Do not publish this as total peak RAM or infer total-memory
  superiority. A total-memory comparison needs a separate, consistently
  isolated cgroup measurement for both tools.
- Native rustroo uses `taskset`; Docker uses `--cpuset-cpus` with the same logical
  CPU IDs. Both have the same thread budget. Rayon and OpenMP are configured;
  unrelated BLAS thread pools are limited to one thread. This does not mean all
  pipeline stages can use the full budget.

This deliberately preserves the current native-rustroo/container-Panaroo setup.
It is not an identical-runtime-environment experiment. Docker Desktop/WSL may
use a separate daemon kernel/VM; CPU-ID equivalence alone does not establish
physical placement equivalence. Record Docker/WSL versions, host hardware,
Docker Desktop CPU/RAM limits, storage, power mode, and any competing workloads
in the manuscript. The harness captures available software versions, logical
CPU topology, affinity, host memory, and per-run host load.

Use an otherwise idle machine on AC power with a stable power policy. There
is no global cache flushing, governor modification, or process killing. This
is a warmed, uncontrolled-OS-cache protocol, not a cold-cache experiment.
Successful warmups are skipped on resume. After a long interruption, those
warmups may no longer represent the cache state. For publication, restart the
affected subset block including its warmups and all five repetitions, or use
a fresh experiment; do not silently pool measurements across changed machine
conditions. The execution lock prevents overlapping stages in one output
directory, not interference from unrelated workloads or other experiments.

## Tool settings and concordance

Panaroo defaults to the previous run's `--clean-mode strict`,
`--remove-invalid-genes`, and the selected thread count. rustroo explicitly uses
0.98/0.70 clustering identities, 0.95 length ratio, minimum CDS 90 bp, and refind
enabled. The algorithms and effective CDS inclusion rules differ. Invalid-CDS
log messages, reference counts, group counts, and core/soft-core group counts
are retained. Strict cleaning is not a matched analogue of rustroo's behavior.
For a cleaning-mode sensitivity analysis, repeat in a separate directory with
`--clean-mode sensitive`; do not combine the results with strict runs.

BLASTN runs in both directions with `-task blastn`, `-dust no`, and e-value
`1e-10`. The target limit equals the database's reference count, avoiding the
previous comparison's 100-target truncation. The highest-bit-score single HSP is
used. Queries tied across different subjects are marked ambiguous and excluded.
Unique reciprocal best hits are accepted with **>=70% nucleotide identity and
>=90% coordinate-span coverage of both sequences**. Coverage includes the span
of the aligned sequence, not a stitched collection of disjoint HSPs.

Per-subset results include accepted matches, exact full-length matches, mean
and alignment-column-weighted identity, no-hit sequences, ambiguous hits,
nonreciprocal unique hits, and reciprocal pairs rejected by the filters. These
categories distinguish unmatched references from demonstrably absent genes.
The accepted subset is selected for similarity: high identity in that subset
must not be described as whole-tool accuracy. This protocol does not evaluate
ortholog membership, paralog resolution, or sample-level presence accuracy.

Repeat-output hashes detect reference or matrix variability. Reference hashes
normalize FASTA record order; matrix hashes are byte-level and can change solely
because of ordering. Concordance uses repetition 1, not the most favorable run.
Investigate variability before drawing reproducibility claims.

## Publication outputs

The `report/` directory contains:

- `runs.tsv`: all 450 measured executions, no warmups.
- `subset_summary.tsv`: five-repeat medians and ranges per subset/tool.
- `scaling_summary.tsv`: medians and ranges across the five subset medians.
- `paired_speedups.tsv`, `speedup_summary.tsv`: per-subset ratios of median
  Panaroo runtime to median rustroo runtime, then their cross-subset summaries.
- `concordance.tsv`, `concordance_summary.tsv`: all 45 subset comparisons and
  summaries at each size. Undefined identity for zero matches is blank, not 0%.
- `output_variability.json`: output-hash differences between timing repetitions.
- Editable SVG figures for runtime, maximum-process RSS, group counts, speedup,
  and matched-reference fractions.

Figure error bars show the **minimum--maximum across five subset medians**,
not 95% confidence intervals. Timing repetitions are not pooled as independent
biological replicates. Use the raw tables for other justified analyses; do not
treat the 25 timings per point as 25 independent sampled cohorts.

Archive `config.json`, `environment.json`, `manifests/`, `protocol/`, logs,
raw results, input checksums, and the image digest with the paper. Cite the
data's origin and licensing; this harness does not supply those. A single
species/dataset and 3--40 closely related genomes do not demonstrate general
scaling or accuracy across diverse bacterial populations.

## Harness tests

These use synthetic files and mocked command construction only. They never run
Panaroo, rustroo, Docker, or the real GFF dataset:

```bash
python3 -m unittest discover -s benchmarks -p 'test_*.py' -v
```
