# Standalone Rustroo Apptainer image for MACS

This packages Rustroo separately from the shared pipeline tool image.
The binary is `rustroo` even though the repository is named `rustaroo`.
No image build, Cargo tests/builds, container execution, or MACS submission
was performed when adding this recipe.

## What is packaged

[rustroo.def](rustroo.def) uses two digest-pinned Linux **amd64** stages:

| Stage | Originating image tag | Pinned platform manifest |
|---|---|---|
| Builder | `rust:1.90.0-bookworm` | `sha256:3f6e6f8d8725a65a2db964bb828850f888d430c68784d661f753144e5d787207` |
| Runtime | `debian:bookworm-slim` | `sha256:a4672c0cb26fbdde88e38fa2dfb6c681942306680e41e4378b28770b6e79ee91` |

The builder runs the existing library and CLI tests in release mode, then
compiles the executable with `Cargo.lock` enforced by `--locked`.
The runtime contains the stripped executable, build metadata and Debian's
runtime libraries. It does not contain Cargo, the Rust toolchain, source,
datasets, or benchmark outputs. Runtime analysis requires no network access
or external bioinformatics executables.

The CPU target is baseline `x86-64`, not `target-cpu=native`, so the build does
not assume the build machine's optional CPU instructions. Both stages use
Debian Bookworm to keep the glibc runtime compatible.

**Architecture check:** run `uname -m` on the actual target compute node.
This recipe is for `x86_64`, following the requested MACS target. If the node
reports `aarch64`, stop: this image is not suitable. An ARM image needs its own
platform digests, target and validation; this recipe does not cross-compile.

## Build once on Linux x86-64

Use WSL Ubuntu or another approved Linux x86-64 build machine with Apptainer
1.2 or newer, networking, and configured fakeroot support. The initial build
downloads the builder/runtime images and locked crates; its temporary disk
requirements are much larger than the final SIF. Do not build on a cluster
login node unless site policy explicitly permits it.

First commit the intended source **and these packaging files**, reviewing the
staged contents rather than using `git add .`. Then select that commit.
The following commands are for you to run; they have not been executed.

```bash
(
  set -euo pipefail
  cd ~/newdev/rustaroo

  test "$(uname -m)" = x86_64
  apptainer --version
  if ! git diff --quiet HEAD -- Cargo.toml Cargo.lock src; then
    echo "Rustroo source differs from HEAD; review and commit it before building." >&2
    git --no-pager diff --stat HEAD -- Cargo.toml Cargo.lock src
    exit 1
  fi
  REV=$(git rev-parse --verify HEAD)
  git cat-file -e "$REV:containers/rustroo.def"

  mkdir -p "$PWD/out/containers"
  IMAGE="$PWD/out/containers/rustroo-${REV:0:12}-linux-amd64.sif"
  test ! -e "$IMAGE"
  SNAPSHOT=$(mktemp -d "${TMPDIR:-/tmp}/rustroo-apptainer.XXXXXXXX")
  printf 'Building commit: %s\nSource snapshot: %s\nImage: %s\n' \
    "$REV" "$SNAPSHOT" "$IMAGE"

  git archive --format=tar "$REV" | tar -xf - -C "$SNAPSHOT"
  cd "$SNAPSHOT"

  apptainer build --fakeroot \
    --build-arg "SOURCE_REVISION=$REV" \
    "$IMAGE" containers/rustroo.def \
    2>&1 | tee "$IMAGE.build.log"

  apptainer test "$IMAGE"
  apptainer exec --cleanenv "$IMAGE" rustroo --version
  apptainer exec --cleanenv "$IMAGE" rustroo validate --help
  apptainer exec --cleanenv "$IMAGE" \
    cat /opt/rustroo-build/source-revision.txt

  cd "$(dirname "$IMAGE")"
  sha256sum "$(basename "$IMAGE")" > "$(basename "$IMAGE").sha256"
  printf 'Ready to transfer: %s\n' "$IMAGE"
  printf 'Retained source snapshot: %s\n' "$SNAPSHOT"
)
```

The source comes from `git archive`, not the potentially dirty working tree.
The full revision is recorded inside the image; manifest/source/binary SHA-256
checksums are also included under `/opt/rustroo-build`.
These pins identify inputs, but do not promise byte-identical SIF files across
different build hosts or build times.

The repository's `.gitattributes` normalizes text to LF so Windows and WSL Git
agree about source contents. If an older checkout shows whole files as changed,
compare with `git diff --ignore-cr-at-eol --stat -- Cargo.toml Cargo.lock src`
to distinguish line endings from real edits. Keep the source guard: review and
commit real changes rather than removing it to build an outdated revision.

Keep the snapshot until the build is verified, then remove only that specific
temporary directory if desired. If fakeroot is unavailable, use the site's
approved build service or a privileged build on a machine you administer;
do not bypass MACS permissions. If tests or compilation fail, retain the build
log and fix the failure rather than using `--notest` or changing the lockfile
silently.

To audit a recipe later, use `apptainer inspect --deffile IMAGE.sif` and inspect
`/opt/rustroo-build/toolchain.txt` and `source.sha256` within the image.

## Transfer to MACS

Transfer the SIF and its matching `.sha256` file to the pipeline's existing
container directory. Use your actual SSH alias and directory:

```bash
scp out/containers/rustroo-<commit>-linux-amd64.sif \
    out/containers/rustroo-<commit>-linux-amd64.sif.sha256 \
    YOUR_MACS_SSH_ALIAS:/YOUR/PIPELINE/CONTAINER_DIRECTORY/
```

On MACS, in that directory:

```bash
sha256sum --check rustroo-<commit>-linux-amd64.sif.sha256
```

Use a versioned, immutable filename in the workflow configuration. Do not
replace a SIF underneath running or resumed Nextflow jobs: configure a new
filename for each new build. No private repository credentials need to be
embedded in the image.

## Standalone execution

This example runs the existing 90 whole-isolate Bakta GFFs. It does not also add
their overlapping plasmid-bin annotations. Edit paths to match MACS storage.
Run analysis through the site's normal scheduler allocation.

```bash
(
  set -euo pipefail
  IMAGE="/YOUR/PIPELINE/CONTAINER_DIRECTORY/rustroo-<commit>-linux-amd64.sif"
  INPUT="/YOUR/DATA/macs_PRJNA1078342_full_20261012_132112_selected_results"
  OUTPUT="/YOUR/RESULTS/rustroo_$(date +%Y%m%d_%H%M%S)"
  THREADS="${SLURM_CPUS_PER_TASK:?Run inside an allocated task with cpus-per-task set}"
  mkdir "$OUTPUT"

  apptainer exec --cleanenv \
    --bind "$INPUT:/input:ro,$OUTPUT:/output" \
    "$IMAGE" env RAYON_NUM_THREADS="$THREADS" \
    rustroo -i /input/bakta -o /output/whole_isolate \
      --threads "$THREADS" --audit-members
)
```

Input and output binds are explicit; the input is read-only. Apptainer runs
with the calling user's identity, so output files belong to that user.
To analyse plasmid bins instead, change the input to `/input/plasmid_bakta`
and use a new output directory. This treats each bin as a sample; for
family-specific analysis, provide only the bins of that family.

## Nextflow integration

No Nextflow pipeline source was present to modify. Add a process/module to
your existing workflow, using its existing scheduler, project/account, queue
and resource conventions. The container is selected by Nextflow; **do not**
put another `apptainer exec` inside the process script.

Example process:

```nextflow
process RUSTROO {
    tag "$cohort"

    input:
    tuple val(cohort), path(gff_dir)

    output:
    tuple val(cohort), path('rustroo'), emit: results

    script:
    """
    export RAYON_NUM_THREADS=${task.cpus}
    rustroo -i "${gff_dir}" -o rustroo \
      --threads ${task.cpus} --audit-members
    """
}
```

Example settings to merge into the pipeline's existing configuration:

```groovy
apptainer.enabled = true
apptainer.autoMounts = true

process {
    withName: RUSTROO {
        container = '/YOUR/PIPELINE/CONTAINER_DIRECTORY/rustroo-<commit>-linux-amd64.sif'
        cpus = 8
        memory = '8 GB'
        time = '2h'
    }
}
```

These resource numbers are starting examples, not measured requirements for
your dataset. Keep the site's existing Slurm profile rather than inventing a
second executor configuration. Disable another container backend for this
profile if the existing configuration enables one.

Important staging and sample rules:

- One task receives **all sample GFFs in the intended cohort**, not one task
  per isolate. For a plasmid family, one task receives all its member bins.
- Rustroo's directory input is **not recursive**. Stage a flat directory of
  `.gff`/`.gff3` files, optionally gzip-compressed.
- Features-only GFFs need their matching sibling FASTA staged in the same
  directory, with the same stem and a supported FASTA extension. Embedded
  `##FASTA` GFFs need no additional sequence file.
- Preserve unique sample filenames. Duplicate stems can produce duplicate
  sample names; do not rely on directory names to distinguish them.
- Nextflow normally stages symlinks. Ensure the original targets are visible
  through supported automatic mounts or explicit site-configured binds.
- Set **both** `RAYON_NUM_THREADS` and `--threads` to `task.cpus`.
  The existing CLI's thread option alone does not cover every Rayon stage.
- Keep cohort outputs separate, and retain the exact SIF checksum with the
  pipeline run. The existing `validate` subcommand is also available in the
  same image, but is not automatically run.

## Acceptance checks before production

1. The build completes, including all existing Rust tests and the container
   smoke tests; revision and checksums match the selected source.
2. `--version` and both build/validation help commands execute on an allocated
   MACS **x86-64** compute node.
3. Run a small known cohort once natively and once in the container with the
   same source, inputs and thread count. Compare `genes.fa`, `presence.tsv`,
   `lengths.tsv`, `reference_report.tsv` and the two audit exports byte-for-byte.
4. Run one small Nextflow cohort, confirm the expected sample count and exported
   files, and check the scheduler's CPU/memory accounting before scaling up.

No output-equivalence, image-size, performance or MACS compatibility result is
claimed until those checks are executed.

## Monocle GFF benchmark on MACS

One Slurm job runs the flat `monocle_gff` directory: 1,514 GFF3 files, the
same cohort as the local `monocle_run`. Refind stays on and `--audit-members`
stays on, matching that run. The job prints `cpu_min` and `actual_min`.
Default allocation is 16 CPUs, 64 GB, and 4 hours on `genoa` for account
`brics.b36ad`. That memory and time are headroom, not a measured requirement.
Change `--cpus-per-task` at submission; the script sets both Rayon and
`--threads` from `SLURM_CPUS_PER_TASK`.

Build the image with the block above, then run this from local WSL. It uploads
the SIF, its checksum, the job script, and the GFFs. Filenames containing `#`
are passed through `rsync` as files, not as a shell list.

```bash
cd ~/newdev/rustaroo
REV=$(git rev-parse --verify HEAD)
SHORT=${REV:0:12}
SIF="$PWD/out/containers/rustroo-${SHORT}-linux-amd64.sif"
test -s "$SIF"
test -s "$SIF.sha256"

ssh -F ~/.ssh/config_clifton b36ad.macs3.isambard 'mkdir -p "$SCRATCHDIR/rustaroo-macs/sif" "$SCRATCHDIR/rustaroo-macs/monocle_gff" "$SCRATCHDIR/rustaroo-macs/logs" "$SCRATCHDIR/rustaroo-macs/tmp" "$SCRATCHDIR/rustaroo-macs/apptainer-cache"'
REMOTE=$(ssh -F ~/.ssh/config_clifton b36ad.macs3.isambard 'printf %s "$SCRATCHDIR"')

scp -F ~/.ssh/config_clifton \
  "$SIF" "$SIF.sha256" \
  "b36ad.macs3.isambard:$REMOTE/rustaroo-macs/sif/"
scp -F ~/.ssh/config_clifton \
  containers/macs_monocle.sbatch \
  "b36ad.macs3.isambard:$REMOTE/rustaroo-macs/macs_monocle.sbatch"

rsync -a --info=progress2 \
  -e "ssh -F $HOME/.ssh/config_clifton" \
  "$PWD/monocle_gff/" \
  "b36ad.macs3.isambard:$REMOTE/rustaroo-macs/monocle_gff/"
```

On MACS, after `ssh -F ~/.ssh/config_clifton b36ad.macs3.isambard`:

```bash
cd "$SCRATCHDIR/rustaroo-macs"
( cd sif && sha256sum --check rustroo-*-linux-amd64.sif.sha256 )
find monocle_gff -maxdepth 1 -type f -name '*.gff3' | wc -l
sbatch --job-name=rustroo_monocle macs_monocle.sbatch
```

Expect `1514` from the file count. The Slurm log is `logs/rustroo_monocle_<jobid>.out`.
Results are `runs/monocle_<jobid>/`, including `rustroo.log`. After the job leaves
the queue, `sacct -j JOBID --format=JobID,State,ExitCode,Elapsed,TotalCPU,MaxRSS,ReqMem --units=G`
is the scheduler's elapsed time and CPU total. `cpu_min` and `actual_min` in the
log are the process measurement.
