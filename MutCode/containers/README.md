# MutCode Apptainer image

Build once on **Linux x86-64**, then copy the resulting SIF and checksum to
MACS or another Linux x86-64 machine. Analyses do not compile Rust or require
Snippy, Cargo, Docker or internet access inside the image.

## Build

From this repository, with Git and Apptainer installed:

```bash
bash MutCode/containers/build.sh
```

The helper captures Git-trackable product files, including uncommitted files,
in an isolated source snapshot. It does not commit, push or publish anything.
It computes a SHA-256 source-manifest fingerprint rather than falsely
claiming that the current Git commit contains uncommitted MutCode changes.
The image builder validates every captured file against that manifest.
The snapshot is retained for inspection/reproduction.

The build uses rootless `--fakeroot`, which requires a working Apptainer
fakeroot configuration and subordinate UID/GID mappings. If your cluster
does not permit image builds, build on a supported workstation and transfer
the SIF instead. No fakeroot is needed to **run** this image.

Outputs go to `MutCode/out/containers/` by default:

- `mutcode-<source-fingerprint-prefix>-linux-amd64.sif`
- its `.sha256`, `.def` and `.build.log` sidecars.

Optionally pass an output directory as the first build-script argument.
Existing images are not overwritten. Download access to Docker Hub and
crates.io is required during a fresh build.

The recipe pins the Linux amd64 digests originally resolved for
`rust:1.90.0-bookworm` and `debian:bookworm-slim`. It uses the locked Cargo
dependencies, runs all Rust tests, and builds for baseline x86-64 rather
than this workstation's instruction set. The small runtime contains the
binary and provenance, not the Rust compiler or datasets.

## Transfer and verify

Copy the SIF and matching `.sha256` file together. In their destination
directory:

```bash
sha256sum --check mutcode-<fingerprint>-linux-amd64.sif.sha256
apptainer test mutcode-<fingerprint>-linux-amd64.sif
```

Replace the placeholder with the actual generated name. Record the full
image checksum in pipeline provenance. Source checksums, toolchain version
and binary checksum are also inside `/opt/mutcode-build`.

## Run

```bash
IMAGE=/absolute/path/mutcode-<fingerprint>-linux-amd64.sif
DATA=/absolute/path/analysis_inputs
OUT=/absolute/path/mutcode_results
mkdir -p "$OUT"

apptainer exec --cleanenv \
  --bind "$DATA:/inputs:ro,$OUT:/outputs" \
  "$IMAGE" mutcode encode \
  --reference /inputs/pan_genome_reference.fa \
  --snippy-dir /inputs/snippy \
  --output /outputs/mutation_codes.tsv.gz \
  --threads 8

apptainer exec --cleanenv \
  --bind "$DATA:/inputs:ro,$OUT:/outputs" \
  "$IMAGE" mutcode decode \
  --reference /inputs/pan_genome_reference.fa \
  --codes /outputs/mutation_codes.tsv.gz \
  --verify-only --threads 8
```

`apptainer run IMAGE.sif ...` is also supported. Nextflow uses
`apptainer exec` and finds `mutcode` on the image's PATH.

This is a portable local image, **not a published registry image**.
See [the reusable module for your existing Nextflow pipeline](../nextflow/README.md).

## Image built and tested in this checkout

```text
File:   mutcode-ce511d013e03-linux-amd64.sif
Size:   29,020,160 bytes (27.7 MiB)
SHA256: c683b7be83491bca849336cca8e85000ee7079751811c70bd1dce1a28ae61de4
```

It is under `MutCode/out/containers/` from the repository root. The source
manifest fingerprint is
`ce511d013e03831eb4c6d3cbd75aa53f93cd9cbb199507eb63b505588e74bf15`.
All 13 automated tests passed inside the builder. Runtime checks, native/image
output equivalence and an actual Nextflow 24.10.5 run over the supplied
90 samples passed. This checks container execution locally, not MACS's
scheduler or site policies.
