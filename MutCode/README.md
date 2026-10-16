# MutCode

Standalone Rust encoder/decoder for **lossless gene-reference consensus codes**
from final Snippy VCFs. Rustroo's analysis and outputs are unchanged.

## Exact reconstruction target

For each VCF sample and each FASTA gene, MutCode:

1. Starts with the supplied reference gene sequence.
2. Applies the literal allele selected by the sample's genotype at every final
   VCF record.
3. Keeps the reference sequence at positions not changed by those records.

The result is the **called consensus**, not necessarily the complete original
assembly or the unknown true biological sequence. There are no coverage
decisions, confidence scores, consequence annotations or read processing.
An absent gene cannot be inferred from a variant-only VCF.

Panaroo's `pan_genome_reference.fa` and Rustroo's `genes.fa` are both accepted
FASTA formats. They are **not interchangeable references for an existing VCF**:
use the exact IDs and nucleotide sequences against which Snippy was run.
Matching gene names alone is insufficient.

Every VCF REF interval is checked against the FASTA; declared contig IDs and
lengths are checked too. These cannot establish the identity of unreported
reference bases, because most VCFs do not provide whole-reference checksums.
Decoding does verify every reference sequence against the SHA-256 recorded
when encoding.

## Build

Run in WSL/Linux from the repository root:

```bash
cargo test --manifest-path MutCode/Cargo.toml --locked --all-targets
cargo build --manifest-path MutCode/Cargo.toml --release --locked
```

The executable is `MutCode/target/release/mutcode`. This is a separate Cargo
package with its own lockfile, not a new Rustroo subcommand or dependency.
Use Rust 1.85 or newer.

For a compiled Linux x86-64 Apptainer image, use
[the dedicated build helper and container guide](containers/README.md).
A [reusable Nextflow process module](nextflow/README.md) is included for
integration into your own pipeline. MutCode is a standalone CLI/container,
not a Nextflow pipeline; no workflow or pipeline configuration is supplied.

## Encode your PRJNA1078342 results

```bash
(
  set -euo pipefail
  cd ~/rustaroo

  BASE="$PWD/macs_PRJNA1078342_full_20261012_132112_selected_results"
  RUN="$PWD/MutCode/out/PRJNA1078342_$(date +%Y%m%d_%H%M%S)"
  mkdir -p "$RUN"

  ./MutCode/target/release/mutcode encode \
    --reference "$BASE/panaroo/chromosome/panaroo_output/pan_genome_reference.fa" \
    --snippy-dir "$BASE/snippy" \
    --output "$RUN/mutation_codes.tsv.gz" \
    --threads 8

  ./MutCode/target/release/mutcode decode \
    --reference "$BASE/panaroo/chromosome/panaroo_output/pan_genome_reference.fa" \
    --codes "$RUN/mutation_codes.tsv.gz" \
    --verify-only --threads 8

  printf '\nMutation codes: %s\n' "$RUN/mutation_codes.tsv.gz"
)
```

`--snippy-dir` recursively finds final `.vcf`/`.vcf.gz` files, excluding
`.raw.vcf`, `.filt.vcf` and `.subs.vcf`. If both plain and compressed final VCFs
exist, the plain one is used once. Explicit `--vcf FILE [FILE ...]` inputs are
also supported instead of directory discovery.

Sample IDs are taken from the **VCF header**, not guessed from filenames.
For these inputs, IDs include the `_snippy_results` suffix.
Samples are sorted by ID; genes follow FASTA record order.
All reference genes get a row, including genes with no variant records.

The supplied dataset contains 90 final VCFs and 327,986 literal variant records:
260,020 SNPs, 3,895 MNPs, 1,157 insertions, 3,204 deletions and 59,710 complex
replacements. These counts describe the inspected inputs, not a general
Snippy output guarantee.

## The one main output

`mutation_codes.tsv` or `mutation_codes.tsv.gz`:

```text
sample  gene  reference_length  reference_sha256  target_length  target_sha256  codes
```

Columns are tab-separated. The first line is `#mutcode` followed by a JSON
metadata object containing the format version, reconstruction policy,
reference gene lengths/hashes, sample list, source VCF fingerprints and
haplotype selection. SHA-256 sequence hashes cover sequence bytes only, not
FASTA headers or line wrapping; nucleotide letter case is preserved.

`codes` is a compact JSON array:

```json
[[12,1,"T"],[20,0,"AC"],[30,3,""]]
```

Each edit is:

```text
[start, delete_length, inserted_sequence]
```

- `start` is a **0-based offset in the original reference**, never the
  progressively edited sequence.
- Delete the half-open interval `[start, start + delete_length)`.
- Insert the literal string at that point.
- Copy all unchanged reference intervals between edits.
- `[]` reconstructs the reference exactly under the called-consensus policy.

In this example, replace base 12 with T, insert AC before base 20 and delete
bases 30-32. Insertions at reference length are terminal insertions.

VCF's 1-based coordinates and anchoring bases are converted automatically.
Unchanged prefix/suffix bases are removed from each replacement; nearby
different records are not interpreted as independent evolutionary events.
No alignment, protein translation or lossy SNP-only conversion is used.
Repeat indels retain a deterministic representation derived from the VCF;
MutCode does not claim a biologically unique indel placement.

## Reconstruct sequences

```bash
./MutCode/target/release/mutcode decode \
  --reference EXACT_REFERENCE.fa \
  --codes mutation_codes.tsv.gz \
  --output reconstructed --threads 8
```

The output directory must not exist. It contains one full-variant FASTA per
sample and `samples.tsv` mapping original sample IDs to safe filenames
(`s000001.fa`, etc.). The file IDs avoid path traversal and filename collisions.
FASTA records include inserted bases and omit deleted bases; they are not
reference-length alignments.

Use `--sample EXACT_VCF_SAMPLE_ID` to write just one sample. Every code row is
still checked. `--verify-only` performs all reconstruction/hash checks without
writing FASTAs.

## Optional independent consensus check

If a full-variant consensus FASTA is available, verify it while encoding one
selected sample:

```bash
./MutCode/target/release/mutcode encode \
  --reference EXACT_REFERENCE.fa \
  --vcf FINAL_SAMPLE.vcf \
  --output sample_codes.tsv.gz \
  --check-consensus FULL_VARIANT_CONSENSUS.fa
```

The check compares sequence bytes by gene ID, not FASTA wrapping/order.
The consensus must contain exactly those gene IDs and use the same calling
policy. A substitutions-only consensus, masked consensus or reference-length
aligned FASTA is **not** the same reconstruction target.

The supplied Snippy result folders contain VCFs but no consensus FASTAs.
MutCode does not pretend to recover sequences or masks absent from those files.

## Alleles and genotypes

- SNP, MNP, insertion, deletion and complex replacement records are all
  handled by their literal REF/ALT strings. INFO `TYPE`, SnpEff `ANN`, depth
  and quality annotations do not change the encoded sequence.
- Multi-allelic ALT records select the genotype's actual allele index.
- Haploid genotypes and homozygous genotypes (`1`, `1/1`, `2/2`, etc.) work.
- Phased mixed genotypes need `--haplotype 1` or `--haplotype 2`.
  The file metadata records that choice.
- Unphased mixed or missing genotypes fail: they do not define one exact
  sequence without additional policy/evidence.
- `FILTER=PASS` and `FILTER=.` are accepted. The latter is used throughout
  these final Snippy files. Other FILTER values cause an error instead of
  silently dropping records.
- Selected symbolic alleles (`<DEL>`, `<INS>`, etc.), `*` spanning deletions
  and breakends fail explicitly. They are not the five literal Snippy
  mutation classes and do not supply a standalone literal target sequence.
- Conflicting/overlapping effective edits and duplicate calls fail explicitly.
  Unsorted, nonoverlapping VCF records are sorted per gene.

MutCode encodes every literal nucleotide change, including IUPAC ambiguity
symbols and non-triplet changes. It does not trim terminal stops, incomplete
codons or start sites. It does not invent unknown-base masks or convert
zero coverage into a deletion.

## Reliability and efficiency

- Each nonempty edit list is serialized, decoded and compared byte-for-byte
  with its constructed target before the table is published.
- Every row has reference/target lengths and SHA-256 hashes.
- Decode rejects a wrong reference, altered edits, missing/duplicate/reordered
  rows, truncation or unexpected trailing data.
- Output files must be fresh. Encoding uses a same-directory staging file
  and publishes it without overwriting an existing destination only after
  all samples succeed.
- Inputs are read-only. VCF fingerprints are checked after parsing;
  reference/verification-consensus files are rechecked before publication.
- VCF parsing is parallel across samples. Gene-code creation, reconstruction,
  round-trip checks and hashing are parallel across genes, including for a
  single sample. Decode and verify-only also parallelize reconstruction/hashing.
- Output remains ordered and deterministic: identical inputs produce
  byte-identical tables, including gzip, irrespective of thread count.
- Gzip and file writing are streamed through a single ordered writer.
  References are shared, unchanged alleles avoid sequence copies, VCF calls
  are bounded to one sample batch, and row work is bounded to 512 genes.
  Compression/input decompression are not parallelized.
- Default worker count is the available CPU count capped at eight;
  `--threads` controls an explicit fixed-size pool for either command.
  No nested pools or globally configured Rayon workers are used.
  This tool requires no external
  mapper, Snippy runtime, SnpEff installation or alignment program.

## Verified on the supplied data

The development validation used the supplied 90 final Snippy VCFs and the
7,005-record Panaroo reference:

- Encoded and decoded/verified **630,450 sample/gene rows**, applying all
  **327,986 variant records**.
- After adding per-gene threading, native encoding took **2.88 seconds**
  with one worker and **2.50 seconds** with eight. Maximum-process RSS was
  **18,808 KiB** and **26,384 KiB**, respectively.
- Native verify-only took **0.91 seconds** with one worker and **1.30 seconds**
  with eight. Parallelism does not guarantee a speedup: streamed decompression,
  writing and worker overhead can dominate. These are single-run observations
  on this WSL machine, not a cross-machine performance guarantee.
- One- and eight-worker compressed tables matched byte-for-byte, including
  the original pre-threading validation output.
- The compressed code table is approximately **34 MiB**.
- Independently compared all **7,005 genes** for `SRR31742274_snippy_results`
  with `bcftools consensus` 1.10.2 from the installed Snippy 4.6.0 container.
  Every sequence matched byte-for-byte.
- That independent sample contains 3,024 SNPs, 51 MNPs, 16 insertions,
  40 deletions and 808 complex replacements.
- The automated tests additionally cover multiallelic/phased selection,
  unchanged genes, ambiguous symbols, gzip, deterministic thread output,
  invalid/missing genotypes, overlaps, corruption and wrong references.
- All **13 automated tests** passed natively and inside the image builder;
  formatting and Clippy with warnings denied also passed.
- The dedicated **27.7 MiB Apptainer image** decoded/verified the full table.
  A development-only Nextflow/Apptainer integration test exercised the module
  on all 90 samples; its compressed output matched native encoding
  byte-for-byte. Single-sample
  staging and two different VCFs sharing the same filename were also tested.
  MACS scheduler/account configuration has not been tested locally.

Local validation output is at
`out/PRJNA1078342_validation/mutation_codes.tsv.gz` relative to this directory.
Outputs and build artifacts remain ignored by Git.
The threading/container/Nextflow validation outputs are under
`out/threaded_validation/`; see [the container guide](containers/README.md)
and [the Nextflow guide](nextflow/README.md).

The guarantee is lossless reconstruction of the stated **called consensus**.
It is not a claim about absent genes, distinct physical copies, confidence,
orthology, or the biological truth of the upstream variant calls.
