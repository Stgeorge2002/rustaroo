# rustroo

Gene-level bacterial pangenome from Bakta or Prokka GFF3.

Given annotated genomes, rustroo builds ortholog families and writes:

1. **One intact nucleotide reference** per gene (for SNP / LoF mapping)
2. **Which samples carry that gene** (presence matrix)

The reference is the longest intact coding sequence that is still that gene — not merely the most common allele — so losses map as variants instead of disappearing into the reference.

## Requirements

- [Rust](https://rustup.rs/) (stable)
- Input: GFF3 with contig sequence in a `##FASTA` block, or a sibling FASTA

## Install

```bash
cargo install --git https://github.com/Stgeorge2002/rustaroo --locked
```

Or clone and build:

```bash
git clone https://github.com/Stgeorge2002/rustaroo.git
cd rustaroo
cargo build --release --locked
```

Binary: `~/.cargo/bin/rustroo` or `./target/release/rustroo`.

## Usage

```bash
rustroo -i path/to/gffs -o out/
```

Useful options:

| Flag | Meaning |
|------|---------|
| `--threads N` | Worker threads (`0` = default pool) |
| `--audit-members` | Also write `members.tsv` and `member_alleles.fa` |
| `--quality strict\|high\|standard` | Filtering strictness (default `strict`) |
| `--no-refind` | Skip searching intergenic gaps for missed genes |
| `--strict-input` | Fail on GFF/FASTA repairs the reader would otherwise apply |
| `--length-ratio R` | Minimum shorter/longer protein ratio in both homology passes (default `0.95`) |
| `--min-cds-bp N` | Minimum annotated and refound CDS length (default `90`) |

```bash
rustroo -i gffs/ -o out/ --threads 16 --audit-members
```

### Repair safeguards

Refinding accepts only closed ORFs meeting `--min-cds-bp`, at least 80%
shorter/longer protein length ratio, and at least 80% protein identity.
Candidates are ranked by identity weighted by length ratio. Tiny perfect
peptide matches cannot substitute for a whole gene. Refound copies do not
establish a finished gene end; refound-only alleles shorter than 95% of the
longest intact annotated allele are excluded from reference selection, but
their accepted memberships remain in the presence matrix.

Same-locus reconciliation uses replicated strand-normalized neighborhoods,
up to three frequent intact annotated protein representatives, and disjoint
sample sets. Comparable-length alleles require at least 90% identity; shorter
fragments require at least 30 amino acids and at least 98% terminal identity.
Competing hosts and distinct co-occurring loci are not arbitrarily merged.
Neighborhood evidence tolerates overlaps of up to 30 bp.

Refinding uses a sparse context graph keyed by the two flanking families.
For each flank pair, it considers up to eight middle-family candidates
ranked by support in distinct samples. A one-sample family can be considered,
but only a qualifying ORF sequence match in a target sample's matching gap
adds presence. Equal-scoring assignments of one ORF to different families
are rejected. This bounded lookup avoids comparing every family with every
gap; it can recover missed annotations for rare families without making
context alone count as gene presence.

Split/merge repair runs to stability with a 12-round safety limit, both before
refinding and, when genes are refound, afterward. Limit warnings and
`reconcile` proposal/rejection records are written to `rustroo.log`.
These safeguards intentionally change inference; use fresh output directories
and revalidate biological accuracy rather than treating Panaroo as ground truth.

### Validate a completed audit run

```bash
rustroo validate -i out/ -o out/validation/ --threads 16
```

Requires `genes.fa`, `presence.tsv`, `members.tsv`, and `member_alleles.fa`. Use a **fresh** output directory.

### Near-core split experiment

The standalone `near-core-splits` subcommand ranks possible splits in a completed
`--audit-members` run. It reports sample-presence unions, exact protein identity
and coverage of both references, replicated strand-normalized neighborhoods,
and co-occurrence at separate physical loci. It does not change clustering.

Reporting is the default. The opt-in `--test-merges` flag writes a separate,
presence-only simulation for at most `--max-merges` (default 10) strongest,
unambiguous pairs. Sample overlap, multiple loci, uncertain locus bundles,
insufficient sequence or neighborhood evidence, and competing partners veto a
merge. Original memberships and reference sequences are never rewritten.
An optional `--comparison-core-count` measures the remaining core-count gap;
the comparator must use the identical samples and core threshold.

See the [near-core split methods](benchmarks/README.md#near-core-split-method)
for thresholds, ranking, outputs, and interpretation limits.

## Outputs

| File | Content |
|------|---------|
| `genes.fa` | One reference CDS per gene |
| `presence.tsv` | Gene × sample presence (0/1) |
| `lengths.tsv` | Allele length per present sample |
| `reference_report.tsv` | Why each reference was chosen |
| `rustroo.log` | Full run record |
| `members.tsv` / `member_alleles.fa` | Optional membership audit (`--audit-members`) |

## Develop

```bash
cargo test --lib --bin rustroo
cargo build --release --locked
```

Container packaging: [containers/README.md](containers/README.md).  
Benchmark protocol vs Panaroo: [benchmarks/README.md](benchmarks/README.md).

## License

MIT — see [LICENSE](LICENSE).
