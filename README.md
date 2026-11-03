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

```bash
rustroo -i gffs/ -o out/ --threads 16 --audit-members
```

### Validate a completed audit run

```bash
rustroo validate -i out/ -o out/validation/ --threads 16
```

Requires `genes.fa`, `presence.tsv`, `members.tsv`, and `member_alleles.fa`. Use a **fresh** output directory.

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
