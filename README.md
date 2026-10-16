# rustroo

Gene-level bacterial pangenome from Bakta or Prokka GFF3.

## Goal

Build one ortholog group per gene across a set of genomes, then emit:

1. **A single nucleotide reference** for each gene — the sequence that downstream SNP and loss-of-function calling should map against.
2. **Which samples carry that gene** — presence of the gene family, not of one allele.

The reference is the longest intact coding sequence that is still that gene: in frame, a bacterial start, no internal stop. It is not the most common allele. A rare full-length copy beats a truncation that is present in every other sample, so a loss is called as a variant instead of being baked into the reference. Internal insertions and fusion tails that match another pangenome gene are rejected.

## Why

PanGLoFi and similar LoF pipelines need a stable gene reference and a sample × gene matrix. Panaroo provides that today. rustroo aims at the same job in Rust: cluster annotated CDS into families, repair fragments and fusions with synteny, optionally refind missed genes in intergenic gaps, then write clean outputs without a graph file.

## Independent Rust validation

`rustroo validate` examines a completed `--audit-members` run without rebuilding
the pangenome or changing its memberships, references, presence, or lengths.
Ordinary builds do not run any validation alignments or trees. Validation
itself costs additional CPU, memory, and I/O; no zero-cost or biological
correctness guarantee is implied.

```bash
cargo test --lib --bin rustroo && cargo build --release --locked

RUN="$PWD/full44_allele_rescue_20261010_001601"
CHECK="$RUN/validation_$(date +%Y%m%d_%H%M%S)"
./target/release/rustroo validate \
  -i "$RUN/rustroo" -o "$CHECK" --threads 26
```

The inputs are `genes.fa`, `presence.tsv`, `members.tsv`, and
`member_alleles.fa`. Sample/contig/locus identity is preserved. Duplicate rows,
missing sequences, conflicting locus assignments, and presence/membership
disagreements fail explicitly. Validation requires a **fresh** output directory.
A failed attempt is retained, but only successful completion writes
`summary.tsv` with `complete=true`.

### All-member sequence and structure evidence

Every distinct member DNA allele is examined and its result is linked to all
occurrences. Exact DNA aliases are pooled; identical translated proteins reuse
the same within-family alignment while retaining their distinct nucleotide
variation. Protein comparisons use exact global affine-gap traceback
(match 2, mismatch -1, gap open -5, extension -1), not the clustering fitter.

- `alleles.tsv`: full allele inventory, occurrence counts, paired canonical
  amino-acid identity, canonical coverage of both sequences, reading-frame and
  start/stop status, unknown residues, and structural gap totals.
- `structural_blocks.tsv`: locations of terminal reference gaps, internal
  reference gaps, and allele extensions/insertions. Amino-acid intervals are
  1-based inclusive; insertion boundaries are 0-based reference boundaries.
  An aligned gap is **not proof** of an evolutionary deletion or a lost domain.
  A shared initiator methionine can anchor residue 1 of both proteins and make
  an N-terminal truncation appear as a gap starting at reference residue 2.
  Interpret these near-start gaps alongside terminal sequence placements and
  genomic boundary evidence; do not automatically call them internal deletions.
- `families.tsv`: links gene IDs to filesystem-safe `fNNNNNN` alignment IDs.
- `alignments/*.faa`: all successfully aligned unique alleles, anchored on the
  family reference. Insertion blocks are padded for display; their columns are
  not assumed to be homologous across different alleles.
- `alignments/*.core.fna`: codon projections onto reference amino-acid
  positions, preserving synonymous variation and excluding inserted columns.
- `allele_aliases.tsv`: lossless mapping of input allele IDs to pooled DNA tips.

Named protein domains are explicitly **unassessed**: no Pfam/InterPro models or
domain annotations are supplied. These are sequence-structure checks, not
named-domain predictions or a general progressive multiple-sequence alignment.
Equally scoring traceback choices are flagged as ambiguous.

### Physical loci and copies

`loci.tsv` reports every exported physical locus, its coordinates and strand,
all member models, normalized observed flanks, origins, and coordinate/DNA
length disagreements. Bundles with inconsistent gaps/overlaps or touching
full-length models are flagged for review. Missing flanks are unknown, not
asserted contig boundaries: membership exports do not include contig lengths.

`copies.tsv` distinguishes multiple same-sample loci from overlapping/adjacent
models of one exported locus. It reports intact near-full-length model counts
and spatial separation, including copies on different contigs with identical
coordinates. This is copy evidence, not automatic proof of paralogy; source
contig sequences/reads and independent annotation review remain useful.

### Diagnostic gene phylogenies

The default `--tree-mode uncertain` requests trees for sequence/copy review
cases and candidate homologous family sets. Cross-family discovery is
restricted to shared oriented flanks, >=50% reference protein length ratio,
and >=10% shared distinct 4-mers of the smaller set. Exact pairwise assessment
then requires >=70% paired canonical identity and >=80% shorter coverage.
`related_families.tsv` reports tested pairs and outcomes. These filters are
**not exhaustive homology detection**, and connected candidate sets are not
new ortholog groups.

For each selected set, unique DNA haplotypes are aligned to its longest
reference. Trees use nucleotide sites that are unambiguous and covered in
**every** allele, excluding insertion columns and terminal stops. This avoids
inventing distances for nonoverlapping fragments. Out-of-frame members,
ambiguous alignments, insufficient common sites, no variation, saturation,
and resource limits produce explicit no-tree reasons.

- `trees.tsv`: every requested set, method outcomes, common/variable site
  counts, bootstrap availability, and negative-branch warnings.
- `trees/*.nwk`: **unrooted neighbor-joining** trees using JC69 nucleotide
  distances; negative branches are retained, not silently clamped.
- `trees/*.core.fna`: nucleotide alignment underlying a written tree.
- `tree_branches.tsv`: canonical internal bipartitions and site-bootstrap
  support. The denominator is all requested replicates; uninformative or
  saturated resamples do not create arbitrary supported splits.
- `tree_members.tsv`: maps each haplotype tip to every gene/sample/contig/locus,
  original annotation ID, and origin, including identical copies.

The default is 100 deterministic site bootstraps. This simple distance model
and bootstrap do not account for recombination, compositional bias, or
gene/species-tree discordance. Rooted duplication/speciation reconciliation is
explicitly **unassessed**; no independent species/strain tree is provided.
Identical copies cannot be resolved by a sequence tree. No tree changes any
assignment or certifies orthology.

### Resource limits and completion

`--tree-mode all` also requests ordinary single-copy family trees;
`--tree-mode none` disables tree construction but keeps sequence/copy reports.
`--bootstrap-replicates 0` disables support estimation. Defaults are:

- `--max-alignment-cells 5000000`: approximately 10 MB of traceback per active
  nonidentical protein pair, plus rolling rows. Larger pairs are explicitly
  unassessed; there is no heuristic substitute.
- `--max-tree-alleles 128`: larger candidate sets get an explicit no-tree result.
- `--min-tree-sites 90`: minimum complete-case nucleotide sites.

Setting either maximum to `0` removes that limit and can be expensive. Inspect
unassessed counts before treating an audit as covering all alleles.
`complete=true` means report generation completed, **not that all checks passed
or all biology was resolved**.

Allele scoring, candidate-pair scoring, and tree sets use the requested Rayon
pool; decisions/report ordering are deterministic. Outputs include settings,
method labels, and before/after FNV-1a input fingerprints to detect accidental
input changes. These fingerprints are noncryptographic, not a substitute for
archiving SHA-256 checksums and the exact executable/source version.
Semicolon-delimited identifier lists percent-encode `%`, `;`, and `:`.

## Reconciling shortened alleles

Strict homology length filtering is unchanged. During synteny repair, a shorter
representative longer than 80% of a candidate can now join that candidate when:

- Both groups have the same two known modal flanking families.
- No sample carries both groups.
- At least 98% of residues match over the entire shorter protein in an ungapped
  placement within two residues of either end of the longer protein.
- Exactly one longer candidate satisfies those rules.

Sample overlap is checked again against accumulated groups during application,
so a chain of individually compatible proposals cannot join co-occurring copies.
Equal-length proteins, interior domains, and internal indels are not handled by
this original terminal rule. Shorter fragments at or below 80% retain the
existing containment and neighborhood rules, except that multiple qualifying
containment hosts are now rejected instead of choosing the host present in
the most samples. Ambiguous terminal hosts also cannot fall through to a
containment choice. This is conservative allele reconciliation, not a target
gene count or an attempt to reproduce Panaroo's strict filtering.

Flank and adjacency indices restrict merge searches to the same context-eligible
families as the previous exhaustive scan. New terminal comparisons inspect at
most six placements in linear sequence time, without dynamic-programming
alignment. Large repeated-flank buckets can still produce many candidates;
real-data runtime and accuracy must be measured rather than assumed unchanged.
Merge logs report total group unions and the guarded shortened-allele subset
each round.

### Strand-aware and contig-edge allele rescue

An additional indexed rescue handles terminally shortened annotations whose
neighborhoods were reversed with their contigs or lost a flank at a contig
boundary. It does not change homology identity or length-ratio thresholds.

- Flanks are normalized into transcriptional upstream/downstream order using
  the central CDS strand. This applies to the new rescue, not recursive
  refinement of paralog context anchors.
- Every occurrence in the shorter family must be one whole annotated/refound
  CDS and at most one physical locus per sample. Bundled fragment pieces and
  unresolved multicopy families are not whole-allele rescue evidence.
- Every shorter protein must have at least 30 residues, cover at least 50% but
  less than 100% of the longer representative, and match >=98% over its
  **entire** length in an ungapped placement within two residues of a terminus.
  Unknown residues and internal stops are not accepted by this rescue.
- Both known oriented flanks must match. One known flank is sufficient only
  when the missing side is within 90 bases of the corresponding physical
  contig end and the locus is outermost on that side of the contig. The
  protein must also be shortened on that same side: loss of
  upstream sequence requires an N-terminal deletion, not a C-terminal one.
  A missing annotation or large intergenic gap inside a contig does not
  substitute for boundary evidence.
- The target must have an observed intact longest representative and at least
  two intact, >=95%-representative-length occurrences in distinct samples with
  the supporting complete flank pair. Every shorter occurrence's context and
  protein must qualify; modal representatives alone cannot hide a conflicting
  minority allele or neighborhood.
  Intactness uses the existing in-frame, valid-start, no-internal-stop rules:
  absence of a terminal stop alone does not invalidate a contig-ending gene.
- Exactly one sequence/context host may qualify. Singleton alternative hosts,
  damaged full-length hosts, bundled or unresolved multicopy hosts, and
  co-occurring hosts still count as ambiguity evidence rather than being
  silently ignored in favor of a frequent target.
  The accepted target itself must be single-locus per sample, and accumulated
  sample sets must be disjoint when the union is actually applied.

This targets exact contig-edge suffix matches such as the audited `abrB_2` and
`acrA_6` cases. It deliberately does not claim to repair the repeat-related
`aceF` discrepancy, internal deletions, two-ended domain fragments, or
low-identity alleles. Those require different structural evidence rather than
relaxing this rule. Ambiguous decisions and proposed strand-aware adoptions are
explicitly logged. A proposal is not an accepted union: accumulated sample
overlap can still block it during application.

Searches use observed flank-pair/side indices, reuse the physical-locus scan
already needed for merge voting, cache intactness by pooled DNA ID, and cache
terminal comparisons by pooled protein IDs and placement side. No new
dynamic-programming alignment, all-family protein search, or external program
is added. There is additional linear metadata work; large repeated-flank
buckets can still cost more and runtime must be measured.

## Stable fragment and paralog repair

Splitting, merge application, modal-flank voting and refinding now share a
physical-locus model. Consecutive annotations in one family on the same strand
and contig form one locus when they overlap at least 95% of the shorter
annotation, or abut within three bases and have fragment evidence. Fragment
evidence is either complementary proteins (global identity below 40%), or a
>=98%-identity full-coverage terminal match of their joined proteins to a longer
observed representative, covering at least 95% of that representative's length.
The latter protects repetitive fragments without treating two full tandem
copies as a broken gene. Different strands, contigs and intervening annotations
remain separate. External flanks and one vote per locus prevent fragment halves
from voting for artificial self-neighbor contexts.

Homology context anchors are captured after fusion repair and are independent
of the mutable ortholog-group IDs. Splitting a neighboring family does not
automatically become new evidence for splitting this family on the next pass.
Accepted fragment unions may coarsen these anchors across their descendants,
but never merge the descendants' ortholog assignments just because their
anchors agree. Merge candidate selection and refinding still use current
ortholog IDs, not these coarser anchors.

Supported physical bundles are retained independently of the current group's
representative. A repetitive fragment pair does not become two loci merely
because a subsequent paralog split moved the supporting full allele into
another group. Only previously supported bundles receive this protection;
unproven touching copies remain separate. Refound occurrences inherit their
target's context anchor but do not inherit another occurrence's physical locus.

Paralog splitting starts only with multiple physical loci in a sample.
Shared-flank context variants can be reconciled when their samples do not
overlap; accumulated sample constraints prevent transitive consolidation of
co-occurring loci. If candidate hosts themselves co-occur, placement is
ambiguous and is not assigned arbitrarily. Single-copy samples without a
resolvable host remain together in an unresolved family, not one family per
flank variant. Copies with indistinguishable contexts remain unresolved in
their family. Counts of unresolved placements are explicitly logged, so output
groups are not a claim of fully resolved orthology.

After flank reconciliation, exact whole-protein alleles can bridge contexts
whose two flanks have both changed. This second pass uses interned protein IDs,
not additional alignments. Eligible roots must have disjoint sample sets;
co-occurring candidate hosts make the assignment ambiguous, and accumulated
overlap is checked again during union. Individual pieces of a bundled fragment
are not whole-allele evidence. This does not override an earlier flank-based
assignment or resolve genuinely indistinguishable copies by iteration order.
Copy-specific neighbor labels are not recursively promoted into orthology
evidence: unresolved neighborhoods remain a documented limitation rather than
an excuse to force a smaller gene count.

All proposed fragment unions, including the older containment and abutting
rules, are checked for distinct co-occurring loci after earlier unions have
been applied. Co-occurring annotations may join only with the same physical
locus evidence; a low overall co-occurrence fraction no longer licenses merging
distant copies. Near-identical whole-protein paralogs (>=98% identity and >=95%
length ratio) are also excluded from the tail-only fusion-reference test.
Independent shorter genes matching a tail still reject the extended reference.

Synthetic reference joins now require a single supported physical locus on
one contig, on one strand, with non-overlapping pieces separated by at most
three bases in coding order. Coincident coordinates on different contigs,
unbundled tandem copies, and overlapping alternative annotations are not
concatenated into a reference. Per-sample covered lengths are unioned within
each contig and then summed, rather than conflating equal coordinates on
different contigs. The final locus-ID scan is shared by reference joining and
optional membership export.

The common single-copy split path is linear in its locus count. Reconciliation
is indexed by flanks within each multicopy family, not across all families.
Local protein comparisons are cached per protein-ID pair, joined-protein
checks allocate no sequence copy, and union-find uses path compression without
a fixed chain-depth limit. Large repeated-flank buckets can still have
quadratic candidate counts, as can repeated exact-allele buckets. Context
aliases and retained loci add linear scans after accepted merges and a small
amount of per-call metadata, not copies of protein or DNA sequences.
Shared-sample merge checks also add work; no unchanged-runtime guarantee is
implied.

Logs show family counts after fusion repair, additions from paralog splitting,
exact-allele context unions, counts after each round's merges, blocked unions,
and post-refind splits. With stable anchors and retained locus evidence,
splitting is idempotent until a real merge or a refound occurrence changes the
evidence. A pass without accepted merges therefore stops even if that pass
split families; it does not spend further rounds on relabeling-induced splits.
The four-round budget for merge-driven repair remains a safety limit, with an
explicit notice if accepted merges are still occurring at that limit.

## Usage

```bash
rustroo -i path/to/gffs -o out/
```

Inputs are GFF3 files (or a directory of them), with contig sequence in a `##FASTA` block or a sibling FASTA. Outputs:

| File | Content |
|------|---------|
| `genes.fa` | One reference CDS per gene (header is the gene id only) |
| `presence.tsv` | Gene × sample presence (0/1) |
| `lengths.tsv` | Allele length in each sample that has the gene |
| `reference_report.tsv` | Why each reference was chosen |
| `members.tsv` | Optional original CDS membership, coordinates, origin, locus and allele ID |
| `member_alleles.fa` | Optional deduplicated coding sequences referenced by memberships |

### Auditing excess families

Use `--audit-members` to preserve the information needed to compare actual
annotation memberships with Panaroo:

```bash
rustroo -i gff_test/chromosome -o audit_run/ --audit-members
```

The two additional files are opt-in and do not change clustering, references,
presence or repair thresholds. A fusion piece retains its original GFF CDS ID;
refound occurrences have no original ID and are marked explicitly. Multipart
CDS retain the decoded original ID. Coordinates describe the current repaired
annotation, and physical-locus IDs are local to this run, not stable IDs
across runs. Allele sequences are oriented coding DNA, not genomic intervals.
Missing GFF IDs are not replaced with guessed matches.

The [family audit workflow](benchmarks/README.md#annotation-level-family-audit)
joins sample-scoped annotation IDs, checks member-sequence agreement against
Panaroo's gene data, and reports sparse family relationships and reference
alignments. It retains ambiguous sequence hits and isolates fusion/refind
exceptions; neither presence overlap nor a similar reference is used to invent
a shared annotation. A lower group count or agreement with Panaroo alone does
not establish biological accuracy.

Membership export costs an extra linear locus scan, provenance storage, and
optional sequence/TSV I/O. The comparative audit and BLAST run are separate from
the normal pangenome pipeline. Use fresh output directories: a run without
`--audit-members` refuses existing membership files rather than leaving stale
audit evidence beside new outputs.

For library callers, `Params::audit_members` defaults to false,
`ParsedGene::annotation_id` supplies provenance, and
`Pangenome::membership_audit` holds the optional export data.

## Build

```bash
cargo build --release
```

For a dedicated Linux x86-64 Apptainer image and MACS/Nextflow integration,
see the [container packaging guide](containers/README.md). Compilation and tests
run only when building the image, never inside analysis tasks.

Regression tests, including moderate truncations, full-reference/presence
preservation, competing targets and co-occurring paralogs:

```bash
cargo test --lib
```

### Repository hygiene

Keep `Cargo.lock` in version control: Rustroo is an application, and its lockfile
pins dependencies for reproducible builds.

The `.gitignore` ignores everything by default and explicitly allows the
application manifests, Rust source and tests, benchmark scripts/documentation,
container recipes/documentation, and GitHub build/release workflows.
Local datasets, downloads, run outputs, SIF images and downloaded build tools
stay excluded regardless of their root-directory names.

Keep deliberately maintained small biological test fixtures under
`tests/fixtures/`, not among local datasets. New product assets or directories
outside the allowlist need an explicit rule before they can be committed.
`git status --short` and `git check-ignore -v --no-index PATH` help check this.

Python caches, local virtual environments, editor backups and `.env` settings
remain excluded even inside allowed directories. Sanitized root-level
`.env.example`, `.env.sample` and `.env.template` files remain trackable;
Cargo registry credentials do not. Ignore rules do not remove already tracked
files or prevent `git add -f`, so inspect staged changes before publishing.

The `.gitattributes` file normalizes text to LF across Windows, WSL and Linux.
This prevents line-ending-only changes from failing the container build's
committed-source check.

## Comparing with Panaroo

The [publication benchmark protocol](benchmarks/README.md) provides seeded,
nested chromosome subsets at 3, 5, 10, 15, 20, 25, 30, 35, and 40 genomes,
repeated randomized runs against `staphb/panaroo:1.5.1`, resource measurements,
and reciprocal reference-sequence comparisons. Preparing the experiment does
not execute the benchmark; running it is a separate explicit command.

## Standalone nucleotide reconstruction codes

[MutCode](MutCode/README.md) is a separate Rust encoder/decoder for final Snippy
VCFs and the exact gene-reference FASTA used for calling. It produces one
lossless edit-code row per sample/reference gene, with reconstruction hashes.
It does not change Rustroo, call variants, or infer gene absence or coverage.
