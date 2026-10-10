# Panaroo's published testing and benchmarking

This document describes the testing reported in Tonkin-Hill et al. (2020),
**Producing polished prokaryotic pangenomes with the Panaroo pipeline**.
It is a reproduction guide, not an executable benchmark. No datasets, tools,
simulations, or tests were run to prepare it.

The paper and public scripts were inspected. Numerical simulation parameters
from Supplementary Table 1 have **not** been transcribed or verified here:
retrieve Additional file 1 through the article before configuring those runs.
Do not substitute invented rates and call the result an exact reproduction.

## 1. Sources

- [Published article and figures](https://pmc.ncbi.nlm.nih.gov/articles/PMC7376924/).
- [DOI and publisher page](https://doi.org/10.1186/s13059-020-02090-4).
- [Full article XML, including Methods](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7376924/fullTextXML).
- Supplementary Table 1 and supplementary methods: **Additional file 1**,
  `13059_2020_2090_MOESM1_ESM.pdf`, linked from the article.
- [Panaroo repository](https://github.com/gtonkinhill/panaroo).
- [Public scripts inspected at commit
  e96b4978327536256b5c663f741138fc7aabed64](https://github.com/gtonkinhill/panaroo/tree/e96b4978327536256b5c663f741138fc7aabed64/scripts).

The inspected repository commit is a public source snapshot, **not verified as
the commit used in the 2020 paper**. Current regression tests are also separate
from the paper's experiments.

## 2. Overview of the evidence

| Experiment | Data | Main measurements | Evidence type |
|---|---|---|---|
| Simulated pangenomes, Fig. 3 | Seven scenarios, three replicates of 100 genomes each | Errors relative to a known simulated presence/absence matrix | Known-truth accuracy within the simulation model |
| Clonal outbreak, Fig. 2 | 413 London tuberculosis outbreak assemblies | Core/accessory counts and annotation/assembly discrepancies | Biological negative control, not a complete truth-labelled matrix |
| Diverse collection, Fig. 4 | 328 Klebsiella genomes; nine QC outliers removed | Core/accessory counts and conflicting annotations within clusters | Internal-consistency evidence, not ground truth |
| Computational performance, Fig. 7 | 10, 100, and 1,000 gonococcal isolates | CPU time and memory with five CPUs per tool | Performance/scaling |
| Association analysis, Fig. 5 | 1,054 Euro-GASP gonococcal isolates and antibiotic phenotypes | Gene/triplet associations with resistance | Downstream biological utility |
| Evolutionary analysis, Fig. 6 | 51 major pneumococcal sequence clusters with dated trees | Gene gain/loss and effective pangenome-size estimates | Downstream model/biological consistency |

The simulations provide the most direct published accuracy test. The other
experiments support usefulness and plausibility, but do not establish that
every individual inferred cluster contains only the correct orthologs.

## 3. Known-truth simulation benchmark

### 3.1 Reference and model

The starting reference was:

- Organism: *Escherichia coli*.
- Assembly: **ASM584v2**.
- Accession stated in the paper: **NC000913**.

Gene presence/absence variation was simulated using the **Infinitely Many
Genes (IMG) model**. A population phylogeny was generated with the **Kingman
coalescent** using DendroPy.

The authors varied:

1. Accessory gene gain rates.
2. Accessory gene loss rates.
3. Within-gene codon substitution rates.

Lower gain/loss settings produced a larger core and smaller accessory genome;
higher settings produced more accessory variation. Mutation settings tested
the effects of gene sequence divergence.

### 3.2 Scenarios and replication

The Results describe:

- **Five simpler simulation scenarios** varying gain/loss and mutation.
- **One increased-fragmentation scenario**.
- **One contamination scenario**.
- **Three dataset replicates per scenario**.
- **100 sampled genomes per dataset replicate**.

This is 21 dataset replicates, each containing 100 genomes. These are biological
simulation replicates, not three repeated timings of the same dataset.

The exact five parameter combinations, fragmentation settings, and any
additional scenario-specific choices must come from **Supplementary Table 1
and the supplementary methods**. This document does not infer them from the
figures or from current script defaults.

### 3.3 Assembly and annotation workflow

The Methods describe the following simulation-to-annotation pipeline:

1. Generate population gene-content and codon variation from the reference.
2. Preserve the simulated true gene presence/absence matrix.
3. For the sequencing/assembly experiments, simulate short reads using
   **Mason or ART**.
4. Assemble the reads using **SPAdes**.
5. Annotate the assemblies using **Prokka**.
6. Supply Prokka with a custom protein database containing correctly assigned
   proteins from the simulation **before assembly**.
7. Run each pangenome tool on the resulting annotations.
8. Compare inferred clusters and presence/absence with the simulated truth.

For each scenario, confirm the exact input branch and simulator from the
supplement rather than assuming every simpler scenario uses an identical
read-generation path.

The custom protein database is important: generic Prokka annotation without
that database is not the same annotation experiment.

### 3.4 Fragmentation scenario

The simulated complete genomes were broken into fragments **before read
simulation**, followed by assembly and annotation.

The paper reports a resulting **mean assembly N50 of 23,112 bases**.

That N50 is an observed outcome, not a sufficient specification of the
fragmentation procedure. Retrieve the cut-generation settings from the
supplement or archived workflow; do not select arbitrary cuts merely to hit
the reported N50.

### 3.5 Contamination scenario

The contaminating organism was:

- *Staphylococcus epidermidis*.
- Reference assembly **ASM764v1**.

The Methods describe randomly adding approximately **10-kb reference segments**
to simulated genomes before read simulation. The number of added segments was
sampled from a **Poisson distribution with mean 1**.

These foreign segments are not part of the intended *E. coli* pangenome truth.
Report their retention/removal separately from recovery of genuine rare genes.
Do not count all rare sequences as contamination simply because they are rare.

### 3.6 Software versions reported in the Methods

| Component | Reported version |
|---|---|
| DendroPy | 4.4.0 |
| Mason | 2.0.9 |
| ART | 2.5.8 |
| SPAdes | 3.13.0 |
| Prokka | 1.13.3 |
| Panaroo | 1.0.0 |
| PanX | 1.5.1 |
| Roary | 1.007002, as printed in the paper |
| PIRATE | 1.0 |
| COGsoft | v201204 |

The Results also compare PPanGGoLiN and explain that `--defrag` was enabled in
subsequent analyses after its default behavior inflated cluster counts.
Its exact version/settings should be recovered from the original materials;
they are not established by the Methods version list above.

The paper evaluates Panaroo in **strict and sensitive modes**. Sensitive mode
does not remove candidate contaminating gene clusters, so contamination
retention must not be interpreted identically across modes.

### 3.7 Scoring in the paper

The Methods state that inferred orthologous clusters containing an error were
counted by comparison with the simulated matrix.

Figure 3 and its surrounding text discuss:

- Erroneous clusters / accessory-genome inflation.
- Missing genes / false-negative calls.
- Incorrectly clustered or annotated genes.

The article alone does **not** provide enough detail to reconstruct every
cluster-matching and error-counting convention. Recover those rules from the
supplement and original analysis code before claiming to reproduce Figure 3.

In particular:

- Do not compare matrices by row order or family name.
- Do not use Panaroo's own groups as the truth for scoring Rustroo.
- Reference similarity alone cannot detect an incorrect merge or split.
- A true gene never annotated in the shared GFF input should be distinguishable
  from a gene incorrectly removed or grouped by a pangenome tool.
- Report refound, ambiguous, fused, and unmapped annotations explicitly.

As additional metrics for a new paired benchmark, report false merges, false
splits, presence/absence errors, contamination retention, and copy-resolution
errors. Label these as **additional metrics**, not necessarily the exact
original paper's definitions.

## 4. London tuberculosis outbreak: biological negative control

The authors analysed **413 Prokka-annotated assemblies** from a highly clonal,
isoniazid-resistant London *Mycobacterium tuberculosis* outbreak.

The maximum pairwise SNP distance was **nine**. Very little genuine pangenome
variation was therefore expected; Figure 2 describes an expected core of
approximately 4,000 genes.

They compared inferred core and accessory genome sizes across tools. Panaroo
reported the largest core and smallest accessory genome, consistent with the
clonal outbreak interpretation.

The paper also examined causes of disagreement:

- Assembly-driven gene fragmentation accounted for a large fraction.
- Inconsistent annotation could make nearly identical underlying sequence
  appear present in only a minority of samples.
- Some residual apparent absences arose when genes could not be refound.

**CheckM** was run as a separate assembly-QC assessment. The authors reported
that filtering the problematic assemblies would discard approximately 12% of
the collection and argued that correction allowed more assemblies to be kept.

This is a valuable negative control, but "clonal" does not mean every sample
must have precisely identical gene content. Do not force zero accessory genes
or optimize Rustroo until it reaches Panaroo's count.

Reproduction needs the original outbreak accession/sample list, assembly and
annotation settings, and the same QC inclusion decisions.

## 5. Diverse Klebsiella collection: internal consistency

The starting collection comprised **328 globally sourced
Klebsiella pneumoniae genomes** from human and animal hosts.

The Panaroo QC procedure identified **nine outliers** based on contig or gene
counts. They were removed **before running each algorithm**, leaving 319
genomes for that comparison.

Figure 4 compares:

1. Total pangenome size.
2. Core and accessory genome sizes.
3. Conflicting annotations within inferred gene clusters.

Core genes used a **99% sample-presence threshold**.

The authors treated clusters containing different informative annotations as
possible evidence of overmerging. Conflicts involving fragments or
hypothetical annotations were excluded in the paper's description.

This is an internal-consistency metric:

- Different annotations can describe homologous or multifunctional genes.
- Matching annotations do not prove orthology.
- Removing more genes or splitting every ambiguous family can artificially
  reduce conflicts.

The public annotation-conflict script differs in potentially relevant details:
its inspected version skips some multi-ID cells and refound IDs, normalizes
annotation strings, and has a commented-out length-based fragment filter.
Therefore, do not assume the current script implements every exclusion
described in the 2020 paper. Retrieve the historical scoring implementation.

For Rustroo, reuse the same source annotation descriptions and exclusions
rather than comparing only the single chosen reference description per family.

## 6. Computational scaling benchmark

Figure 7 used subsets of the **Euro-GASP Neisseria gonorrhoeae collection**:

- 10 isolates.
- 100 isolates.
- 1,000 isolates.
- **Five CPUs per tool**.

The reported measurements were **CPU time and memory**. Panaroo was compared
with competing pangenome methods; COGsoft did not complete the largest dataset
within a week.

The passages inspected do not establish the exact hardware, memory-accounting
method, subset seeds, or number of timing repetitions. Recover these where
available; otherwise identify them explicitly as reconstruction choices.

For a Rustroo/Panaroo comparison:

- Give both tools the same GFF files and the same CPU allocation.
- Keep annotation and simulation time separate from pangenome-tool time.
- Record wall, user, and system times separately.
- State whether container startup is included.
- State how memory is measured. GNU time's maximum-process RSS is **not**
  aggregate peak memory across all concurrent subprocesses.
- Use multiple timings and a saved order as an extension if desired, but do
  not describe invented repetitions as the original publication protocol.

The existing 44-genome, 26-thread comparison is useful, but it is **not**
the paper's 10/100/1,000-genome, five-CPU workload.

## 7. Downstream biological demonstrations

These were part of the paper, but are not direct known-truth ortholog tests.

### 7.1 Gene and structural-variant association studies

The authors used **1,054 Euro-GASP gonococcal isolates** from 20 European
countries, collected in 2013, with susceptibility measurements for seven
antibiotics.

Panaroo gene presence/absence and gene-triplet structural patterns were
analysed using **pyseer**. The paper reports 67 associated genes and 138
associated triplets at an adjusted significance threshold of 0.05 or lower,
with biological interpretation of resistance-related loci and plasmids.

Reproduction requires phenotypes, population-structure covariates, association
settings, multiple-testing conventions, and structural-pattern definitions.
Matching the number of significant associations is not proof of grouping
accuracy or of a causal resistance mechanism.

Rustroo's current gene matrix can support gene-content association work, but
it does not currently reproduce Panaroo's graph-based structural-triplet
workflow. Treat that as a separate implementation task.

### 7.2 Pneumococcal pangenome dynamics

The authors ran Panaroo separately on **51 major Global Pneumococcal Sequence
Clusters** for which reliable dated phylogenies could be constructed.

They used a **Finitely Many Genes (FMG) model** to infer gene gain/loss rates
and effective pangenome sizes, then examined relationships with recombination,
serotype diversity, and invasive-disease associations.

This requires the original cluster/sample lists, dated trees, model settings,
and metadata. It tests downstream consistency, not known family membership.

## 8. Relevant public regression tests

The following tests were inspected in the public repository snapshot. They are
software regression fixtures, **not substitutes for the publication's accuracy
benchmark**, and their data dependencies must be obtained separately.

### Identical-input test

[tests/test_clean.py](https://github.com/gtonkinhill/panaroo/blob/e96b4978327536256b5c663f741138fc7aabed64/tests/test_clean.py):

- Runs strict-mode Panaroo with `paralog.gff` supplied three times.
- Expects 5,782 presence-matrix rows.
- Expects every inferred family present in all three samples.

Rustroo deduplicates repeated input paths. An equivalent Rustroo fixture needs
three distinctly named file copies with identical content and recorded hashes,
not the exact same path repeated three times. The expected family count should
not automatically be treated as biological truth.

### Messy-annotation test

[tests/test_messy.py](https://github.com/gtonkinhill/panaroo/blob/e96b4978327536256b5c663f741138fc7aabed64/tests/test_messy.py):

- Runs strict-mode Panaroo on `aa1.gff` through `aa4.gff`.
- Expects 5,118 family rows.
- Expects 20,431 positive presence entries in total.

These assertions test reproducibility of Panaroo's expected output on those
fixtures. For Rustroo, inspect actual membership differences rather than simply
requiring those counts.

The current test configuration obtains fixture locations via `--datafolder`.
Without supplied fixture paths, tests parametrized by that option may not
execute on any data. Check test collection and executed-test counts explicitly.

Other public tests cover alignment, cleaning, integration, merging, MSA, QC,
plotting, and Spydrpick. They test Panaroo-specific APIs and outputs; they cannot
all be run unchanged against Rustroo.

## 9. Public scripts and reproduction hazards

The inspected snapshot contains:

| Script | Intended role |
|---|---|
| `pseudo_full_pangenome.py` | Coalescent/IMG gene-content and codon simulation; exports genome FASTAs, a tree, a truth matrix, and a custom Prokka protein database |
| `simulate_reads_assem_prok.sh` | ART simulation, SPAdes assembly, and Prokka annotation |
| `split_genome.py` | Random fragmentation of a single-contig genome |
| `pseudo_contam.py` | Addition of contaminating reference segments |
| `pseudo_diversity.py` | Artificial gene-sequence variation |
| `count_annotation_conflicts.py` | Annotation-conflict scoring with tool-specific matrix parsing |

Do not execute these blindly:

- `split_genome.py` uses Python-2 syntax and a historical dependency.
- The shell assembly script assumes directory names and a globbed Prokka
  database. It also deletes intermediate files. Use controlled, dedicated
  experiment directories and preserve provenance.
- The shell script's active ART command uses HS25 paired-end 150-base reads,
  20-fold coverage, fragment mean 200, and standard deviation 10.
  Its Mason command is commented out. These are **snapshot settings**, not
  verified settings for every paper scenario.
- The scripts use multiple random-number sources without a single explicit
  experiment seed interface. Reproduction needs recorded, coordinated seeds.
- Exact historical dependency versions and APIs may be necessary.
- The truth presence matrix does not by itself provide a complete
  sample-scoped mapping from every newly annotated CDS to its true gene.
  Preserve or reconstruct that mapping and flag ambiguities before scoring.

Archive the exact upstream commit, patches, tool versions, input hashes,
random seeds, command lines, and truth outputs. Modernizing a script can be
reasonable, but disclose changes instead of calling the result bit-for-bit
historical reproduction.

## 10. Practical reproduction checklist

### Recover the original specifications

- [ ] Obtain Additional file 1 and transcribe Supplementary Table 1.
- [ ] Recover exact reference files and assembly/accession versions.
- [ ] Recover the paper-era scripts and analysis/scoring conventions.
- [ ] Recover datasets, sample lists, QC exclusions, and tool command lines.
- [ ] Confirm which scenarios used which read simulator and input branch.

### Construct a fair paired experiment

- [ ] Generate each simulated replicate once.
- [ ] Annotate once and give both tools the exact same GFF inputs.
- [ ] Freeze independent gene-content truth and annotation-to-truth mappings.
- [ ] Score Panaroo and Rustroo independently against that truth.
- [ ] Evaluate Panaroo strict and sensitive modes separately.
- [ ] Record invalid-CDS, fusion, refound, and unmapped-call treatment.
- [ ] Keep biological simulation replication separate from timing replication.
- [ ] Preserve failed or interrupted runs; do not replace errors with zeros.

### Interpret the results

- [ ] Distinguish historical reproduction from a modern paired extension.
- [ ] Distinguish cluster errors from presence errors and contamination.
- [ ] Report uncertainty across independently simulated datasets.
- [ ] Do not optimize Rustroo thresholds merely to match Panaroo counts.
- [ ] Add duplication/copy-history tests separately: the published IMG
  benchmark is not a comprehensive test of all paralog evolutionary histories.

For a historical reproduction, use the paper-era Panaroo version. For a useful
present-day comparison, run the modern version as an additional arm. The local
`staphb/panaroo:1.5.1` image is **not** the Panaroo 1.0.0 version used in the
paper's Methods.
