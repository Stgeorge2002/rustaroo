# MACS pneumococcal benchmark

Neither tool is a gold standard. A pair is shared only when it is a reciprocal best BLASTN hit at ≥70% identity and ≥90% coverage of both sequences (evalue 1e-10). Capture is shared pairs divided by that tool's references. Mean identity is over shared pairs only.

## Methods, step by step

Both start from the same annotated CDS and emit one nucleotide reference per gene plus a presence matrix. rustroo does that by clustering and local synteny. Panaroo does it by clustering and then cleaning a gene-neighborhood graph. Defaults below are the ones this run used.

1. **Filter the input calls.** rustroo drops a CDS shorter than 90 bp. Panaroo `--remove-invalid-genes` drops calls that fail its gene checks.
2. **Cluster by sequence.** rustroo clusters proteins at 98%, then at 70%, and keeps a pair only when the lengths are within 5%. Panaroo clusters at 98% (`--threshold`), then links families at 70% (`--family_threshold`).
3. **Use gene order.** rustroo looks at immediate neighbors: it splits paralogs that sit in different neighborhoods, and merges a fragment into the one gene it belongs to. Two intact genes are not merged, including through a scrap between them. Panaroo builds a graph whose edges are gene adjacency, then `--clean-mode strict` collapses split annotations and removes rare nodes with weak context.
4. **Recover missed genes.** Both search the assemblies for a copy the annotation omitted. rustroo searches intergenic gaps. Panaroo searches where the graph expects the gene.
5. **Choose the reference.** rustroo publishes the longest intact allele. A rare tail past a stop the rest of the collection shares is refused. Panaroo publishes one sequence per cleaned graph node. A node that absorbed more than one gene call is named `A~~~B`.
6. **Last pass in this run.** rustroo `strict` also drops a gene with no intact allele or an unconfirmed short length, and folds a short near-copy into the single longer gene that contains it. A full-length paralog below 90% identity stays. Panaroo's strict clean is the graph pass in step 3. rustroo labels a gene `core`, `accessory`, `unsupported`, or `refound`, and does not drop `unsupported`.

## What was run

- *Streptococcus pneumoniae* GFF assemblies (Monocle), nested cohorts from 2 to 1,500 genomes. 32 threads.
- rustroo `--quality strict`. Panaroo `--clean-mode strict --remove-invalid-genes`.
- Intact neighbouring genes are kept separate. Joined-fragment references: 0 at n=15 and n=1,500.
- UniProt check is a gene-name lookup for *S. pneumoniae* (taxon 1313). A length matches when reference bp equals amino acids × 3, or that plus the stop codon. Names UniProt does not return are excluded.

## Shared references

n=15. 2,514 shared pairs, mean identity 99.3%. Panaroo 2,739 genes, 91.8% captured. rustroo 2,716 genes, 92.6% captured. Wall time 0.7 min vs 0.8 min.

n=1,500. 3,390 shared pairs, mean identity 98.4%, 1,839 identical full length. Panaroo 3,969 genes, 85.4% captured. rustroo 4,688 genes, 72.3% captured. Wall time 8.8 min vs 22.9 min. CPU time 70 min vs 192 min.

rustroo publishes more genes, so a smaller share of its references are captured. Identity of the shared genes stays about 98.4%.

## Lengths at n=1,500

Of the 3,390 shared pairs, 3,143 are the same length. Mean absolute difference 1.6 bp. Panaroo is 1.5 bp longer on average. Where the 247 lengths differ, the mean gap is 22 bp.

Across all 3,573 reciprocal best hits, including pairs too different in length to count as shared, the mean absolute difference is 14 bp. Panaroo is longer by 12 bp on average.

## UniProt lengths at n=1,500

rustroo: 879 names found, 674 length matches (76.7%). Panaroo: 854 names found, 579 length matches (67.8%).

On 817 names both tools share and UniProt knows: both match 549, rustroo only 88, Panaroo only 17, neither 163.

Large rustroo misses are the plain name on a short copy. The UniProt length is also published: `polC_2` 4,392 bp, `nrdD_2` 2,208 bp, `yesN_2` 1,287 bp, `comA_2` 2,154 bp. Panaroo's `polC` is that 4,392 bp sequence; rustroo's `polC` is what Panaroo calls `polC~~~dnaQ`. A few-codon gap (`murB` 906 vs 951, `dcm` 165 vs 1,092) is the allele present in nearly every genome.
