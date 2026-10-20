//! Build the pangenome from parsed samples.
//!
//! Homology runs first. Synteny then splits paralogs and fusions and adopts
//! fragments, because those are annotation boundaries rather than different
//! genes. The reference for each finished group is chosen only after that, so
//! a fragment cannot outvote the intact copy it was broken off from.

use std::collections::{HashMap, HashSet};
use std::sync::atomic::{AtomicUsize, Ordering};

use rayon::prelude::*;

use crate::align::fit_short_in_long;
use crate::cluster::cluster_families;
use crate::error::{Error, Result};
use crate::genetic_code::{is_intact, orf_status, translate};
use crate::log::{self, progress};
use crate::model::{Call, Pools};
use crate::reference::{choose_reference, join_fragments, Allele};
use crate::synteny::{
    cluster_facts, compress_clusters, families_bounded_by_neighbors, freeze_contexts,
    length_evidence, membership_loci, merge_fragments, refind, screen_unsupported_families,
    split_fusions, split_paralogs, ClusterFacts, ALLELE_EDGE_BP, UNCONFIRMED_SHORT,
};
use crate::types::{Gene, GeneMember, MembershipAudit, Pangenome, Params, ParsedSample, Quality};

pub fn build(mut samples: Vec<ParsedSample>, params: &Params) -> Result<Pangenome> {
    params.validate()?;
    if samples.is_empty() {
        return Err(Error::new("no samples"));
    }
    samples.sort_by(|a, b| a.name.cmp(&b.name));
    for pair in samples.windows(2) {
        if pair[0].name == pair[1].name {
            return Err(Error::new(format!(
                "two inputs resolved to the sample name '{}'",
                pair[0].name
            )));
        }
    }

    let sample_names: Vec<String> = samples.iter().map(|sample| sample.name.clone()).collect();
    let mut pools = Pools::default();
    let mut calls = Vec::new();
    let mut annotation_ids = Vec::new();
    let mut contigs = Vec::with_capacity(samples.len());
    let mut logged_counts = Vec::with_capacity(samples.len());
    for (sample_index, sample) in samples.iter_mut().enumerate() {
        logged_counts.push((sample.contigs.len(), sample.genes.len()));
        // The contig string moves into the refind table. Names stay on the
        // sample for the membership audit. Gene strings move into the pools,
        // which keep one copy of each distinct sequence.
        contigs.push(
            sample
                .contigs
                .iter_mut()
                .map(|contig| std::mem::take(&mut contig.sequence))
                .collect::<Vec<_>>(),
        );
        for gene in sample.genes.drain(..) {
            if gene.dna.len() < params.min_cds_bp {
                continue;
            }
            let source = if params.audit_members {
                let id = annotation_ids.len() as u32;
                annotation_ids.push(gene.annotation_id);
                Some(id)
            } else {
                None
            };
            calls.push(Call {
                sample: sample_index as u32,
                contig: gene.contig as u32,
                start: gene.start,
                end: gene.end,
                strand: gene.strand,
                dna: pools.intern_dna(gene.dna),
                protein: pools.intern_protein(gene.protein),
                cluster: 0,
                context: 0,
                locus: None,
                source,
                refound: false,
                fusion_piece: false,
                symbol: gene.symbol,
                product: gene.product,
            });
        }
    }
    if calls.is_empty() {
        return Err(Error::new("no coding sequences passed the minimum length"));
    }
    if log::enabled() {
        for (sample, &(contig_count, gene_count)) in samples.iter().zip(&logged_counts) {
            log::record(&format!(
                "sample\t{}\tcontigs\t{}\tgenes\t{}",
                sample.name, contig_count, gene_count
            ));
        }
    }
    progress!(
        "rustroo: {} samples, {} coding sequences, {} unique proteins",
        sample_names.len(),
        calls.len(),
        pools.protein.len()
    );

    // `pools.protein` is already deduplicated by content (interning hashes on
    // insert), so clustering it directly gives one family id per unique
    // protein. Cloning a protein per call first (235k clones on a 167-sample
    // run, for only ~8.6k distinct proteins) made `cluster_families` redo
    // that same deduplication from scratch over the expanded copy.
    progress!("rustroo: homology clustering");
    let families = cluster_families(
        &pools.protein,
        params.precluster_identity,
        params.family_identity,
        params.length_ratio,
    );
    for call in calls.iter_mut() {
        call.cluster = families[call.protein as usize];
    }
    compress_clusters(&mut calls);
    let family_count = calls
        .iter()
        .map(|call| call.cluster)
        .max()
        .map_or(0, |id| id + 1);
    progress!("rustroo: {family_count} families after homology");

    progress!("rustroo: splitting fusions");
    let fusions = split_fusions(&mut calls, &mut pools);
    progress!("rustroo: split {fusions} fusion annotations");
    compress_clusters(&mut calls);
    progress!(
        "rustroo: {} groups after fusion repair",
        dense_family_count(&calls)
    );
    freeze_contexts(&mut calls, &pools);
    let contig_lengths: Vec<Vec<usize>> = contigs
        .iter()
        .map(|sample| sample.iter().map(String::len).collect())
        .collect();
    for round in 1..=4 {
        progress!("rustroo: synteny round {round}");
        let before = dense_family_count(&calls);
        let split = split_paralogs(&mut calls, &pools);
        let merged = merge_fragments(&mut calls, &pools, &contig_lengths);
        compress_clusters(&mut calls);
        let after = dense_family_count(&calls);
        progress!(
            "rustroo: synteny round {round}: {before} -> {} after {split} paralog splits -> {after} after merges",
            before + split
        );
        // Stable anchors and retained loci make splitting idempotent until a
        // real merge changes the evidence; a split-only pass is already settled.
        if !merged {
            progress!("rustroo: synteny settled on round {round}");
            break;
        }
        if round == 4 {
            progress!(
                "rustroo: synteny reached the 4-round limit while groups were still changing"
            );
        }
    }
    if params.refind {
        progress!("rustroo: searching intergenic gaps");
        let found = refind(&mut calls, &mut pools, &contigs);
        progress!("rustroo: refound {found} genes in intergenic gaps");
        if found > 0 {
            let split = split_paralogs(&mut calls, &pools);
            compress_clusters(&mut calls);
            progress!(
                "rustroo: {split} paralog splits after refinding, {} groups",
                dense_family_count(&calls)
            );
        }
    }
    // Gap search is the last reader of contig bases. Lengths were already copied.
    drop(contigs);

    // Contig-end tails and wholly invalid families are annotation, not genes.
    // This runs after repair so a fragment that synteny adopted is not judged
    // on its own, and a refound ORF can still supply the intact allele.
    screen_unsupported_families(&mut calls, &pools, sample_names.len());
    if calls.is_empty() {
        return Err(Error::new(
            "no supported gene family remained after removing invalid and contig-end families",
        ));
    }
    compress_clusters(&mut calls);

    let bounded = families_bounded_by_neighbors(&mut calls, &pools);
    let mut groups: HashMap<u32, Vec<usize>> = HashMap::new();
    for (index, call) in calls.iter().enumerate() {
        groups.entry(call.cluster).or_default().push(index);
    }
    let reps = cluster_reps(&groups, &calls, &pools);
    let locus_ids = membership_loci(&calls, &pools);
    let mut genes = Vec::with_capacity(groups.len());
    let mut ids: Vec<u32> = groups.keys().copied().collect();
    ids.sort_by(|&a, &b| groups[&b].len().cmp(&groups[&a].len()).then(a.cmp(&b)));
    let total = ids.len();
    progress!("rustroo: choosing references for {total} groups");
    // Everything about one group's reference - its alleles, the fusion/
    // extension check, presence and length stats - only reads `calls`,
    // `pools`, and `reps`, none of which are mutated here, so every group is
    // independent of every other and the search can run across threads.
    // Only turning a symbol into a final, de-duplicated gene id has to stay
    // sequential, since a later group's name depends on which names earlier
    // groups (in this exact, already-sorted order) have already taken.
    let progress = AtomicUsize::new(0);
    let partials: Vec<Option<PartialGene>> = ids
        .par_iter()
        .map(|&cluster| {
            let partial = compute_partial_gene(
                cluster,
                &groups[&cluster],
                &calls,
                &pools,
                &reps,
                &locus_ids,
                sample_names.len(),
                bounded.contains(&cluster),
                &contig_lengths,
            );
            let done = progress.fetch_add(1, Ordering::Relaxed) + 1;
            if total >= 100 && (done == 1 || done == total || done % 250 == 0) {
                progress!("rustroo: reference {done}/{total}");
            }
            partial
        })
        .collect();
    let mut used_names = HashSet::new();
    let mut name_of_cluster = HashMap::new();
    let mut anonymous = 0u32;
    for (&cluster, partial) in ids.iter().zip(partials) {
        let Some(partial) = partial else { continue };
        let id = unique_name(&partial.symbol, &mut anonymous, &mut used_names);
        name_of_cluster.insert(cluster, id.clone());
        genes.push(Gene {
            id,
            product: partial.product,
            reference_dna: partial.reference_dna,
            intact: partial.intact,
            selection: partial.selection,
            samples_present: partial.samples_present,
            length_evidence: partial.length_evidence,
            samples_with_exact_reference: partial.samples_with_exact_reference,
            distinct_alleles: partial.distinct_alleles,
            intact_alleles: partial.intact_alleles,
            modal_end_bp: partial.modal_end_bp,
            modal_end_samples: partial.modal_end_samples,
            extension_samples: partial.extension_samples,
            presence: partial.presence,
            lengths: partial.lengths,
            support: SUPPORT_ACCESSORY,
        });
    }
    genes.sort_by(|a, b| {
        b.samples_present
            .cmp(&a.samples_present)
            .then(a.id.cmp(&b.id))
    });
    let mut loci: HashMap<String, Vec<Locus>> = HashMap::new();
    for call in &calls {
        let Some(name) = name_of_cluster.get(&call.cluster) else {
            continue;
        };
        loci.entry(name.clone()).or_default().push(Locus {
            sample: call.sample,
            contig: call.contig,
            start: call.start,
            end: call.end,
            strand: call.strand,
        });
    }
    let mut removed = apply_quality(&mut genes, params.quality, &loci);
    let cluster_of_gene: HashMap<String, u32> = name_of_cluster
        .iter()
        .map(|(&cluster, id)| (id.clone(), cluster))
        .collect();
    assign_support(
        &mut genes,
        &cluster_of_gene,
        &cluster_facts(&calls, &pools),
        sample_names.len(),
    );
    if params.quality == Quality::Strict {
        removed += drop_unsupported(&mut genes);
    }
    if params.quality != Quality::Standard {
        progress!(
            "rustroo: quality {label} removed {removed} genes",
            label = quality_label(params.quality)
        );
        if genes.is_empty() {
            return Err(Error::new("no gene remained after the quality filter"));
        }
    }
    let membership_audit = if params.audit_members {
        let gene_indexes: HashMap<&str, usize> = genes
            .iter()
            .enumerate()
            .map(|(index, gene)| (gene.id.as_str(), index))
            .collect();
        let mut members = Vec::with_capacity(calls.len());
        for (call, &locus) in calls.iter().zip(&locus_ids) {
            let Some(name) = name_of_cluster.get(&call.cluster) else {
                continue;
            };
            let Some(&gene_index) = gene_indexes.get(name.as_str()) else {
                continue;
            };
            members.push(GeneMember {
                gene: gene_index,
                sample: call.sample,
                contig: call.contig,
                annotation_id: call
                    .source
                    .map(|id| annotation_ids[id as usize].clone())
                    .unwrap_or_default(),
                start: call.start,
                end: call.end,
                strand: call.strand,
                locus,
                dna: call.dna,
                origin: if call.refound {
                    "refound"
                } else if call.fusion_piece {
                    "fusion_piece"
                } else {
                    "annotated"
                },
            });
        }
        Some(MembershipAudit {
            contig_names: samples
                .iter()
                .map(|sample| {
                    sample
                        .contigs
                        .iter()
                        .map(|contig| contig.name.clone())
                        .collect()
                })
                .collect(),
            members,
            sequences: pools.dna,
        })
    } else {
        None
    };
    progress!("rustroo: {} genes", genes.len());
    if log::enabled() {
        log::record("gene\tid\tproduct\treference_bp\tintact\tselection\tsamples_present\tlength_evidence\tsupport\tmodal_end_bp\tmodal_end_samples\textension_samples\tsamples_with_exact_reference\tdistinct_alleles\tintact_alleles");
        for gene in &genes {
            log::record(&format!(
                "gene\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                gene.id,
                log_field(&gene.product),
                gene.reference_dna.len(),
                if gene.intact { "yes" } else { "no" },
                gene.selection,
                gene.samples_present,
                gene.length_evidence,
                gene.support,
                gene.modal_end_bp,
                gene.modal_end_samples,
                gene.extension_samples,
                gene.samples_with_exact_reference,
                gene.distinct_alleles,
                gene.intact_alleles,
            ));
        }
        let mut tally: HashMap<&str, usize> = HashMap::new();
        for gene in &genes {
            *tally.entry(gene.selection).or_default() += 1;
        }
        let mut labels: Vec<_> = tally.into_iter().collect();
        labels.sort_unstable();
        for (label, count) in labels {
            log::record(&format!("summary\tselection\t{label}\t{count}"));
        }
        let mut support_tally: HashMap<&str, usize> = HashMap::new();
        for gene in &genes {
            *support_tally.entry(gene.support).or_default() += 1;
        }
        let mut support_labels: Vec<_> = support_tally.into_iter().collect();
        support_labels.sort_unstable();
        for (label, count) in support_labels {
            log::record(&format!("summary\tsupport\t{label}\t{count}"));
        }
        log::flush();
    }
    Ok(Pangenome {
        sample_names,
        genes,
        membership_audit,
    })
}

fn quality_label(quality: Quality) -> &'static str {
    match quality {
        Quality::Standard => "standard",
        Quality::High => "high",
        Quality::Strict => "strict",
    }
}

/// `high` removes a family with no intact copy, and a short family the
/// collection does not support by frequency or by neighboring genes.
/// `strict` then removes a second reference of a gene already kept. A nested
/// or deep fragment, or any shorter piece, can match at 80% identity, and a
/// copy at 90%. A repeat offset does not keep the second copy. The same allele
/// in other genomes is folded into the longer gene so those genomes stay
/// present. The reference stays the longer finished allele. A paralogue below
/// 90% identity at nearly the same length stays. Support labels are applied
/// afterwards. `strict` then drops only `unsupported`.
fn apply_quality(
    genes: &mut Vec<Gene>,
    quality: Quality,
    loci: &HashMap<String, Vec<Locus>>,
) -> usize {
    if quality == Quality::Standard {
        return 0;
    }
    let before = genes.len();
    let label = quality_label(quality);
    genes.retain(|gene| {
        let reject = !gene.intact || gene.length_evidence == UNCONFIRMED_SHORT;
        if reject && log::enabled() {
            log::record(&format!(
                "quality\t{label}\tdrop\t{}\tintact\t{}\tlength_evidence\t{}",
                gene.id,
                if gene.intact { "yes" } else { "no" },
                gene.length_evidence,
            ));
        }
        !reject
    });
    if quality == Quality::Strict {
        let redundant = redundant_fragments(genes, loci);
        let mut index = 0;
        genes.retain(|_| {
            let reject = redundant.contains(&index);
            index += 1;
            !reject
        });
    }
    before - genes.len()
}

const SUPPORT_CORE: &str = "core";
const SUPPORT_ACCESSORY: &str = "accessory";
const SUPPORT_UNSUPPORTED: &str = "unsupported";
const SUPPORT_REFOUND: &str = "refound";

/// A gene in at least 95% of genomes is core, whether or not it was annotated.
/// Below that, a gene with no annotated call is refound. An annotated gene
/// whose copies share one neighborhood is accessory. Anything else is
/// unsupported.
fn support_label(
    samples_present: u32,
    n_samples: usize,
    annotated: bool,
    neighbors_agree: bool,
) -> &'static str {
    if n_samples > 0 && u64::from(samples_present) * 20 >= n_samples as u64 * 19 {
        return SUPPORT_CORE;
    }
    if !annotated {
        return SUPPORT_REFOUND;
    }
    if neighbors_agree {
        SUPPORT_ACCESSORY
    } else {
        SUPPORT_UNSUPPORTED
    }
}

fn assign_support(
    genes: &mut [Gene],
    cluster_of_gene: &HashMap<String, u32>,
    facts: &HashMap<u32, ClusterFacts>,
    n_samples: usize,
) {
    for gene in genes {
        let fact = cluster_of_gene
            .get(&gene.id)
            .and_then(|cluster| facts.get(cluster))
            .copied();
        let (annotated, neighbors_agree) = match fact {
            Some(fact) => (fact.annotated, fact.neighbors_agree),
            None => (true, false),
        };
        gene.support = support_label(gene.samples_present, n_samples, annotated, neighbors_agree);
    }
}

fn drop_unsupported(genes: &mut Vec<Gene>) -> usize {
    let before = genes.len();
    genes.retain(|gene| {
        let reject = gene.support == SUPPORT_UNSUPPORTED;
        if reject && log::enabled() {
            log::record(&format!(
                "quality\tstrict\tdrop\t{}\tsupport\t{SUPPORT_UNSUPPORTED}",
                gene.id,
            ));
        }
        !reject
    });
    before - genes.len()
}

struct Locus {
    sample: u32,
    contig: u32,
    start: u32,
    end: u32,
    strand: i8,
}

/// Below this, a nested ORF can be a different gene that only shares a domain.
const NESTED_PIECE_IDENTITY: f64 = 0.80;
/// A fragment this close, covering less than 90% of the longer protein, is
/// the same peptide.
const DEEP_PIECE_IDENTITY: f64 = 0.80;
/// A piece covering less than 98% of the longer protein. Same floor as a
/// deep fragment. A paralogue below 90% identity at nearly full length stays.
const PIECE_IDENTITY: f64 = 0.80;
/// Minimum fraction of the shorter protein that must align.
const COVERED: f64 = 0.90;
/// A second copy this close is the same gene, including at full length.
/// A paralogue below this, at nearly the same length, stays.
const NEAR_COPY_IDENTITY: f64 = 0.90;

/// Indexes of genes that add no presence. The copy present in fewer genomes
/// is the one removed. Equal support keeps the earlier gene.
fn redundant_fragments(genes: &mut [Gene], loci: &HashMap<String, Vec<Locus>>) -> HashSet<usize> {
    let proteins: Vec<Vec<u8>> = genes
        .iter()
        .map(|gene| translate(&gene.reference_dna))
        .collect();
    // Five residues is specific enough to connect a fragment with its long
    // copy. A tetrapeptide shared by dozens of unrelated genes is not.
    let mut postings: HashMap<[u8; 5], Vec<usize>> = HashMap::new();
    let mut kmer_counts = vec![0u32; genes.len()];
    for (index, protein) in proteins.iter().enumerate() {
        if protein.len() < 5 {
            continue;
        }
        let mut seen = HashSet::new();
        for window in protein.windows(5) {
            let key = [window[0], window[1], window[2], window[3], window[4]];
            if seen.insert(key) {
                postings.entry(key).or_default().push(index);
                kmer_counts[index] += 1;
            }
        }
    }
    // A repeated family shares every k-mer, so the posting is large. Dropping
    // it hid the long copy from every fragment. Keep the longest members,
    // and the earliest index on a tie, which is the copy the others fold into.
    for hits in postings.values_mut() {
        if hits.len() <= 256 {
            continue;
        }
        hits.sort_by(|&left, &right| {
            proteins[right]
                .len()
                .cmp(&proteins[left].len())
                .then(
                    genes[right]
                        .samples_present
                        .cmp(&genes[left].samples_present),
                )
                .then(left.cmp(&right))
        });
        hits.truncate(256);
    }
    let mut order: Vec<usize> = (0..genes.len()).collect();
    order.sort_by(|&left, &right| {
        proteins[left]
            .len()
            .cmp(&proteins[right].len())
            .then(left.cmp(&right))
    });
    let mut redundant = HashSet::new();
    for short in order {
        let needed = kmer_counts[short];
        if needed == 0 || redundant.contains(&short) {
            continue;
        }
        let mut scores: HashMap<usize, u32> = HashMap::new();
        let protein = &proteins[short];
        let mut seen = HashSet::new();
        for window in protein.windows(5) {
            let key = [window[0], window[1], window[2], window[3], window[4]];
            if !seen.insert(key) {
                continue;
            }
            let Some(hits) = postings.get(&key) else {
                continue;
            };
            for &long in hits {
                if long != short && is_covered_by(short, long, &proteins, genes) {
                    *scores.entry(long).or_default() += 1;
                }
            }
        }
        // One mismatch wipes five 5-mers, so an 80% piece keeps only a small
        // share of them. The alignment, not this count, decides the fold.
        let mut ranked: Vec<(usize, u32)> = scores
            .into_iter()
            .filter(|&(_, score)| score * 10 >= needed)
            .collect();
        // Longest first, so a pile of equal copies cannot crowd the real parent
        // out of the alignments.
        ranked.sort_by(|a, b| {
            proteins[b.0]
                .len()
                .cmp(&proteins[a.0].len())
                .then(b.1.cmp(&a.1))
                .then(a.0.cmp(&b.0))
        });
        ranked.truncate(32);
        let mut longer_subset: Option<usize> = None;
        let mut equal_subset: Option<usize> = None;
        let mut longer_exclusive: Vec<usize> = Vec::new();
        let mut equal_exclusive: Vec<usize> = Vec::new();
        for (long, _) in ranked {
            if redundant.contains(&long) {
                continue;
            }
            let Some(fit) = fit_short_in_long(protein, &proteins[long]) else {
                continue;
            };
            let nested = calls_nested(loci, &genes[short].id, &genes[long].id);
            let longer = proteins[long].len() > protein.len();
            if presence_covered(&genes[short].presence, &genes[long].presence) {
                if !covered_duplicate(&fit, proteins[long].len(), nested) {
                    continue;
                }
                if longer {
                    longer_subset.get_or_insert(long);
                } else {
                    equal_subset.get_or_insert(long);
                }
            } else if exclusive_allele(&fit, proteins[long].len()) {
                if longer {
                    longer_exclusive.push(long);
                } else {
                    equal_exclusive.push(long);
                }
            }
        }
        // One longer gene explains the fragment. Several longer copies of
        // that same gene do too. Two longer genes that are not the same
        // reference leave the fragment in place. Equal-length copies are
        // one reference and fold even when several of them match.
        let parent = if let Some(long) = longer_subset {
            Some(long)
        } else if let Some(long) = one_gene(&proteins, &longer_exclusive) {
            Some(long)
        } else if !longer_exclusive.is_empty() {
            None
        } else if let Some(long) = equal_subset {
            Some(long)
        } else {
            equal_exclusive
                .into_iter()
                .min_by_key(|&index| (std::cmp::Reverse(genes[index].samples_present), index))
        };
        if let Some(long) = parent {
            if log::enabled() {
                log::record(&format!(
                    "quality\tstrict\tdrop\t{}\tcovered_by\t{}",
                    genes[short].id, genes[long].id,
                ));
            }
            absorb_at(genes, long, short);
            redundant.insert(short);
        }
    }
    redundant
}

/// The longest of `parents` when every parent is that same reference.
/// Two parents that fail the duplicate test are different genes.
fn one_gene(proteins: &[Vec<u8>], parents: &[usize]) -> Option<usize> {
    let best = parents
        .iter()
        .copied()
        .min_by_key(|&index| (std::cmp::Reverse(proteins[index].len()), index))?;
    let same = parents.iter().copied().all(|other| {
        if other == best {
            return true;
        }
        let (shorter, longer) = if proteins[other].len() <= proteins[best].len() {
            (other, best)
        } else {
            (best, other)
        };
        fit_short_in_long(&proteins[shorter], &proteins[longer])
            .is_some_and(|fit| exclusive_allele(&fit, proteins[longer].len()))
    });
    same.then_some(best)
}

fn absorb_at(genes: &mut [Gene], long: usize, short: usize) {
    if long == short {
        return;
    }
    if long < short {
        let (head, tail) = genes.split_at_mut(short);
        absorb(&mut head[long], &tail[0]);
    } else {
        let (head, tail) = genes.split_at_mut(long);
        absorb(&mut tail[0], &head[short]);
    }
}

fn absorb(keeper: &mut Gene, extra: &Gene) {
    if keeper.presence.len() != extra.presence.len() {
        return;
    }
    // Presence moves onto the keeper. Its sequence stays: the shorter allele
    // fails coverage of a longer reference of the same gene.
    for index in 0..keeper.presence.len() {
        if extra.presence[index] && !keeper.presence[index] {
            keeper.presence[index] = true;
            if index < keeper.lengths.len() && keeper.lengths.len() == extra.lengths.len() {
                keeper.lengths[index] = extra.lengths[index];
            }
        }
    }
    keeper.samples_present = keeper.presence.iter().filter(|&&present| present).count() as u32;
}

/// `long` already explains `short`: it is a longer protein, or the same length
/// and present in more genomes. Equal length and equal support keep the
/// earlier gene, so a mutual match cannot remove both.
fn is_covered_by(short: usize, long: usize, proteins: &[Vec<u8>], genes: &[Gene]) -> bool {
    match proteins[long].len().cmp(&proteins[short].len()) {
        std::cmp::Ordering::Greater => true,
        std::cmp::Ordering::Less => false,
        std::cmp::Ordering::Equal => {
            genes[long].samples_present > genes[short].samples_present
                || (genes[long].samples_present == genes[short].samples_present && long < short)
        }
    }
}

/// `Some(true)` when every short call in a shared genome lies inside a long
/// call. `Some(false)` when a short call sits beside the long gene, or the
/// coordinates are missing from a map that otherwise has them. `None` when
/// the map is empty, which keeps the sequence-only rule for tests.
fn calls_nested(loci: &HashMap<String, Vec<Locus>>, short_id: &str, long_id: &str) -> Option<bool> {
    if loci.is_empty() {
        return None;
    }
    let (Some(short), Some(long)) = (loci.get(short_id), loci.get(long_id)) else {
        return Some(false);
    };
    let mut shared = false;
    for short_call in short {
        let same_sample: Vec<&Locus> = long
            .iter()
            .filter(|call| call.sample == short_call.sample)
            .collect();
        if same_sample.is_empty() {
            continue;
        }
        shared = true;
        let inside = same_sample.iter().any(|call| {
            call.contig == short_call.contig
                && call.strand == short_call.strand
                && short_call.start >= call.start
                && short_call.end <= call.end
        });
        if !inside {
            return Some(false);
        }
    }
    Some(shared)
}

fn covered_duplicate(fit: &crate::align::Fit, long_len: usize, nested: Option<bool>) -> bool {
    if fit.coverage + f64::EPSILON < COVERED {
        return false;
    }
    // A second diagonal is a repeat of the same peptide, not a second gene,
    // once this diagonal already covers the shorter protein.
    let span = fit.long_end.saturating_sub(fit.long_start);
    let deep = span * 10 < long_len * 9;
    let piece = span * 50 < long_len * 49;
    let identity = fit.identity + f64::EPSILON;
    if nested == Some(true) && piece && identity >= NESTED_PIECE_IDENTITY {
        return true;
    }
    (deep && identity >= DEEP_PIECE_IDENTITY)
        || (piece && identity >= PIECE_IDENTITY)
        || identity >= NEAR_COPY_IDENTITY
}

/// Same gene in genomes the longer copy does not already cover. A near-full
/// paralogue below 90% identity is not an allele split.
fn exclusive_allele(fit: &crate::align::Fit, long_len: usize) -> bool {
    if fit.coverage + f64::EPSILON < COVERED {
        return false;
    }
    let span = fit.long_end.saturating_sub(fit.long_start);
    let deep = span * 10 < long_len * 9;
    let piece = span * 50 < long_len * 49;
    let identity = fit.identity + f64::EPSILON;
    identity >= NEAR_COPY_IDENTITY
        || (piece && identity >= PIECE_IDENTITY)
        || (deep && identity >= DEEP_PIECE_IDENTITY)
}

fn presence_covered(short: &[bool], long: &[bool]) -> bool {
    short.len() == long.len()
        && short
            .iter()
            .zip(long)
            .all(|(&has_short, &has_long)| !has_short || has_long)
}

fn log_field(value: &str) -> String {
    value
        .chars()
        .map(|character| match character {
            '\t' | '\n' | '\r' => ' ',
            _ => character,
        })
        .collect()
}

fn dense_family_count(calls: &[Call]) -> usize {
    calls
        .iter()
        .map(|call| call.cluster as usize)
        .max()
        .map_or(0, |id| id + 1)
}

/// Everything about one group's chosen reference except its final, unique
/// gene id - that part alone has to be assigned sequentially afterward, so
/// it is deliberately not computed here. Pulled out of `build` so the search
/// for each group can be an independent unit of parallel work.
struct PartialGene {
    symbol: String,
    product: String,
    reference_dna: String,
    intact: bool,
    selection: &'static str,
    presence: Vec<bool>,
    lengths: Vec<u32>,
    samples_present: u32,
    length_evidence: &'static str,
    samples_with_exact_reference: u32,
    distinct_alleles: u32,
    intact_alleles: u32,
    modal_end_bp: u32,
    modal_end_samples: u32,
    extension_samples: u32,
}

fn compute_partial_gene(
    cluster: u32,
    indexes: &[usize],
    calls: &[Call],
    pools: &Pools,
    reps: &HashMap<u32, Vec<u8>>,
    locus_ids: &[u32],
    sample_count: usize,
    between_neighbors: bool,
    contig_lengths: &[Vec<usize>],
) -> Option<PartialGene> {
    let alleles = alleles_of(indexes, calls, pools, contig_lengths, locus_ids);
    let joined = joined_allele(indexes, calls, pools, locus_ids);
    let mut with_join = alleles.clone();
    if let Some(dna) = joined {
        let longer = alleles.iter().all(|allele| dna.len() > allele.dna.len());
        if longer {
            with_join.push(Allele {
                dna,
                count: 0,
                synthetic: true,
                finished: 0,
            });
        }
    }
    let choice = choose_reference(&with_join, reps, cluster, sample_count)?;
    let mut presence = vec![false; sample_count];
    let mut lengths = vec![0u32; sample_count];
    let mut intervals: HashMap<(u32, u32), Vec<(u32, u32)>> = HashMap::new();
    let mut exact = HashSet::new();
    for &index in indexes {
        let call = &calls[index];
        presence[call.sample as usize] = true;
        intervals
            .entry((call.sample, call.contig))
            .or_default()
            .push((call.start, call.end));
        if pools.dna[call.dna as usize] == choice.dna {
            exact.insert(call.sample);
        }
    }
    for ((sample, _), spans) in intervals {
        lengths[sample as usize] += covered_bases(&spans);
    }
    let distinct = alleles.len() as u32;
    let intact_alleles = alleles
        .iter()
        .filter(|allele| is_intact(&allele.dna))
        .count() as u32;
    let samples_present = presence.iter().filter(|present| **present).count() as u32;
    Some(PartialGene {
        symbol: preferred_symbol(indexes, calls),
        product: preferred_product(indexes, calls),
        reference_dna: choice.dna.clone(),
        intact: choice.intact,
        selection: choice.selection,
        samples_present,
        length_evidence: length_evidence(
            choice.dna.len(),
            samples_present,
            sample_count,
            between_neighbors,
        ),
        samples_with_exact_reference: exact.len() as u32,
        distinct_alleles: distinct,
        intact_alleles,
        modal_end_bp: choice.modal_end_bp,
        modal_end_samples: choice.modal_end_samples,
        extension_samples: choice.extension_samples,
        presence,
        lengths,
    })
}

fn alleles_of(
    indexes: &[usize],
    calls: &[Call],
    pools: &Pools,
    contig_lengths: &[Vec<usize>],
    locus_ids: &[u32],
) -> Vec<Allele> {
    // Two abutting pieces of one locus are a broken gene, not a shorter
    // finished stop. Counting either piece would let the break define the
    // reference and throw away the intact copy.
    let mut parts_of_locus: HashMap<u32, u32> = HashMap::new();
    for &index in indexes {
        *parts_of_locus.entry(locus_ids[index]).or_default() += 1;
    }
    let mut samples_for_dna: HashMap<u32, HashSet<u32>> = HashMap::new();
    let mut finished_for_dna: HashMap<u32, HashSet<u32>> = HashMap::new();
    let mut closed: HashMap<u32, bool> = HashMap::new();
    for &index in indexes {
        let call = &calls[index];
        samples_for_dna
            .entry(call.dna)
            .or_default()
            .insert(call.sample);
        let closed = *closed.entry(call.dna).or_insert_with(|| {
            let status = orf_status(&pools.dna[call.dna as usize]);
            status.is_intact() && status.terminal_stop
        });
        let whole = parts_of_locus.get(&locus_ids[index]).copied().unwrap_or(0) == 1;
        // One internal whole copy is enough for that sample to vouch for the
        // end. A second copy parked on a contig break does not cancel it.
        if closed && whole && !contig_clipped(call, contig_lengths) {
            finished_for_dna
                .entry(call.dna)
                .or_default()
                .insert(call.sample);
        }
    }
    let mut alleles: Vec<Allele> = samples_for_dna
        .into_iter()
        .map(|(dna, samples)| Allele {
            dna: pools.dna[dna as usize].clone(),
            count: samples.len() as u32,
            synthetic: false,
            finished: finished_for_dna
                .get(&dna)
                .map(|samples| samples.len() as u32)
                .unwrap_or(0),
        })
        .collect();
    alleles.sort_by(|a, b| b.count.cmp(&a.count).then(a.dna.cmp(&b.dna)));
    alleles
}

/// A copy this close to either contig boundary may have been clipped. It
/// stays in the family, and it does not get to define where the gene ends.
fn contig_clipped(call: &Call, contig_lengths: &[Vec<usize>]) -> bool {
    let Some(length) = contig_lengths
        .get(call.sample as usize)
        .and_then(|contigs| contigs.get(call.contig as usize))
        .copied()
    else {
        return true;
    };
    call.start as usize <= ALLELE_EDGE_BP
        || length.saturating_sub(call.end as usize) <= ALLELE_EDGE_BP
}

fn joined_allele(
    indexes: &[usize],
    calls: &[Call],
    pools: &Pools,
    locus_ids: &[u32],
) -> Option<String> {
    let mut by_locus: HashMap<(u32, u32, u32), Vec<usize>> = HashMap::new();
    for &index in indexes {
        let call = &calls[index];
        by_locus
            .entry((call.sample, call.contig, locus_ids[index]))
            .or_default()
            .push(index);
    }
    let mut joins: HashMap<String, u32> = HashMap::new();
    for group in by_locus.values() {
        if group.len() < 2 {
            continue;
        }
        let mut ordered = group.clone();
        ordered.sort_by_key(|index| calls[*index].start);
        let strand = calls[ordered[0]].strand;
        if ordered.iter().any(|index| calls[*index].strand != strand) {
            continue;
        }
        if strand < 0 {
            ordered.reverse();
        }
        let mut abutting = true;
        for pair in ordered.windows(2) {
            let gap = if strand >= 0 {
                i64::from(calls[pair[1]].start) - i64::from(calls[pair[0]].end) - 1
            } else {
                i64::from(calls[pair[0]].start) - i64::from(calls[pair[1]].end) - 1
            };
            if !(0..=3).contains(&gap) {
                abutting = false;
                break;
            }
        }
        if !abutting {
            continue;
        }
        let parts: Vec<&str> = ordered
            .iter()
            .map(|index| pools.dna[calls[*index].dna as usize].as_str())
            .collect();
        if let Some(joined) = join_fragments(&parts) {
            *joins.entry(joined).or_default() += 1;
        }
    }
    joins
        .into_iter()
        .max_by(|a, b| a.1.cmp(&b.1).then(b.0.cmp(&a.0)))
        .map(|(dna, _)| dna)
}

fn cluster_reps(
    groups: &HashMap<u32, Vec<usize>>,
    calls: &[Call],
    pools: &Pools,
) -> HashMap<u32, Vec<u8>> {
    let mut reps = HashMap::new();
    for (&cluster, indexes) in groups {
        let protein = indexes
            .iter()
            .map(|index| &pools.protein[calls[*index].protein as usize])
            .max_by_key(|protein| protein.len())
            .cloned()
            .unwrap_or_default();
        reps.insert(cluster, protein);
    }
    reps
}

fn preferred_symbol(indexes: &[usize], calls: &[Call]) -> String {
    let mut votes: HashMap<String, u32> = HashMap::new();
    let mut seen = HashSet::new();
    for &index in indexes {
        let call = &calls[index];
        if !seen.insert((call.sample, call.symbol.clone())) {
            continue;
        }
        if usable_symbol(&call.symbol) {
            *votes.entry(call.symbol.clone()).or_default() += 1;
        }
    }
    votes
        .into_iter()
        .max_by(|a, b| a.1.cmp(&b.1).then(b.0.cmp(&a.0)))
        .map(|(symbol, _)| symbol)
        .unwrap_or_default()
}

fn preferred_product(indexes: &[usize], calls: &[Call]) -> String {
    let mut votes: HashMap<String, u32> = HashMap::new();
    for &index in indexes {
        let product = calls[index].product.trim();
        if product.is_empty() {
            continue;
        }
        *votes.entry(product.to_string()).or_default() += 1;
    }
    votes
        .into_iter()
        .max_by(|a, b| a.1.cmp(&b.1).then(b.0.cmp(&a.0)))
        .map(|(product, _)| product)
        .unwrap_or_else(|| "hypothetical protein".to_string())
}

fn usable_symbol(symbol: &str) -> bool {
    let trimmed = symbol.trim();
    if trimmed.is_empty() {
        return false;
    }
    let lower = trimmed.to_ascii_lowercase();
    !matches!(
        lower.as_str(),
        "hypothetical" | "hypothetical protein" | "unknown" | "uncharacterized protein"
    )
}

fn unique_name(symbol: &str, anonymous: &mut u32, used: &mut HashSet<String>) -> String {
    let mut base = sanitize(symbol);
    if !usable_symbol(&base) {
        *anonymous += 1;
        base = format!("group_{anonymous:05}");
    }
    let mut name = base.clone();
    let mut suffix = 2u32;
    while !used.insert(name.clone()) {
        name = format!("{base}_{suffix}");
        suffix += 1;
    }
    name
}

fn sanitize(symbol: &str) -> String {
    let mut out = String::with_capacity(symbol.len());
    for character in symbol.chars() {
        if character.is_ascii_alphanumeric() || matches!(character, '.' | '_' | '-') {
            out.push(character);
        } else {
            out.push('_');
        }
    }
    out
}

fn covered_bases(intervals: &[(u32, u32)]) -> u32 {
    let mut spans = intervals.to_vec();
    spans.sort_by_key(|span| span.0);
    let mut total = 0u32;
    let mut current: Option<(u32, u32)> = None;
    for (start, end) in spans {
        if end < start {
            continue;
        }
        if let Some((open, close)) = current {
            if start <= close.saturating_add(1) {
                current = Some((open, close.max(end)));
            } else {
                total = total.saturating_add(close - open + 1);
                current = Some((start, end));
            }
        } else {
            current = Some((start, end));
        }
    }
    if let Some((open, close)) = current {
        total = total.saturating_add(close - open + 1);
    }
    total
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genetic_code::{back_translate, translate};
    use crate::types::{Contig, Gene, ParsedGene, ParsedSample};

    fn diverse(length: usize, salt: usize) -> Vec<u8> {
        let alphabet = b"ACDEFGHIKLMNPQRSTVWY";
        let mut protein = Vec::with_capacity(length);
        protein.push(b'M');
        for index in 1..length {
            let mut residue = alphabet[(index * 5 + salt) % alphabet.len()];
            if residue == b'M' {
                residue = b'L';
            }
            protein.push(residue);
        }
        protein
    }

    fn random_protein(length: usize, mut state: u64) -> Vec<u8> {
        let alphabet = b"ACDEFGHIKLMNPQRSTVWY";
        let mut sequence = vec![b'M'];
        for _ in 1..length {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let aa = alphabet[(state % 20) as usize];
            sequence.push(if aa == b'M' { b'L' } else { aa });
        }
        sequence
    }

    fn coding(protein: &[u8], contig: usize, start: u32, symbol: &str) -> ParsedGene {
        let dna = back_translate(protein);
        let end = start + dna.len() as u32 - 1;
        ParsedGene {
            annotation_id: format!("cds_{contig}_{start}"),
            contig,
            start,
            end,
            strand: 1,
            protein: translate(&dna),
            dna,
            symbol: symbol.to_string(),
            product: symbol.to_string(),
        }
    }

    fn sample(name: &str, genes: Vec<ParsedGene>) -> ParsedSample {
        let mut end = 1u32;
        for gene in &genes {
            end = end.max(gene.end + 10);
        }
        ParsedSample {
            name: name.to_string(),
            contigs: vec![Contig {
                name: "contig".to_string(),
                sequence: "N".repeat(end as usize + 10),
            }],
            genes,
        }
    }

    fn labelled(id: &str, intact: bool, length_evidence: &'static str) -> Gene {
        Gene {
            id: id.to_string(),
            product: String::new(),
            reference_dna: "ATGAAATAA".to_string(),
            intact,
            selection: if intact {
                "intact_longest"
            } else {
                "no_intact_copy"
            },
            presence: vec![true],
            lengths: vec![9],
            samples_present: 1,
            length_evidence,
            samples_with_exact_reference: 1,
            distinct_alleles: 1,
            intact_alleles: u32::from(intact),
            modal_end_bp: 9,
            modal_end_samples: 1,
            extension_samples: 0,
            support: SUPPORT_ACCESSORY,
        }
    }

    fn params() -> Params {
        Params {
            min_cds_bp: 30,
            refind: false,
            ..Params::default()
        }
    }

    #[test]
    fn support_label_is_core_then_refound_then_neighborhood() {
        assert_eq!(support_label(19, 20, false, false), SUPPORT_CORE);
        assert_eq!(support_label(18, 20, false, true), SUPPORT_REFOUND);
        assert_eq!(support_label(18, 20, true, true), SUPPORT_ACCESSORY);
        assert_eq!(support_label(18, 20, true, false), SUPPORT_UNSUPPORTED);
        assert_eq!(support_label(1, 1, true, false), SUPPORT_CORE);
    }

    #[test]
    fn strict_drops_only_an_unsupported_gene() {
        let left = random_protein(90, 1);
        let middle = random_protein(90, 2);
        let right = random_protein(90, 3);
        let private = random_protein(90, 11);
        let wander = random_protein(90, 17);
        let mut inputs = Vec::new();
        for index in 0..60 {
            inputs.push(sample(
                &format!("g{index:02}"),
                vec![
                    coding(&left, 0, 1_000, "left"),
                    coding(&middle, 0, 2_000, "middle"),
                    coding(&right, 0, 3_000, "right"),
                ],
            ));
        }
        inputs.push(sample(
            "acc",
            vec![
                coding(&left, 0, 1_000, "left"),
                coding(&private, 0, 2_000, "private"),
                coding(&right, 0, 3_000, "right"),
            ],
        ));
        inputs.push(sample(
            "u1",
            vec![
                coding(&left, 0, 1_000, "left"),
                coding(&wander, 0, 2_000, "wander"),
                coding(&middle, 0, 3_000, "middle"),
            ],
        ));
        inputs.push(sample(
            "u2",
            vec![
                coding(&middle, 0, 1_000, "middle"),
                coding(&wander, 0, 2_000, "wander"),
                coding(&right, 0, 3_000, "right"),
            ],
        ));
        let labelled = |genes: &[Gene], id: &str| {
            genes
                .iter()
                .find(|gene| gene.id == id)
                .map(|gene| gene.support)
        };
        let standard = build(inputs.clone(), &params()).unwrap();
        assert_eq!(labelled(&standard.genes, "middle"), Some(SUPPORT_CORE));
        assert_eq!(
            labelled(&standard.genes, "private"),
            Some(SUPPORT_ACCESSORY)
        );
        assert_eq!(
            labelled(&standard.genes, "wander"),
            Some(SUPPORT_UNSUPPORTED)
        );
        let mut high = params();
        high.quality = Quality::High;
        let high = build(inputs.clone(), &high).unwrap();
        assert_eq!(labelled(&high.genes, "wander"), Some(SUPPORT_UNSUPPORTED));
        assert_eq!(high.genes.len(), standard.genes.len());
        let mut strict = params();
        strict.quality = Quality::Strict;
        let strict = build(inputs, &strict).unwrap();
        assert_eq!(labelled(&strict.genes, "wander"), None);
        assert_eq!(labelled(&strict.genes, "middle"), Some(SUPPORT_CORE));
        assert_eq!(labelled(&strict.genes, "private"), Some(SUPPORT_ACCESSORY));
        assert_eq!(strict.genes.len() + 1, standard.genes.len());
    }

    #[test]
    fn membership_audit_does_not_change_references_or_grouping() {
        let protein = diverse(40, 4);
        let other = diverse(40, 7);
        let inputs = vec![
            sample(
                "b",
                vec![
                    coding(&other, 0, 100, "extra"),
                    coding(&protein, 0, 300, "gene"),
                ],
            ),
            sample(
                "a",
                vec![
                    coding(&other, 0, 100, "extra"),
                    coding(&protein, 0, 400, "gene"),
                ],
            ),
        ];
        let normal = build(inputs.clone(), &params()).unwrap();
        assert!(normal.membership_audit.is_none());
        let mut audit_params = params();
        audit_params.audit_members = true;
        let audited = build(inputs, &audit_params).unwrap();
        assert_eq!(normal.sample_names, audited.sample_names);
        assert_eq!(normal.genes.len(), audited.genes.len());
        for (a, b) in normal.genes.iter().zip(&audited.genes) {
            assert_eq!(a.id, b.id);
            assert_eq!(a.reference_dna, b.reference_dna);
            assert_eq!(a.presence, b.presence);
            assert_eq!(a.lengths, b.lengths);
            assert_eq!(a.selection, b.selection);
        }
        let audit = audited.membership_audit.as_ref().unwrap();
        assert_eq!(audit.members.len(), 4);
        for member in &audit.members {
            let extra = member.start == 100;
            assert_eq!(
                audited.genes[member.gene].id,
                if extra { "extra" } else { "gene" }
            );
            let expected = if extra {
                100
            } else if audited.sample_names[member.sample as usize] == "a" {
                400
            } else {
                300
            };
            assert_eq!(member.start, expected);
            assert_eq!(member.annotation_id, format!("cds_0_{expected}"));
            assert_eq!(member.origin, "annotated");
            assert_eq!(
                audit.sequences[member.dna as usize],
                back_translate(if extra { &other } else { &protein })
            );
        }
        assert_ne!(audit.members[0].locus, audit.members[1].locus);
    }

    fn shared_triplet(name: &str, left: &[u8], middle: &[u8], right: &[u8]) -> ParsedSample {
        sample(
            name,
            vec![
                coding(left, 0, 1, "flkL"),
                coding(middle, 0, 400, "core"),
                coding(right, 0, 900, "flkR"),
            ],
        )
    }

    #[test]
    fn contig_end_singleton_is_removed_and_an_internal_rare_gene_stays() {
        let left = diverse(36, 8);
        let core = diverse(40, 1);
        let right = diverse(36, 15);
        let tail = diverse(36, 21);
        let embedded = diverse(36, 28);
        let mut samples: Vec<_> = (0..6)
            .map(|index| shared_triplet(&format!("s{index}"), &left, &core, &right))
            .collect();
        let mut with_tail = shared_triplet("tail", &left, &core, &right);
        let flank_end = with_tail.genes[2].end;
        with_tail
            .genes
            .push(coding(&tail, 0, flank_end + 40, "tail"));
        with_tail.contigs[0].sequence = "N".repeat(with_tail.genes[3].end as usize + 20);
        samples.push(with_tail);
        let mut with_insert = shared_triplet("insert", &left, &core, &right);
        with_insert
            .genes
            .insert(1, coding(&embedded, 0, 200, "rare"));
        samples.push(with_insert);

        let built = build(samples, &params()).unwrap();
        let ids: Vec<_> = built.genes.iter().map(|gene| gene.id.as_str()).collect();
        assert!(ids.contains(&"core") && ids.contains(&"flkL") && ids.contains(&"flkR"));
        assert!(
            ids.contains(&"rare"),
            "a gene between two neighbors is evidence, not a tail"
        );
        assert!(
            !ids.contains(&"tail"),
            "a private gene at a contig end is not a family"
        );
        let rare = built.genes.iter().find(|gene| gene.id == "rare").unwrap();
        assert_eq!(rare.samples_present, 1);
        assert!(rare.intact);
        assert_eq!(rare.length_evidence, crate::synteny::CONSERVED_SHORT);
        let core = built.genes.iter().find(|gene| gene.id == "core").unwrap();
        assert_eq!(core.length_evidence, crate::synteny::CONSERVED_SHORT);
    }

    #[test]
    fn a_private_contig_of_unsupported_genes_is_removed() {
        let left = diverse(36, 8);
        let core = diverse(40, 1);
        let right = diverse(36, 15);
        let mut samples: Vec<_> = (0..4)
            .map(|index| shared_triplet(&format!("s{index}"), &left, &core, &right))
            .collect();
        samples.push(sample(
            "island",
            vec![
                coding(&diverse(36, 31), 0, 1, "c1"),
                coding(&diverse(36, 32), 0, 400, "c2"),
                coding(&diverse(36, 33), 0, 800, "c3"),
            ],
        ));
        let built = build(samples, &params()).unwrap();
        assert!(built
            .genes
            .iter()
            .all(|gene| !matches!(gene.id.as_str(), "c1" | "c2" | "c3")));
        assert_eq!(built.genes.len(), 3);
    }

    #[test]
    fn a_gap_isolated_singleton_in_the_middle_of_a_contig_stays() {
        let left = diverse(36, 8);
        let core = diverse(40, 1);
        let right = diverse(36, 15);
        let isolated = diverse(36, 22);
        let mut samples: Vec<_> = (0..4)
            .map(|index| shared_triplet(&format!("s{index}"), &left, &core, &right))
            .collect();
        samples.push(sample(
            "middle",
            vec![
                coding(&left, 0, 1, "flkL"),
                coding(&core, 0, 400, "core"),
                coding(&isolated, 0, 20_000, "island"),
                coding(&right, 0, 40_000, "flkR"),
            ],
        ));
        let built = build(samples.clone(), &params()).unwrap();
        let gene = built
            .genes
            .iter()
            .find(|gene| gene.id == "island")
            .expect("gap-isolated gene");
        assert_eq!(gene.samples_present, 1);
        assert_eq!(gene.length_evidence, crate::synteny::UNCONFIRMED_SHORT);
        let mut high = params();
        high.quality = Quality::High;
        let filtered = build(samples, &high).unwrap();
        assert!(filtered.genes.iter().all(|gene| gene.id != "island"));
        assert!(filtered.genes.iter().any(|gene| gene.id == "core"));
        assert!(filtered.genes.iter().all(|gene| {
            gene.intact && gene.length_evidence != crate::synteny::UNCONFIRMED_SHORT
        }));
    }

    #[test]
    fn high_quality_keeps_an_intact_conserved_gene_and_drops_a_broken_one() {
        use crate::synteny::{CONSERVED_SHORT, LENGTH_SUFFICIENT, UNCONFIRMED_SHORT};
        let mut genes = vec![
            labelled("core", true, LENGTH_SUFFICIENT),
            labelled("short_ok", true, CONSERVED_SHORT),
            labelled("short_no", true, UNCONFIRMED_SHORT),
            labelled("broken", false, LENGTH_SUFFICIENT),
        ];
        let loci = std::collections::HashMap::new();
        assert_eq!(
            super::apply_quality(&mut genes.clone(), Quality::Standard, &loci),
            0
        );
        assert_eq!(super::apply_quality(&mut genes, Quality::High, &loci), 2);
        let ids: Vec<_> = genes.iter().map(|gene| gene.id.as_str()).collect();
        assert_eq!(ids, vec!["core", "short_ok"]);
    }

    #[test]
    fn strict_folds_a_private_allele_and_keeps_a_diverged_paralog() {
        use crate::genetic_code::back_translate;
        use crate::synteny::LENGTH_SUFFICIENT;
        let long_aa = random_protein(80, 7);
        let piece_aa = long_aa[20..50].to_vec();
        let near_aa = long_aa[..76].to_vec();
        let copy_aa = long_aa[..79].to_vec();
        let mut paralog_aa = long_aa.clone();
        for site in [4, 12, 20, 28, 36, 44, 52, 60, 68] {
            paralog_aa[site] = if paralog_aa[site] == b'A' { b'L' } else { b'A' };
        }
        let mut band_aa = long_aa.clone();
        for site in [4, 12, 20, 28, 36, 44, 52] {
            band_aa[site] = if band_aa[site] == b'A' { b'L' } else { b'A' };
        }
        let mut close_aa = long_aa.clone();
        for site in [4, 12, 20] {
            close_aa[site] = if close_aa[site] == b'A' { b'L' } else { b'A' };
        }
        let mut deep_aa = long_aa[20..50].to_vec();
        for site in [0, 8, 16, 24] {
            deep_aa[site] = if deep_aa[site] == b'A' { b'L' } else { b'A' };
        }
        let mut long = labelled("long", true, LENGTH_SUFFICIENT);
        long.reference_dna = back_translate(&long_aa);
        long.presence = vec![true, true, false];
        long.samples_present = 2;
        let mut piece = labelled("piece", true, LENGTH_SUFFICIENT);
        piece.reference_dna = back_translate(&piece_aa);
        piece.presence = vec![true, false, false];
        piece.samples_present = 1;
        let mut alone = labelled("alone", true, LENGTH_SUFFICIENT);
        alone.reference_dna = back_translate(&piece_aa);
        alone.presence = vec![false, false, true];
        alone.samples_present = 1;
        let mut near = labelled("near", true, LENGTH_SUFFICIENT);
        near.reference_dna = back_translate(&near_aa);
        near.presence = vec![true, false, false];
        near.samples_present = 1;
        let mut copy = labelled("copy", true, LENGTH_SUFFICIENT);
        copy.reference_dna = back_translate(&copy_aa);
        copy.presence = vec![true, false, false];
        copy.samples_present = 1;
        let mut twin = labelled("twin", true, LENGTH_SUFFICIENT);
        twin.reference_dna = back_translate(&long_aa);
        twin.presence = vec![true, true, false];
        twin.samples_present = 2;
        let mut paralog = labelled("paralog", true, LENGTH_SUFFICIENT);
        paralog.reference_dna = back_translate(&paralog_aa);
        paralog.presence = vec![true, true, false];
        paralog.samples_present = 2;
        let mut close = labelled("close", true, LENGTH_SUFFICIENT);
        close.reference_dna = back_translate(&close_aa);
        close.presence = vec![true, false, false];
        close.samples_present = 1;
        let mut deep = labelled("deep", true, LENGTH_SUFFICIENT);
        deep.reference_dna = back_translate(&deep_aa);
        deep.presence = vec![true, false, false];
        deep.samples_present = 1;
        let mut band = labelled("band", true, LENGTH_SUFFICIENT);
        band.reference_dna = back_translate(&band_aa);
        band.presence = vec![true, false, false];
        band.samples_present = 1;
        let mut genes = vec![
            long, piece, alone, near, copy, twin, paralog, close, deep, band,
        ];
        let loci = std::collections::HashMap::new();
        assert_eq!(
            super::apply_quality(&mut genes.clone(), Quality::High, &loci),
            0
        );
        assert_eq!(super::apply_quality(&mut genes, Quality::Strict, &loci), 7);
        let ids: Vec<_> = genes.iter().map(|gene| gene.id.as_str()).collect();
        // `alone` is an exact piece of `long` and also of the diverged paralogue,
        // so it is not assigned to either gene.
        assert_eq!(ids, vec!["long", "alone", "paralog"]);
        assert_eq!(genes[0].presence, vec![true, true, false]);
        assert_eq!(genes[0].reference_dna, back_translate(&long_aa));
    }

    #[test]
    fn strict_folds_a_diverged_piece_and_a_repeat_and_keeps_a_full_paralog() {
        use crate::genetic_code::back_translate;
        use crate::synteny::LENGTH_SUFFICIENT;
        let long_aa = random_protein(100, 21);
        let mut piece_aa = long_aa[..95].to_vec();
        for site in [0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84, 90] {
            piece_aa[site] = if piece_aa[site] == b'A' { b'L' } else { b'A' };
        }
        let mut paralog_aa = long_aa.clone();
        for site in [4, 12, 20, 28, 36, 44, 52, 60, 68, 76, 84, 92] {
            paralog_aa[site] = if paralog_aa[site] == b'A' { b'L' } else { b'A' };
        }
        let head = long_aa[..40].to_vec();
        let tail = long_aa[60..].to_vec();
        let mut split_aa = head.clone();
        split_aa.extend_from_slice(&tail);
        let mut inserted = head.clone();
        inserted.extend(std::iter::repeat(b'W').take(20));
        inserted.extend_from_slice(&tail);
        let mut long = labelled("long", true, LENGTH_SUFFICIENT);
        long.reference_dna = back_translate(&long_aa);
        long.presence = vec![true, false, false, false];
        long.samples_present = 1;
        let mut piece = labelled("piece", true, LENGTH_SUFFICIENT);
        piece.reference_dna = back_translate(&piece_aa);
        piece.presence = vec![false, true, false, false];
        piece.samples_present = 1;
        let mut paralog = labelled("paralog", true, LENGTH_SUFFICIENT);
        paralog.reference_dna = back_translate(&paralog_aa);
        paralog.presence = vec![true, true, false, false];
        paralog.samples_present = 2;
        let mut host = labelled("host", true, LENGTH_SUFFICIENT);
        host.reference_dna = back_translate(&inserted);
        host.presence = vec![false, false, true, false];
        host.samples_present = 1;
        let mut split = labelled("split", true, LENGTH_SUFFICIENT);
        split.reference_dna = back_translate(&split_aa);
        split.presence = vec![false, false, false, true];
        split.samples_present = 1;
        let mut genes = vec![long, piece, paralog, host, split];
        let loci = std::collections::HashMap::new();
        assert_eq!(super::apply_quality(&mut genes, Quality::Strict, &loci), 1);
        let ids: Vec<_> = genes.iter().map(|gene| gene.id.as_str()).collect();
        assert_eq!(ids, vec!["long", "paralog", "host", "split"]);
        assert_eq!(genes[0].presence, vec![true, true, false, false]);
        assert_eq!(genes[0].reference_dna, back_translate(&long_aa));
    }

    #[test]
    fn strict_folds_a_repeated_family_and_keeps_a_fragment_of_two_genes() {
        use crate::genetic_code::back_translate;
        use crate::synteny::LENGTH_SUFFICIENT;
        let long_aa = random_protein(80, 13);
        let fragment = long_aa[..30].to_vec();
        let mut long = labelled("long", true, LENGTH_SUFFICIENT);
        long.reference_dna = back_translate(&long_aa);
        let mut presence = vec![false; 71];
        presence[0] = true;
        long.presence = presence;
        long.samples_present = 1;
        let mut genes = vec![long];
        for index in 1..71 {
            let mut copy = labelled(&format!("copy{index}"), true, LENGTH_SUFFICIENT);
            copy.reference_dna = back_translate(&fragment);
            let mut where_present = vec![false; 71];
            where_present[index] = true;
            copy.presence = where_present;
            copy.samples_present = 1;
            genes.push(copy);
        }
        let loci = std::collections::HashMap::new();
        assert_eq!(super::apply_quality(&mut genes, Quality::Strict, &loci), 70);
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].id, "long");
        assert_eq!(genes[0].samples_present, 71);
        assert_eq!(genes[0].reference_dna, back_translate(&long_aa));

        let mut other = long_aa.clone();
        for site in 40..80 {
            other[site] = if other[site] == b'A' { b'G' } else { b'A' };
        }
        let mut left = labelled("left", true, LENGTH_SUFFICIENT);
        left.reference_dna = back_translate(&long_aa);
        left.presence = vec![true, false, false];
        left.samples_present = 1;
        let mut right = labelled("right", true, LENGTH_SUFFICIENT);
        right.reference_dna = back_translate(&other);
        right.presence = vec![false, true, false];
        right.samples_present = 1;
        let mut shared = labelled("shared", true, LENGTH_SUFFICIENT);
        shared.reference_dna = back_translate(&fragment);
        shared.presence = vec![false, false, true];
        shared.samples_present = 1;
        let mut ambiguous = vec![left, right, shared];
        assert_eq!(
            super::apply_quality(&mut ambiguous, Quality::Strict, &loci),
            0
        );
        let ids: Vec<_> = ambiguous.iter().map(|gene| gene.id.as_str()).collect();
        assert_eq!(ids, vec!["left", "right", "shared"]);

        let mut same = Vec::new();
        for index in 0..5 {
            let mut gene = labelled(&format!("same{index}"), true, LENGTH_SUFFICIENT);
            gene.reference_dna = back_translate(&long_aa);
            let mut where_present = vec![false; 5];
            where_present[index] = true;
            gene.presence = where_present;
            gene.samples_present = 1;
            same.push(gene);
        }
        assert_eq!(super::apply_quality(&mut same, Quality::Strict, &loci), 4);
        assert_eq!(same[0].id, "same0");
        assert_eq!(same[0].samples_present, 5);
        assert_eq!(same[0].reference_dna, back_translate(&long_aa));
    }

    #[test]
    fn strict_folds_the_shorter_allele_and_keeps_the_longer_sequence() {
        use crate::genetic_code::back_translate;
        use crate::synteny::LENGTH_SUFFICIENT;
        let long_aa = random_protein(80, 9);
        let short_aa = long_aa[..76].to_vec();
        let mut long = labelled("long", true, LENGTH_SUFFICIENT);
        long.reference_dna = back_translate(&long_aa);
        long.presence = vec![true, false, false];
        long.samples_present = 1;
        long.samples_with_exact_reference = 1;
        let mut short = labelled("short", true, LENGTH_SUFFICIENT);
        short.reference_dna = back_translate(&short_aa);
        short.presence = vec![false, true, true];
        short.samples_present = 2;
        short.samples_with_exact_reference = 2;
        let mut genes = vec![long, short];
        let loci = std::collections::HashMap::new();
        assert_eq!(super::apply_quality(&mut genes, Quality::Strict, &loci), 1);
        assert_eq!(genes[0].id, "long");
        assert_eq!(genes[0].reference_dna, back_translate(&long_aa));
        assert_eq!(genes[0].presence, vec![true, true, true]);
        assert_eq!(genes[0].selection, "intact_longest");
    }

    #[test]
    fn strict_drops_a_nested_piece_a_neighboring_copy_and_folds_an_exclusive_allele() {
        use crate::genetic_code::back_translate;
        use crate::synteny::LENGTH_SUFFICIENT;
        use std::collections::HashMap;
        let long_aa = random_protein(80, 11);
        let mut nested_aa = long_aa[10..40].to_vec();
        for site in [0, 1, 28, 29] {
            nested_aa[site] = if nested_aa[site] == b'A' { b'L' } else { b'A' };
        }
        let mut weak_aa = long_aa[10..40].to_vec();
        for site in [0, 3, 6, 9, 12, 15, 18, 21] {
            weak_aa[site] = if weak_aa[site] == b'A' { b'L' } else { b'A' };
        }
        let beside_aa = long_aa[10..40].to_vec();
        let mut long = labelled("long", true, LENGTH_SUFFICIENT);
        long.reference_dna = back_translate(&long_aa);
        long.presence = vec![true, false];
        long.samples_present = 1;
        let mut nested = labelled("nested", true, LENGTH_SUFFICIENT);
        nested.reference_dna = back_translate(&nested_aa);
        nested.presence = vec![true, false];
        nested.samples_present = 1;
        let mut weak = labelled("weak", true, LENGTH_SUFFICIENT);
        weak.reference_dna = back_translate(&weak_aa);
        weak.presence = vec![true, false];
        weak.samples_present = 1;
        let mut beside = labelled("beside", true, LENGTH_SUFFICIENT);
        beside.reference_dna = back_translate(&beside_aa);
        beside.presence = vec![true, false];
        beside.samples_present = 1;
        let mut loci: HashMap<String, Vec<super::Locus>> = HashMap::new();
        loci.insert(
            "long".to_string(),
            vec![
                super::Locus {
                    sample: 0,
                    contig: 0,
                    start: 100,
                    end: 900,
                    strand: 1,
                },
                super::Locus {
                    sample: 1,
                    contig: 0,
                    start: 100,
                    end: 900,
                    strand: 1,
                },
            ],
        );
        loci.insert(
            "nested".to_string(),
            vec![super::Locus {
                sample: 0,
                contig: 0,
                start: 200,
                end: 400,
                strand: 1,
            }],
        );
        loci.insert(
            "weak".to_string(),
            vec![super::Locus {
                sample: 0,
                contig: 0,
                start: 200,
                end: 400,
                strand: 1,
            }],
        );
        loci.insert(
            "beside".to_string(),
            vec![super::Locus {
                sample: 0,
                contig: 0,
                start: 1_000,
                end: 1_200,
                strand: 1,
            }],
        );
        let split_aa = long_aa.clone();
        let mut diverged_aa = long_aa.clone();
        for site in [4, 12, 20, 28, 36, 44, 52, 60, 68, 76] {
            diverged_aa[site] = if diverged_aa[site] == b'A' {
                b'L'
            } else {
                b'A'
            };
        }
        let mut split = labelled("split", true, LENGTH_SUFFICIENT);
        split.reference_dna = back_translate(&split_aa);
        split.presence = vec![false, true];
        split.lengths = vec![0, 240];
        split.samples_present = 1;
        long.lengths = vec![240, 0];
        let mut diverged = labelled("diverged", true, LENGTH_SUFFICIENT);
        diverged.reference_dna = back_translate(&diverged_aa);
        diverged.presence = vec![false, true];
        diverged.samples_present = 1;
        let mut genes = vec![long, nested, weak, beside, split, diverged];
        assert_eq!(super::apply_quality(&mut genes, Quality::Strict, &loci), 3);
        let ids: Vec<_> = genes.iter().map(|gene| gene.id.as_str()).collect();
        assert_eq!(ids, vec!["long", "weak", "diverged"]);
        assert_eq!(genes[0].presence, vec![true, true]);
        assert_eq!(genes[0].samples_present, 2);
    }

    #[test]
    fn a_one_off_neighbor_does_not_create_a_third_paralog() {
        let dup = random_protein(40, 4);
        let left_a = random_protein(36, 101);
        let right_a = random_protein(36, 202);
        let left_b = random_protein(36, 303);
        let right_b = random_protein(36, 404);
        let junk = random_protein(36, 505);
        let mut samples = Vec::new();
        for index in 0..6 {
            samples.push(sample(
                &format!("s{index}"),
                vec![
                    coding(&left_a, 0, 1, "leftA"),
                    coding(&dup, 0, 300, "dup"),
                    coding(&right_a, 0, 600, "rightA"),
                    coding(&left_b, 0, 900, "leftB"),
                    coding(&dup, 0, 1200, "dup"),
                    coding(&right_b, 0, 1500, "rightB"),
                ],
            ));
        }
        samples.push(sample(
            "odd",
            vec![
                coding(&left_a, 0, 1, "leftA"),
                coding(&dup, 0, 300, "dup"),
                coding(&right_a, 0, 600, "rightA"),
                coding(&left_b, 0, 900, "leftB"),
                coding(&dup, 0, 1200, "dup"),
                coding(&right_b, 0, 1500, "rightB"),
                coding(&left_b, 0, 1800, "leftB"),
                coding(&dup, 0, 2100, "dup"),
                coding(&junk, 0, 2400, "junk"),
            ],
        ));
        let built = build(samples, &params()).unwrap();
        let copies: Vec<_> = built
            .genes
            .iter()
            .filter(|gene| gene.reference_dna == back_translate(&dup))
            .collect();
        assert_eq!(copies.len(), 2, "the private neighbor is not a third copy");
        assert!(copies.iter().all(|gene| gene.samples_present == 7));
        assert!(built.genes.iter().all(|gene| gene.id != "junk"));
    }

    #[test]
    fn a_family_of_only_internal_stops_is_removed_and_one_intact_copy_keeps_it() {
        let left = diverse(36, 4);
        let right = diverse(36, 9);
        let body = diverse(40, 2);
        let mut broken = back_translate(&body);
        broken.replace_range(6..9, "TGA");
        let broken_protein = translate(&broken);
        let pseudo = |name: &str, dna: &str, protein: &[u8]| {
            let mut gene = coding(protein, 0, 400, "pseudo");
            gene.dna = dna.to_string();
            gene.protein = protein.to_vec();
            gene.end = gene.start + gene.dna.len() as u32 - 1;
            let mut sample = sample(
                name,
                vec![
                    coding(&left, 0, 1, "flkL"),
                    gene,
                    coding(&right, 0, 900, "flkR"),
                ],
            );
            sample.contigs[0].sequence = "N".repeat(1000);
            sample
        };
        let only_broken = build(
            (0..4)
                .map(|index| pseudo(&format!("b{index}"), &broken, &broken_protein))
                .collect(),
            &params(),
        )
        .unwrap();
        assert!(only_broken.genes.iter().all(|gene| gene.id != "pseudo"));
        assert_eq!(only_broken.genes.len(), 2);

        let mut mixed: Vec<_> = (0..3)
            .map(|index| pseudo(&format!("b{index}"), &broken, &broken_protein))
            .collect();
        mixed.push(shared_triplet("intact", &left, &body, &right));
        let kept = build(mixed, &params()).unwrap();
        let gene = kept.genes.iter().find(|gene| gene.id == "pseudo").unwrap();
        assert!(gene.intact);
        assert_eq!(gene.reference_dna, back_translate(&body));
        assert_eq!(gene.samples_present, 4);
    }

    #[test]
    fn one_sample_collection_is_not_eaten_from_the_contig_ends() {
        let built = build(
            vec![sample(
                "only",
                vec![
                    coding(&diverse(36, 1), 0, 1, "left"),
                    coding(&diverse(40, 2), 0, 400, "middle"),
                    coding(&diverse(36, 3), 0, 800, "right"),
                ],
            )],
            &params(),
        )
        .unwrap();
        assert_eq!(built.genes.len(), 3);
        let evidence = |id: &str| {
            built
                .genes
                .iter()
                .find(|gene| gene.id == id)
                .unwrap()
                .length_evidence
        };
        assert_eq!(evidence("middle"), crate::synteny::CONSERVED_SHORT);
        assert_eq!(evidence("left"), crate::synteny::UNCONFIRMED_SHORT);
        assert_eq!(evidence("right"), crate::synteny::UNCONFIRMED_SHORT);
    }

    #[test]
    fn shared_stop_is_the_reference_and_the_private_extension_stays_present() {
        let full = diverse(70, 1);
        let short = full[..30].to_vec();
        let flank_l = diverse(40, 8);
        let flank_r = diverse(40, 15);
        let mut samples = Vec::new();
        samples.push(sample(
            "intact",
            vec![
                coding(&flank_l, 0, 1, "flkL"),
                coding(&full, 0, 400, "coreA"),
                coding(&flank_r, 0, 900, "flkR"),
            ],
        ));
        for index in 0..4 {
            samples.push(sample(
                &format!("trunc{index}"),
                vec![
                    coding(&flank_l, 0, 1, "flkL"),
                    coding(&short, 0, 400, "coreA"),
                    coding(&flank_r, 0, 900, "flkR"),
                ],
            ));
        }
        let built = build(samples, &params()).unwrap();
        let gene = built
            .genes
            .iter()
            .find(|gene| gene.id == "coreA")
            .expect("coreA");
        assert_eq!(gene.reference_dna, back_translate(&short));
        assert!(gene.intact);
        assert_eq!(gene.selection, "intact_conserved_end");
        assert!(gene.presence.iter().all(|present| *present));
        assert_eq!(gene.samples_with_exact_reference, 4);
        assert_eq!(gene.extension_samples, 1);
        assert_eq!(gene.modal_end_samples, 4);
    }

    #[test]
    fn contig_edge_alleles_keep_the_full_reference_and_every_samples_membership() {
        let full = random_protein(120, 42);
        let upstream = random_protein(36, 101);
        let downstream = random_protein(36, 202);
        let mut inputs = Vec::new();
        let mut expected_lengths = Vec::new();
        for index in 0..10 {
            let reverse = index % 2 == 1;
            let complete = index < 2;
            let protein = if complete {
                full.as_slice()
            } else {
                &full[30..]
            };
            let mut middle = coding(
                protein,
                0,
                if !complete && !reverse { 19 } else { 400 },
                "gene",
            );
            if !complete {
                middle.dna.replace_range(0..3, "ATG");
                middle.protein = translate(&middle.dna);
            }
            middle.strand = if reverse { -1 } else { 1 };
            expected_lengths.push(middle.dna.len() as u32);
            let mut genes = Vec::new();
            if complete {
                genes.push(coding(
                    if reverse { &downstream } else { &upstream },
                    0,
                    100,
                    "flank",
                ));
                genes.push(coding(
                    if reverse { &upstream } else { &downstream },
                    0,
                    1100,
                    "flank",
                ));
            } else if reverse {
                genes.push(coding(&downstream, 0, 100, "flank"));
            } else {
                genes.push(coding(&downstream, 0, middle.end + 80, "flank"));
            }
            genes.push(middle);
            inputs.push(sample(&format!("s{index:02}"), genes));
        }
        let mut params = params();
        params.audit_members = true;
        let built = build(inputs, &params).unwrap();
        assert_eq!(built.genes.len(), 3);
        let gene = built.genes.iter().find(|gene| gene.id == "gene").unwrap();
        assert!(gene.intact);
        assert_eq!(
            gene.reference_dna,
            back_translate(&full),
            "common truncations must not define the reference"
        );
        assert_eq!(gene.presence, vec![true; 10]);
        assert_eq!(gene.lengths, expected_lengths);
        assert_eq!(gene.samples_with_exact_reference, 2);
        let audit = built.membership_audit.as_ref().unwrap();
        let members: Vec<_> = audit
            .members
            .iter()
            .filter(|member| built.genes[member.gene].id == "gene")
            .collect();
        assert_eq!(members.len(), 10);
        assert!(members
            .iter()
            .all(|member| !member.annotation_id.is_empty() && member.origin == "annotated"));
        assert_eq!(
            members
                .iter()
                .map(|member| member.locus)
                .collect::<HashSet<_>>()
                .len(),
            10
        );
    }

    fn pooled_test_call(
        pools: &mut Pools,
        protein: &[u8],
        contig: u32,
        start: u32,
        strand: i8,
    ) -> Call {
        let dna = back_translate(protein);
        Call {
            sample: 0,
            contig,
            start,
            end: start + dna.len() as u32 - 1,
            strand,
            protein: pools.intern_protein(translate(&dna)),
            dna: pools.intern_dna(dna),
            cluster: 0,
            context: 0,
            locus: None,
            source: Some(0),
            refound: false,
            fusion_piece: false,
            symbol: String::new(),
            product: String::new(),
        }
    }

    #[test]
    fn reference_joins_use_one_supported_locus_on_either_strand() {
        let full = random_protein(80, 42);
        for reverse in [false, true] {
            let mut pools = Pools::default();
            let strand = if reverse { -1 } else { 1 };
            let first = if reverse { &full[40..] } else { &full[..40] };
            let second = if reverse { &full[..40] } else { &full[40..] };
            let a = pooled_test_call(&mut pools, first, 0, 100, strand);
            let b = pooled_test_call(&mut pools, second, 0, a.end + 1, strand);
            let calls = vec![a, b];
            let loci = membership_loci(&calls, &pools);
            assert_eq!(loci[0], loci[1]);
            assert_eq!(
                joined_allele(&[0, 1], &calls, &pools, &loci),
                Some(back_translate(&full))
            );
        }
    }

    #[test]
    fn reference_joins_reject_cross_contig_copies_full_tandems_and_overlaps() {
        let full = random_protein(80, 42);
        for other_contig in [false, true] {
            let mut pools = Pools::default();
            let a = pooled_test_call(&mut pools, &full[..40], 0, 100, 1);
            let b = pooled_test_call(
                &mut pools,
                &full[40..],
                u32::from(other_contig),
                a.end + 1,
                1,
            );
            let calls = vec![a, b];
            if other_contig {
                let loci = membership_loci(&calls, &pools);
                assert_ne!(loci[0], loci[1]);
                assert!(joined_allele(&[0, 1], &calls, &pools, &loci).is_none());
            } else {
                assert!(
                    joined_allele(&[0, 1], &calls, &pools, &[0, 1]).is_none(),
                    "abutting coordinates alone cannot license a join"
                );
            }
        }
        for overlap in [false, true] {
            let mut pools = Pools::default();
            let a = pooled_test_call(&mut pools, &full, 0, 100, 1);
            let b = pooled_test_call(
                &mut pools,
                &full,
                0,
                if overlap { 101 } else { a.end + 1 },
                1,
            );
            let calls = vec![a, b];
            let loci = membership_loci(&calls, &pools);
            assert!(joined_allele(&[0, 1], &calls, &pools, &loci).is_none());
        }
    }

    #[test]
    fn unresolved_copies_on_different_contigs_keep_their_covered_lengths() {
        let protein = random_protein(40, 42);
        let first = coding(&protein, 0, 100, "dup");
        let second = coding(&protein, 1, 100, "dup");
        let bp = first.dna.len() as u32;
        let input = ParsedSample {
            name: "one".to_string(),
            contigs: (0..2)
                .map(|index| Contig {
                    name: format!("contig{index}"),
                    sequence: "N".repeat(1000),
                })
                .collect(),
            genes: vec![first, second],
        };
        let mut params = params();
        params.audit_members = true;
        let built = build(vec![input], &params).unwrap();
        assert_eq!(
            built.genes.len(),
            1,
            "indistinguishable copies stay explicitly unresolved"
        );
        assert_eq!(built.genes[0].reference_dna, back_translate(&protein));
        assert_eq!(built.genes[0].lengths, vec![bp * 2]);
        let members = &built.membership_audit.as_ref().unwrap().members;
        assert_eq!(members.len(), 2);
        assert_ne!(members[0].locus, members[1].locus);
    }

    #[test]
    fn moderate_terminal_alleles_preserve_full_reference_and_sample_lengths() {
        let full = diverse(100, 1);
        let short = full[..90].to_vec();
        let flank_l = diverse(36, 8);
        let flank_r = diverse(36, 15);
        let make_sample = |name: &str, protein: &[u8]| {
            sample(
                name,
                vec![
                    coding(&flank_l, 0, 1, "flkL"),
                    coding(protein, 0, 400, "coreA"),
                    coding(&flank_r, 0, 900, "flkR"),
                ],
            )
        };
        let built = build(
            vec![make_sample("intact", &full), make_sample("short", &short)],
            &params(),
        )
        .unwrap();
        assert_eq!(built.genes.len(), 3);
        let gene = built
            .genes
            .iter()
            .find(|gene| gene.id == "coreA")
            .expect("coreA");
        assert_eq!(gene.reference_dna, back_translate(&full));
        assert!(gene.intact);
        assert_eq!(gene.length_evidence, crate::synteny::LENGTH_SUFFICIENT);
        assert_eq!(gene.presence, vec![true, true]);
        assert_eq!(gene.samples_with_exact_reference, 1);
        assert_eq!(
            gene.lengths,
            vec![
                back_translate(&full).len() as u32,
                back_translate(&short).len() as u32,
            ]
        );
    }

    #[test]
    fn fragments_are_the_same_gene_as_the_intact_copy() {
        let full = diverse(80, 3);
        let left = full[..28].to_vec();
        let right = full[28..].to_vec();
        let flank_l = diverse(36, 8);
        let flank_r = diverse(36, 16);
        let intact = sample(
            "intact",
            vec![
                coding(&flank_l, 0, 1, "flkL"),
                coding(&full, 0, 500, "geneA"),
                coding(&flank_r, 0, 1100, "flkR"),
            ],
        );
        let broken_left = coding(&left, 0, 500, "geneA");
        let mut broken_right = coding(&right, 0, broken_left.end + 1, "geneA");
        broken_right.start = broken_left.end + 1;
        broken_right.end = broken_right.start + broken_right.dna.len() as u32 - 1;
        let broken = sample(
            "broken",
            vec![
                coding(&flank_l, 0, 1, "flkL"),
                broken_left,
                broken_right,
                coding(&flank_r, 0, 1100, "flkR"),
            ],
        );
        let built = build(vec![intact, broken], &params()).unwrap();
        let gene = built.genes.iter().find(|gene| gene.id == "geneA").unwrap();
        assert_eq!(gene.reference_dna, back_translate(&full));
        assert_eq!(gene.samples_present, 2);
    }

    #[test]
    fn paralogs_become_two_genes() {
        let paralog = diverse(40, 4);
        let left = diverse(36, 6);
        let mid = diverse(36, 9);
        let right = diverse(36, 12);
        let copy_a = coding(&paralog, 0, 300, "dup");
        let middle = coding(&mid, 0, copy_a.end + 200, "mid");
        let copy_b = coding(&paralog, 0, middle.end + 200, "dup");
        let built = build(
            vec![sample(
                "one",
                vec![
                    coding(&left, 0, 1, "left"),
                    copy_a,
                    middle,
                    copy_b,
                    coding(&right, 0, 2000, "right"),
                ],
            )],
            &params(),
        )
        .unwrap();
        let copies: Vec<_> = built
            .genes
            .iter()
            .filter(|gene| gene.reference_dna == back_translate(&paralog))
            .collect();
        assert_eq!(copies.len(), 2, "paralogs stayed in one group");
        assert!(copies.iter().all(|gene| gene.samples_present == 1));
    }

    #[test]
    fn touching_full_length_paralogs_are_not_joined_as_fragments() {
        let protein = diverse(40, 4);
        let first = coding(&protein, 0, 300, "dup");
        let second = coding(&protein, 0, first.end + 1, "dup");
        let built = build(
            vec![sample(
                "one",
                vec![
                    coding(&diverse(36, 6), 0, 1, "left"),
                    first,
                    second,
                    coding(&diverse(36, 12), 0, 1000, "right"),
                ],
            )],
            &params(),
        )
        .unwrap();
        let copies: Vec<_> = built
            .genes
            .iter()
            .filter(|gene| gene.reference_dna == back_translate(&protein))
            .collect();
        assert_eq!(copies.len(), 2);
        assert!(copies.iter().all(|gene| gene.samples_present == 1));
        assert!(copies
            .iter()
            .all(|gene| gene.lengths == vec![back_translate(&protein).len() as u32]));
    }

    #[test]
    fn fragmented_paralog_preserves_two_full_references_and_presence() {
        let full = random_protein(80, 42);
        let left_a = random_protein(36, 101);
        let right_a = random_protein(36, 202);
        let left_b = random_protein(36, 303);
        let right_b = random_protein(36, 404);
        let mut samples = Vec::new();
        for index in 0..3 {
            let mut genes = vec![coding(&left_a, 0, 1, "leftA")];
            if index == 1 {
                genes.push(coding(&full, 0, 400, "dupA"));
            } else {
                let head = coding(&full[..40], 0, 400, "dupA");
                let tail = coding(&full[40..], 0, head.end + 1, "dupA");
                genes.extend([head, tail]);
            }
            genes.extend([
                coding(&right_a, 0, 900, "rightA"),
                coding(&left_b, 0, 1100, "leftB"),
                coding(&full, 0, 1400, "dupB"),
                coding(&right_b, 0, 2000, "rightB"),
            ]);
            samples.push(sample(&format!("sample{index}"), genes));
        }
        let built = build(samples, &params()).unwrap();
        assert_eq!(built.genes.len(), 6);
        let copies: Vec<_> = built
            .genes
            .iter()
            .filter(|gene| gene.reference_dna == back_translate(&full))
            .collect();
        assert_eq!(copies.len(), 2);
        assert!(copies
            .iter()
            .all(|gene| gene.intact && gene.presence == vec![true; 3]));
        let lengths: Vec<_> = copies.iter().map(|gene| gene.lengths.clone()).collect();
        let full_bp = back_translate(&full).len() as u32;
        assert!(lengths.contains(&vec![full_bp; 3]));
        assert!(lengths.contains(&vec![full_bp + 3, full_bp, full_bp + 3]));
    }

    #[test]
    fn a_fusion_is_recorded_as_both_genes() {
        let gene_a = diverse(40, 2);
        let gene_b = diverse(36, 14);
        let mut fusion = gene_a.clone();
        fusion.extend_from_slice(&gene_b[1..]);
        let mut samples = Vec::new();
        for index in 0..3 {
            samples.push(sample(
                &format!("a{index}"),
                vec![
                    coding(&diverse(30, 5), 0, 1, "flk"),
                    coding(&gene_a, 0, 400, "aaa"),
                ],
            ));
            samples.push(sample(
                &format!("b{index}"),
                vec![coding(&gene_b, 0, 400, "bbb")],
            ));
        }
        samples.push(sample(
            "fused",
            vec![
                coding(&diverse(30, 5), 0, 1, "flk"),
                coding(&fusion, 0, 400, "aaa"),
            ],
        ));
        let mut audit_params = params();
        audit_params.audit_members = true;
        let built = build(samples, &audit_params).unwrap();
        let aaa = built.genes.iter().find(|gene| gene.id == "aaa").unwrap();
        let bbb = built.genes.iter().find(|gene| gene.id == "bbb").unwrap();
        assert_eq!(aaa.reference_dna, back_translate(&gene_a));
        assert_eq!(bbb.reference_dna, back_translate(&gene_b));
        let fused = built
            .sample_names
            .iter()
            .position(|name| name == "fused")
            .unwrap();
        assert!(aaa.presence[fused]);
        assert!(bbb.presence[fused]);
        assert!(built
            .genes
            .iter()
            .all(|gene| gene.reference_dna != back_translate(&fusion)));
        let fusion_members: Vec<_> = built
            .membership_audit
            .as_ref()
            .unwrap()
            .members
            .iter()
            .filter(|member| member.sample as usize == fused && member.origin == "fusion_piece")
            .collect();
        assert_eq!(fusion_members.len(), 2);
        assert!(fusion_members
            .iter()
            .all(|member| member.annotation_id == "cds_0_400"));
        assert_ne!(fusion_members[0].gene, fusion_members[1].gene);
    }

    #[test]
    fn a_missed_gene_is_refound_between_the_same_flanks() {
        let flank_l = diverse(36, 7);
        let gene = diverse(40, 3);
        let flank_r = diverse(36, 18);
        let gene_dna = back_translate(&gene);
        let annotated = |name: &str| {
            let left = coding(&flank_l, 0, 1, "flkL");
            let middle = coding(&gene, 0, left.end + 80, "missed");
            let right = coding(&flank_r, 0, middle.end + 80, "flkR");
            let mut sequence = vec![b'N'; right.end as usize + 10];
            paint(&mut sequence, &left);
            paint(&mut sequence, &middle);
            paint(&mut sequence, &right);
            sample_with(name, sequence, vec![left, middle, right])
        };
        let left = coding(&flank_l, 0, 1, "flkL");
        // The CDS sits inside the gap, unannotated. The flanks are the same
        // pair that brackets it in the other samples, which is the only
        // evidence refind is allowed to use.
        let gene_start = left.end + 41;
        let right = coding(&flank_r, 0, gene_start + gene_dna.len() as u32 + 48, "flkR");
        let mut sequence = vec![b'N'; right.end as usize + 10];
        paint(&mut sequence, &left);
        let at = (gene_start - 1) as usize;
        sequence[at..at + gene_dna.len()].copy_from_slice(gene_dna.as_bytes());
        paint(&mut sequence, &right);
        let gap = right.start as i32 - left.end as i32 - 1;
        assert!((60..=12_000).contains(&gap), "{gap}");

        let mut params = params();
        params.refind = true;
        params.audit_members = true;
        let built = build(
            vec![
                annotated("keepA"),
                annotated("keepB"),
                sample_with("gap", sequence, vec![left, right]),
            ],
            &params,
        )
        .unwrap();
        let row = built
            .genes
            .iter()
            .find(|gene| gene.id == "missed")
            .expect("missed");
        let gap_sample = built
            .sample_names
            .iter()
            .position(|name| name == "gap")
            .unwrap();
        assert!(row.presence[gap_sample]);
        assert_eq!(row.lengths[gap_sample], gene_dna.len() as u32);
        assert_eq!(row.reference_dna, gene_dna);
        assert_eq!(row.samples_present, 3);
        let found: Vec<_> = built
            .membership_audit
            .as_ref()
            .unwrap()
            .members
            .iter()
            .filter(|member| member.origin == "refound")
            .collect();
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].sample as usize, gap_sample);
        assert!(found[0].annotation_id.is_empty());
        assert_eq!(built.genes[found[0].gene].id, "missed");
    }

    fn paint(sequence: &mut [u8], gene: &ParsedGene) {
        let start = (gene.start - 1) as usize;
        let dna = gene.dna.as_bytes();
        sequence[start..start + dna.len()].copy_from_slice(dna);
    }

    fn sample_with(name: &str, sequence: Vec<u8>, genes: Vec<ParsedGene>) -> ParsedSample {
        ParsedSample {
            name: name.to_string(),
            contigs: vec![Contig {
                name: "contig".to_string(),
                sequence: String::from_utf8(sequence).expect("contig is ascii"),
            }],
            genes,
        }
    }
}
