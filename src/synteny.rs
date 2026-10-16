//! Ortholog groups are homology plus gene order.
//!
//! Sequence identity alone keeps a fragment out of the full gene and keeps a
//! fused CDS glued to whichever piece is longer. Gene order repairs both.
//! Paralogs that landed in one cluster are split by their flanking genes.
//! Fusions are given back to the two genes they contain. Fragments are adopted
//! by the full gene they are a piece of.

use std::collections::{HashMap, HashSet};

use rayon::prelude::*;

use crate::align::{fit_short_in_long, global_identity, terminal_identity_parts};
use crate::genetic_code::{is_intact, is_start, is_stop, orf_status, translate, translate_bytes};
use crate::log::{self, progress};
use crate::model::{Call, Pools};

const CONTAIN_IDENTITY: f64 = 0.90;
const CONTAIN_COVERAGE: f64 = 0.85;
const MAX_COOCCUR: f64 = 0.25;
const NEIGHBOR_GAP: i32 = 15_000;
const ABUT_GAP: i32 = 3;
const TERMINAL_IDENTITY: f64 = 0.98;
pub(crate) const ALLELE_EDGE_BP: usize = 90;
const ALLELE_MIN_AA: usize = 30;

/// Split fusion CDS into the two genes they contain. Returns the number split.
pub fn split_fusions(calls: &mut Vec<Call>, pools: &mut Pools) -> usize {
    let (reps, counts) = rayon::join(
        || representatives(calls, pools),
        || sample_counts(calls),
    );
    let kmers = kmer_sets(&reps);
    let mut cluster_ids: Vec<u32> = reps.keys().copied().collect();
    // Sorted so the `best` tie-breaks below (strict `>` on a support score)
    // resolve the same way on every run, instead of depending on this
    // HashMap's randomized iteration order.
    cluster_ids.sort_unstable();
    // Finding the best (left, right) split for one `long_id` only reads
    // `reps`/`counts`/`kmers` - nothing here is mutated until the groups are
    // actually applied below - so every family's search is independent of
    // every other family's and can run on its own thread.
    let planned: Vec<(u32, u32, u32, usize)> = cluster_ids
        .par_iter()
        .filter_map(|&long_id| {
            best_fusion_split(long_id, &reps, &counts, &kmers, &cluster_ids)
                .map(|(left_id, right_id, end)| (long_id, left_id, right_id, end))
        })
        .collect();
    if planned.is_empty() {
        return 0;
    }
    let mut twins = Vec::new();
    let mut split = 0usize;
    for call in calls.iter_mut() {
        let Some((_, left_id, right_id, end)) = planned.iter().find(|item| item.0 == call.cluster)
        else {
            continue;
        };
        let dna = pools.dna[call.dna as usize].clone();
        call.fusion_piece = true;
        let protein_len = pools.protein[call.protein as usize].len();
        let cut = (*end).min(protein_len);
        if let Some(suffix) = slice_coding(&dna, cut, protein_len) {
            let suffix_bp = suffix.len() as u32;
            let mut twin = call.clone();
            twin.cluster = *right_id;
            twin.symbol.clear();
            twin.product.clear();
            twin.protein = pools.intern_protein(translate(&suffix));
            twin.dna = pools.intern_dna(suffix);
            if twin.strand >= 0 {
                twin.start = twin.end.saturating_sub(suffix_bp).max(twin.start);
            } else {
                twin.end = twin.start.saturating_add(suffix_bp).min(twin.end);
            }
            twins.push(twin);
        }
        if let Some(prefix) = slice_coding(&dna, 0, cut) {
            let prefix_bp = prefix.len() as u32;
            call.protein = pools.intern_protein(translate(&prefix));
            call.dna = pools.intern_dna(prefix);
            if call.strand >= 0 {
                call.end = call.start.saturating_add(prefix_bp).saturating_sub(1);
            } else {
                call.start = call.end.saturating_sub(prefix_bp).saturating_add(1);
            }
        }
        call.cluster = *left_id;
        split += 1;
    }
    calls.extend(twins);
    split
}

/// The best `(left_id, right_id, cut_point)` fusion split for `long_id`, if
/// any: `left` has to explain the start of `long`, and the leftover suffix
/// has to be explained by `right`. Pulled out of `split_fusions` so the
/// search for each `long_id` can run as an independent unit of parallel work.
fn best_fusion_split(
    long_id: u32,
    reps: &HashMap<u32, Vec<u8>>,
    counts: &HashMap<u32, u32>,
    kmers: &HashMap<u32, HashSet<&[u8]>>,
    cluster_ids: &[u32],
) -> Option<(u32, u32, usize)> {
    let long = reps.get(&long_id)?;
    let long_samples = counts.get(&long_id).copied().unwrap_or(0);
    let mut best: Option<(u32, u32, usize, u32)> = None;
    for &left_id in cluster_ids {
        if left_id == long_id || counts.get(&left_id).copied().unwrap_or(0) <= long_samples {
            continue;
        }
        let Some(left) = reps.get(&left_id) else {
            continue;
        };
        if !shares_kmer_cached(left_id, long_id, left.len(), long.len(), kmers) {
            continue;
        }
        let Some(fit) = fit_short_in_long(left, long) else {
            continue;
        };
        if fit.identity < CONTAIN_IDENTITY
            || fit.coverage < CONTAIN_COVERAGE
            || fit.long_start > 2
            || fit.internal_gap > 5
        {
            continue;
        }
        let suffix = &long[fit.long_end.min(long.len())..];
        for &right_id in cluster_ids {
            if right_id == long_id || right_id == left_id {
                continue;
            }
            if counts.get(&right_id).copied().unwrap_or(0) <= long_samples {
                continue;
            }
            let Some(right) = reps.get(&right_id) else {
                continue;
            };
            let matched = if right.len() <= suffix.len() {
                fit_short_in_long(right, suffix)
            } else if suffix.len() >= 4 {
                fit_short_in_long(suffix, right)
            } else {
                None
            };
            if !matched.is_some_and(|fit| {
                fit.identity >= 0.85 && fit.coverage >= 0.80 && fit.internal_gap <= 5
            }) {
                continue;
            }
            let support = counts.get(&left_id).copied().unwrap_or(0)
                + counts.get(&right_id).copied().unwrap_or(0);
            if best.map(|(_, _, _, score)| support > score).unwrap_or(true) {
                best = Some((left_id, right_id, fit.long_end, support));
            }
        }
    }
    best.map(|(left_id, right_id, end, _)| (left_id, right_id, end))
}

/// Capture homology anchors and supported loci after fusion repair, before splitting.
pub fn freeze_contexts(calls: &mut [Call], pools: &Pools) {
    sort_calls(calls);
    for call in calls.iter_mut() {
        call.context = call.cluster;
        call.locus = None;
    }
    retain_loci(calls, pools);
}

pub fn split_paralogs(calls: &mut [Call], pools: &Pools) -> usize {
    sort_calls(calls);
    let loci = physical_loci(calls, pools);
    let contexts = locus_contexts(calls, &loci, |call| call.context);
    // A neighbor seen once must not split a copy away from the flank the rest
    // of the collection supports.
    let contexts = stabilize_contexts(calls, &loci, &contexts, |call| call.context);
    let mut by_cluster: HashMap<u32, Vec<usize>> = HashMap::new();
    for (index, locus) in loci.iter().enumerate() {
        by_cluster.entry(calls[locus.first].cluster).or_default().push(index);
    }
    let mut next_id = calls.iter().map(|call| call.cluster).max().unwrap_or(0) + 1;
    let mut cluster_list: Vec<u32> = by_cluster.keys().copied().collect();
    cluster_list.sort_unstable();
    let plans: Vec<ParalogPlan> = cluster_list
        .par_iter()
        .map(|&cluster| paralog_split_groups(&by_cluster[&cluster], calls, &loci, &contexts))
        .collect();
    let mut added = 0;
    let mut unresolved = 0;
    let mut allele_unions = 0;
    for plan in plans {
        unresolved += plan.unresolved_loci;
        allele_unions += plan.allele_unions;
        for group in plan.groups {
            let id = next_id;
            next_id += 1;
            for &index in &group {
                calls[index].cluster = id;
            }
            added += 1;
        }
    }
    if unresolved > 0 {
        progress!(
            "rustroo: {unresolved} loci have unresolved paralog placement; no arbitrary flank assignment"
        );
    }
    if allele_unions > 0 {
        progress!("rustroo: reconciled {allele_unions} exact-allele context groups");
    }
    added
}

#[derive(Default)]
struct ParalogPlan {
    groups: Vec<Vec<usize>>,
    unresolved_loci: usize,
    allele_unions: usize,
}

struct ContextBucket {
    context: (Option<u32>, Option<u32>),
    loci: Vec<usize>,
    samples: HashSet<u32>,
    anchored: bool,
}

fn paralog_split_groups(
    indexes: &[usize],
    calls: &[Call],
    loci: &[Locus],
    contexts: &[(Option<u32>, Option<u32>)],
) -> ParalogPlan {
    let mut per_sample: HashMap<u32, usize> = HashMap::new();
    for &index in indexes {
        *per_sample.entry(calls[loci[index].first].sample).or_default() += 1;
    }
    if per_sample.values().all(|&count| count <= 1) {
        return ParalogPlan::default();
    }
    let mut by_context: HashMap<(Option<u32>, Option<u32>), Vec<usize>> = HashMap::new();
    for &index in indexes {
        by_context.entry(contexts[index]).or_default().push(index);
    }
    let mut ordered: Vec<_> = by_context.into_iter().collect();
    ordered.sort_by_key(|&(context, _)| context);
    let mut unresolved_loci = 0;
    let buckets: Vec<ContextBucket> = ordered
        .into_iter()
        .map(|(context, indexes)| {
            let mut samples = HashSet::new();
            for &index in &indexes {
                if !samples.insert(calls[loci[index].first].sample) {
                    // Equal contexts cannot resolve two copies in one sample.
                    unresolved_loci += 1;
                }
            }
            let anchored = samples.iter().any(|sample| per_sample[sample] > 1);
            ContextBucket { context, loci: indexes, samples, anchored }
        })
        .collect();
    let mut left: HashMap<u32, Vec<usize>> = HashMap::new();
    let mut right: HashMap<u32, Vec<usize>> = HashMap::new();
    for (id, bucket) in buckets.iter().enumerate() {
        if let Some(flank) = bucket.context.0 {
            left.entry(flank).or_default().push(id);
        }
        if let Some(flank) = bucket.context.1 {
            right.entry(flank).or_default().push(id);
        }
    }
    let mut pairs = Vec::new();
    for (id, bucket) in buckets.iter().enumerate() {
        let mut candidates = Vec::new();
        if let Some(group) = bucket.context.0.and_then(|flank| left.get(&flank)) {
            candidates.extend_from_slice(group);
        }
        if let Some(group) = bucket.context.1.and_then(|flank| right.get(&flank)) {
            candidates.extend_from_slice(group);
        }
        candidates.sort_unstable();
        candidates.dedup();
        for other in candidates {
            if other <= id || !bucket.samples.is_disjoint(&buckets[other].samples) {
                continue;
            }
            pairs.push((id, other));
        }
    }
    let mut parents: HashMap<u32, u32> = HashMap::new();
    let mut samples: Vec<HashSet<u32>> = buckets.iter().map(|bucket| bucket.samples.clone()).collect();
    reconcile_disjoint_pairs(pairs, &mut parents, &mut samples);
    // Exact whole alleles can bridge rearranged contexts. Use roots after
    // flank reconciliation so weaker sequence evidence cannot undo it.
    let mut by_protein: HashMap<u32, Vec<usize>> = HashMap::new();
    for (id, bucket) in buckets.iter().enumerate() {
        let root = find(&mut parents, id as u32) as usize;
        for &index in &bucket.loci {
            let locus = &loci[index];
            if locus.first == locus.last {
                by_protein.entry(calls[locus.first].protein).or_default().push(root);
            }
        }
    }
    let mut pairs = Vec::new();
    for roots in by_protein.values_mut() {
        roots.sort_unstable();
        roots.dedup();
        for (index, &a) in roots.iter().enumerate() {
            for &b in &roots[index + 1..] {
                if samples[a].is_disjoint(&samples[b]) {
                    pairs.push((a, b));
                }
            }
        }
    }
    pairs.sort_unstable();
    pairs.dedup();
    let allele_unions = reconcile_disjoint_pairs(pairs, &mut parents, &mut samples);
    let mut grouped: HashMap<u32, (Vec<usize>, bool)> = HashMap::new();
    for (id, bucket) in buckets.into_iter().enumerate() {
        let root = find(&mut parents, id as u32);
        let group = grouped.entry(root).or_default();
        group.0.extend(bucket.loci);
        group.1 |= bucket.anchored;
    }
    let mut groups = Vec::new();
    let mut unplaced = Vec::new();
    for (indexes, anchored) in grouped.into_values() {
        if anchored {
            groups.push(indexes);
        } else {
            unplaced.extend(indexes);
        }
    }
    // Unplaced single-copy samples have no evidence of distinct copies from
    // each other. Retain one unresolved family, not one per new context.
    unresolved_loci += unplaced.len();
    if !unplaced.is_empty() {
        groups.push(unplaced);
    }
    let mut groups: Vec<Vec<usize>> = groups
        .into_iter()
        .map(|group| {
            let mut indexes = Vec::new();
            for index in group {
                indexes.extend(loci[index].first..=loci[index].last);
            }
            indexes.sort_unstable();
            indexes
        })
        .collect();
    groups.sort_by(|a, b| b.len().cmp(&a.len()).then(a.cmp(b)));
    ParalogPlan { groups: groups.into_iter().skip(1).collect(), unresolved_loci, allele_unions }
}

fn reconcile_disjoint_pairs(
    mut pairs: Vec<(usize, usize)>,
    parents: &mut HashMap<u32, u32>,
    samples: &mut [HashSet<u32>],
) -> usize {
    if pairs.is_empty() {
        return 0;
    }
    let mut partners = vec![Vec::new(); samples.len()];
    for &(a, b) in &pairs {
        partners[a].push(b);
        partners[b].push(a);
    }
    // Co-occurring eligible hosts make a context or allele assignment ambiguous.
    let ambiguous: Vec<bool> = partners
        .iter()
        .map(|hosts| {
            let mut seen = HashSet::new();
            hosts.iter().any(|&host| {
                samples[host].iter().any(|&sample| !seen.insert(sample))
            })
        })
        .collect();
    pairs.sort_by(|&(a, b), &(c, d)| {
        let support = |left: usize, right: usize| {
            samples[left].len() + samples[right].len()
        };
        support(c, d).cmp(&support(a, b)).then((a, b).cmp(&(c, d)))
    });
    let mut merged = 0;
    for (a, b) in pairs {
        if ambiguous[a] || ambiguous[b] {
            continue;
        }
        let a = find(parents, a as u32) as usize;
        let b = find(parents, b as u32) as usize;
        if a == b || !samples[a].is_disjoint(&samples[b]) {
            continue;
        }
        let (keep, child) = if a < b { (a, b) } else { (b, a) };
        let added = std::mem::take(&mut samples[child]);
        samples[keep].extend(added);
        parents.insert(child as u32, keep as u32);
        merged += 1;
    }
    merged
}

pub fn merge_fragments(calls: &mut [Call], pools: &Pools, contig_lengths: &[Vec<usize>]) -> bool {
    sort_calls(calls);
    // Representatives, sample sets, loci, and adjacency only read the calls.
    // Building them one after another kept every core but one idle.
    let ((reps, members), (loci, adjacency)) = rayon::join(
        || rayon::join(|| representatives(calls, pools), || sample_sets(calls)),
        || rayon::join(|| physical_loci(calls, pools), || adjacency_pairs(calls)),
    );
    let contexts = locus_contexts(calls, &loci, |call| call.cluster);
    let contexts = stabilize_contexts(calls, &loci, &contexts, |call| call.cluster);
    let ((flanks, kmers), alleles) = rayon::join(
        || rayon::join(
            || modal_flanks_of(calls, &loci, &contexts),
            || kmer_sets(&reps),
        ),
        || AlleleIndex::new(calls, pools, &loci, &contexts, contig_lengths),
    );
    let index = MergeIndex::new(&flanks, &adjacency);
    let mut ids: Vec<u32> = reps.keys().copied().collect();
    ids.sort_unstable();
    let adopted: Vec<(u32, AlleleDecision, Option<MergeTarget>)> = ids
        .par_iter()
        .map(|&short_id| {
            let allele = alleles.target(short_id, pools);
            if allele.target.is_some() || allele.ambiguous {
                return (short_id, allele, None);
            }
            let candidates = index.candidates(short_id, &flanks);
            let legacy = best_merge_target(short_id, &reps, &members, &flanks, &kmers, &candidates);
            (short_id, allele, legacy)
        })
        .collect();
    let mut legacy = Vec::new();
    let mut terminal = Vec::new();
    let mut ambiguous = 0;
    let mut rescues = 0;
    for (child, allele, target) in adopted {
        ambiguous += usize::from(allele.ambiguous);
        if let Some(parent) = allele.target {
            rescues += 1;
            terminal.push((child, parent));
        } else if let Some(target) = target {
            if target.terminal {
                terminal.push((child, target.cluster));
            } else {
                legacy.push((child, target.cluster));
            }
        }
    }
    if ambiguous > 0 {
        progress!("rustroo: {ambiguous} families have ambiguous shortened-allele targets; no arbitrary adoption");
    }
    if rescues > 0 {
        progress!("rustroo: proposed {rescues} strand-aware shortened-allele adoptions; unions remain sample-guarded");
    }
    legacy.extend(abutting_pieces(calls, &reps, &members));
    let (merged, terminal_merged, blocked) = apply_adoptions(calls, pools, members, &legacy, &terminal);
    if merged > 0 {
        progress!(
            "rustroo: merged {merged} groups, including {terminal_merged} guarded shortened-allele groups"
        );
    }
    if blocked > 0 {
        progress!("rustroo: blocked {blocked} merges that would combine distinct co-occurring loci");
    }
    merged > 0
}

#[derive(Default)]
struct AlleleDecision {
    target: Option<u32>,
    ambiguous: bool,
}

struct AlleleOccurrence {
    protein: u32,
    flanks: (Option<u32>, Option<u32>),
    upstream_edge: bool,
    downstream_edge: bool,
}

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum AlleleSide {
    Either,
    MissingPrefix,
    MissingSuffix,
}

struct AlleleIndex {
    occurrences: HashMap<u32, Vec<AlleleOccurrence>>,
    representatives: HashMap<u32, u32>,
    intact_representatives: HashSet<u32>,
    pairs: HashMap<(u32, u32), HashSet<u32>>,
    upstream: HashMap<u32, HashSet<u32>>,
    downstream: HashMap<u32, HashSet<u32>>,
    replicated: HashMap<u32, HashSet<(u32, u32)>>,
}

impl AlleleIndex {
    fn new(
        calls: &[Call], pools: &Pools, loci: &[Locus],
        contexts: &[(Option<u32>, Option<u32>)], contig_lengths: &[Vec<usize>],
    ) -> Self {
        let representatives = representative_ids(calls, pools);
        let mut occurrences: HashMap<u32, Vec<AlleleOccurrence>> = HashMap::new();
        let mut samples: HashMap<u32, HashSet<u32>> = HashMap::new();
        let mut ineligible = HashSet::new();
        let mut observed: HashMap<(u32, u32), HashSet<u32>> = HashMap::new();
        let mut votes: HashMap<(u32, u32), HashMap<u32, usize>> = HashMap::new();
        let mut intact_representatives = HashSet::new();
        let mut intact_cache = HashMap::new();
        for (position, (locus, &context)) in loci.iter().zip(contexts).enumerate() {
            let call = &calls[locus.first];
            let cluster = call.cluster;
            if locus.first != locus.last || !samples.entry(cluster).or_default().insert(call.sample) {
                ineligible.insert(cluster);
            }
            let length = contig_lengths[call.sample as usize][call.contig as usize];
            let first_on_contig = position == 0 || {
                let previous = &calls[loci[position - 1].first];
                previous.sample != call.sample || previous.contig != call.contig
            };
            let last_on_contig = position + 1 == loci.len() || {
                let next = &calls[loci[position + 1].first];
                next.sample != call.sample || next.contig != call.contig
            };
            let start_edge = first_on_contig && context.0.is_none() && call.start as usize <= ALLELE_EDGE_BP;
            let end_edge = last_on_contig && context.1.is_none() && locus.end as usize <= length
                && length - locus.end as usize <= ALLELE_EDGE_BP;
            let (flanks, upstream_edge, downstream_edge) = if call.strand >= 0 {
                (context, start_edge, end_edge)
            } else {
                ((context.1, context.0), end_edge, start_edge)
            };
            occurrences.entry(cluster).or_default().push(AlleleOccurrence {
                protein: call.protein, flanks, upstream_edge, downstream_edge,
            });
            let intact = *intact_cache.entry(call.dna)
                .or_insert_with(|| is_intact(&pools.dna[call.dna as usize]));
            let representative_length = pools.protein[representatives[&cluster] as usize].len();
            let full_length = (locus.first..=locus.last).any(|index| {
                pools.protein[calls[index].protein as usize].len() * 100 >= representative_length * 95
            });
            if full_length {
                if let (Some(upstream), Some(downstream)) = flanks {
                    observed.entry((upstream, downstream)).or_default().insert(cluster);
                    if locus.first == locus.last && intact {
                        *votes.entry((upstream, downstream)).or_default().entry(cluster).or_default() += 1;
                    }
                }
            }
            if call.protein == representatives[&cluster] && intact {
                intact_representatives.insert(cluster);
            }
        }
        occurrences.retain(|cluster, _| !ineligible.contains(cluster));
        let mut pairs: HashMap<(u32, u32), HashSet<u32>> = HashMap::new();
        let mut upstream: HashMap<u32, HashSet<u32>> = HashMap::new();
        let mut downstream: HashMap<u32, HashSet<u32>> = HashMap::new();
        let mut replicated: HashMap<u32, HashSet<(u32, u32)>> = HashMap::new();
        for (flanks, groups) in observed {
            for cluster in groups {
                pairs.entry(flanks).or_default().insert(cluster);
                upstream.entry(flanks.0).or_default().insert(cluster);
                downstream.entry(flanks.1).or_default().insert(cluster);
            }
        }
        for (flanks, groups) in votes {
            for (cluster, support) in groups {
                if support >= 2 {
                    replicated.entry(cluster).or_default().insert(flanks);
                }
            }
        }
        Self { occurrences, representatives, intact_representatives, pairs, upstream, downstream, replicated }
    }

    fn hosts(&self, occurrence: &AlleleOccurrence) -> Option<&HashSet<u32>> {
        match occurrence.flanks {
            (Some(upstream), Some(downstream)) => self.pairs.get(&(upstream, downstream)),
            (None, Some(downstream)) if occurrence.upstream_edge => self.downstream.get(&downstream),
            (Some(upstream), None) if occurrence.downstream_edge => self.upstream.get(&upstream),
            _ => None,
        }
    }

    fn target(&self, short_id: u32, pools: &Pools) -> AlleleDecision {
        let Some(occurrences) = self.occurrences.get(&short_id) else {
            return AlleleDecision::default();
        };
        let mut hosts: Option<HashSet<u32>> = None;
        for occurrence in occurrences {
            let Some(eligible) = self.hosts(occurrence) else {
                return AlleleDecision::default();
            };
            if let Some(hosts) = &mut hosts {
                hosts.retain(|host| eligible.contains(host));
            } else {
                hosts = Some(eligible.clone());
            }
        }
        let mut candidates: Vec<_> = hosts.unwrap_or_default().into_iter().collect();
        candidates.sort_unstable();
        let mut matched = None;
        let mut cache: HashMap<(u32, u32, AlleleSide), bool> = HashMap::new();
        for candidate in candidates {
            if candidate == short_id {
                continue;
            }
            let long_id = self.representatives[&candidate];
            let long = &pools.protein[long_id as usize];
            let compatible = occurrences.iter().all(|occurrence| {
                let side = match occurrence.flanks {
                    (None, Some(_)) => AlleleSide::MissingPrefix,
                    (Some(_), None) => AlleleSide::MissingSuffix,
                    _ => AlleleSide::Either,
                };
                *cache.entry((occurrence.protein, long_id, side)).or_insert_with(|| {
                    shortened_identity(&pools.protein[occurrence.protein as usize], long, side)
                        .is_some_and(|identity| identity + f64::EPSILON >= TERMINAL_IDENTITY)
                })
            });
            if !compatible {
                continue;
            }
            // Count all sequence/context hosts before sample-overlap rejection.
            // A co-occurring host must not make a competing host look unique.
            if matched.is_some() {
                return AlleleDecision { target: None, ambiguous: true };
            }
            matched = Some(candidate);
        }
        let Some(candidate) = matched else {
            return AlleleDecision::default();
        };
        if !self.occurrences.contains_key(&candidate) || !self.intact_representatives.contains(&candidate) {
            return AlleleDecision { target: None, ambiguous: true };
        }
        let supported = self.replicated.get(&candidate).is_some_and(|witnesses| {
            occurrences.iter().all(|occurrence| match occurrence.flanks {
                (Some(upstream), Some(downstream)) => witnesses.contains(&(upstream, downstream)),
                (None, Some(downstream)) => witnesses.iter().any(|&(_, flank)| flank == downstream),
                (Some(upstream), None) => witnesses.iter().any(|&(flank, _)| flank == upstream),
                _ => false,
            })
        });
        AlleleDecision { target: supported.then_some(candidate), ambiguous: false }
    }
}

fn shortened_identity(short: &[u8], long: &[u8], side: AlleleSide) -> Option<f64> {
    if short.len() < ALLELE_MIN_AA || short.len() >= long.len() || short.len() * 2 < long.len()
        || short.iter().chain(long).any(|&residue| residue == b'X' || residue == b'*')
    {
        return None;
    }
    let window = match side {
        AlleleSide::MissingPrefix => &long[long.len().saturating_sub(short.len() + 2)..],
        AlleleSide::MissingSuffix => &long[..(short.len() + 2).min(long.len())],
        AlleleSide::Either => long,
    };
    terminal_identity_parts(short, &[], window)
}

fn apply_adoptions(
    calls: &mut [Call],
    pools: &Pools,
    mut members: HashMap<u32, HashSet<u32>>,
    legacy: &[(u32, u32)],
    terminal: &[(u32, u32)],
) -> (usize, usize, usize) {
    let mut parent: HashMap<u32, u32> = HashMap::new();
    let mut locations: HashMap<u32, Vec<usize>> = HashMap::new();
    for (index, call) in calls.iter().enumerate() {
        locations.entry(call.cluster).or_default().push(index);
    }
    let mut pair_cache = HashMap::new();
    let mut merged = 0;
    let mut terminal_merged = 0;
    let mut blocked = 0;
    for (child, parent_id, guarded) in legacy
        .iter()
        .map(|&(child, parent)| (child, parent, false))
        .chain(terminal.iter().map(|&(child, parent)| (child, parent, true)))
    {
        let root = find(&mut parent, parent_id);
        let child_root = find(&mut parent, child);
        if child_root == root {
            continue;
        }
        // Pairwise-disjoint proposals can still form a chain that joins two
        // copies in one sample. Check the accumulated groups, not just inputs.
        if !members[&child_root].is_disjoint(&members[&root]) {
            if guarded || !cooccurring_loci_compatible(
                &locations[&child_root], &locations[&root], &members[&child_root],
                &members[&root], calls, pools, &mut pair_cache,
            ) {
                blocked += 1;
                continue;
            }
        }
        let child_samples = members.remove(&child_root).expect("child root exists");
        members.get_mut(&root).expect("parent root exists").extend(child_samples);
        parent.insert(child_root, root);
        let child_locations = locations.remove(&child_root).expect("child locations exist");
        locations.get_mut(&root).expect("parent locations exist").extend(child_locations);
        merged += 1;
        terminal_merged += usize::from(guarded);
    }
    for call in calls.iter_mut() {
        call.cluster = find(&mut parent, call.cluster);
    }
    if merged > 0 {
        coarsen_contexts(calls);
        retain_loci(calls, pools);
    }
    (merged, terminal_merged, blocked)
}

fn coarsen_contexts(calls: &mut [Call]) {
    let mut first_context = HashMap::new();
    let mut parents = HashMap::new();
    for call in calls.iter() {
        let first = *first_context.entry(call.cluster).or_insert(call.context);
        let a = find(&mut parents, first);
        let b = find(&mut parents, call.context);
        if a != b {
            parents.insert(a.max(b), a.min(b));
        }
    }
    // Propagate sequence-family aliases to all descendants, not ortholog IDs.
    for call in calls {
        call.context = find(&mut parents, call.context);
    }
}

fn retain_loci(calls: &mut [Call], pools: &Pools) {
    for (id, locus) in physical_loci(calls, pools).into_iter().enumerate() {
        for call in &mut calls[locus.first..=locus.last] {
            call.locus = Some(id as u32);
        }
    }
}

pub fn membership_loci(calls: &[Call], pools: &Pools) -> Vec<u32> {
    let mut ids = vec![0; calls.len()];
    for (id, locus) in physical_loci(calls, pools).into_iter().enumerate() {
        ids[locus.first..=locus.last].fill(id as u32);
    }
    ids
}

fn cooccurring_loci_compatible(
    left: &[usize],
    right: &[usize],
    left_samples: &HashSet<u32>,
    right_samples: &HashSet<u32>,
    calls: &[Call],
    pools: &Pools,
    cache: &mut HashMap<(u32, u32), bool>,
) -> bool {
    let template_id = left.iter().chain(right).map(|&index| calls[index].protein)
        .max_by(|a, b| {
            pools.protein[*a as usize].len().cmp(&pools.protein[*b as usize].len())
                .then(b.cmp(a))
        }).expect("merge groups contain calls");
    let template = &pools.protein[template_id as usize];
    for &sample in left_samples.intersection(right_samples) {
        let mut indexes: Vec<usize> = left.iter().chain(right)
            .copied().filter(|&index| calls[index].sample == sample).collect();
        indexes.sort_unstable();
        if !indexes.windows(2).all(|pair| same_physical_locus(pair[0], pair[1], calls, pools, template, cache)) {
            return false;
        }
    }
    true
}

struct MergeTarget {
    cluster: u32,
    terminal: bool,
}

fn best_merge_target(
    short_id: u32,
    reps: &HashMap<u32, Vec<u8>>,
    members: &HashMap<u32, HashSet<u32>>,
    flanks: &HashMap<u32, (Option<u32>, Option<u32>)>,
    kmers: &HashMap<u32, HashSet<&[u8]>>,
    candidates: &[u32],
) -> Option<MergeTarget> {
    let short = reps.get(&short_id)?;
    let mut legacy_target = None;
    let mut legacy_matches = 0;
    let mut terminal_target = None;
    let mut terminal_matches = 0;
    for &long_id in candidates {
        if long_id == short_id {
            continue;
        }
        let Some(long) = reps.get(&long_id) else {
            continue;
        };
        if short.len() >= long.len() {
            continue;
        }
        if short.len() * 10 > long.len() * 8 {
            let Some(&(Some(left), Some(right))) = flanks.get(&short_id) else {
                continue;
            };
            if flanks.get(&long_id) != Some(&(Some(left), Some(right)))
                || !members[&short_id].is_disjoint(&members[&long_id])
            {
                continue;
            }
            if terminal_identity(short, long)
                .is_some_and(|identity| identity + f64::EPSILON >= TERMINAL_IDENTITY)
            {
                terminal_target = Some(long_id);
                terminal_matches += 1;
            }
            continue;
        }
        if cooccurrence(members.get(&short_id), members.get(&long_id)) > MAX_COOCCUR {
            continue;
        }
        if !shares_kmer_cached(short_id, long_id, short.len(), long.len(), kmers) {
            continue;
        }
        let Some(fit) = fit_short_in_long(short, long) else {
            continue;
        };
        if fit.identity < CONTAIN_IDENTITY || fit.coverage < CONTAIN_COVERAGE || fit.internal_gap > 5
        {
            continue;
        }
        legacy_target = Some(long_id);
        legacy_matches += 1;
    }
    if terminal_matches == 1 {
        terminal_target.map(|cluster| MergeTarget { cluster, terminal: true })
    } else if terminal_matches == 0 && legacy_matches == 1 {
        legacy_target.map(|cluster| MergeTarget { cluster, terminal: false })
    } else {
        None
    }
}

fn terminal_identity(short: &[u8], long: &[u8]) -> Option<f64> {
    if short.is_empty() || short.len() >= long.len() {
        return None;
    }
    terminal_identity_parts(short, &[], long)
}

fn abutting_pieces(
    calls: &[Call],
    reps: &HashMap<u32, Vec<u8>>,
    members: &HashMap<u32, HashSet<u32>>,
) -> Vec<(u32, u32)> {
    let mut gaps: HashMap<(u32, u32), Vec<i32>> = HashMap::new();
    for pair in successive(calls) {
        let gap = pair.gap;
        if gap.abs() > ABUT_GAP {
            continue;
        }
        let key = if pair.left <= pair.right {
            (pair.left, pair.right)
        } else {
            (pair.right, pair.left)
        };
        gaps.entry(key).or_default().push(gap);
    }
    // Sorted so the merge list below is built in the same order on every
    // run rather than following this `HashMap`'s randomized iteration
    // order; each pair's decision is independent and read-only, so the
    // actual check can still run in parallel.
    let mut keys: Vec<(u32, u32)> = gaps.keys().copied().collect();
    keys.sort_unstable();
    keys.par_iter()
        .filter_map(|&(left, right)| {
            let gap_list = &gaps[&(left, right)];
            if left == right || gap_list.len() < 2 {
                return None;
            }
            let median = {
                let mut sorted = gap_list.clone();
                sorted.sort();
                sorted[sorted.len() / 2]
            };
            if median.abs() > ABUT_GAP {
                return None;
            }
            if cooccurrence(members.get(&left), members.get(&right)) < 0.8 {
                return None;
            }
            let left_prot = reps.get(&left)?;
            let right_prot = reps.get(&right)?;
            let identity = global_identity(left_prot, right_prot)?;
            if identity >= 0.4 {
                return None;
            }
            // Keep the longer piece's id so the reference sees both alleles.
            let child = if left_prot.len() <= right_prot.len() {
                left
            } else {
                right
            };
            let parent = if child == left { right } else { left };
            Some((child, parent))
        })
        .collect()
}

pub fn refind(calls: &mut Vec<Call>, pools: &mut Pools, contigs: &[Vec<String>]) -> usize {
    sort_calls(calls);
    let loci = physical_loci(calls, pools);
    let contexts = locus_contexts(calls, &loci, |call| call.cluster);
    let mut between: HashMap<(u32, u32), HashMap<u32, u32>> = HashMap::new();
    for (locus, context) in loci.iter().zip(contexts) {
        if let (Some(left), Some(right)) = context {
            between
                .entry((left, right))
                .or_default()
                .entry(calls[locus.first].cluster)
                .and_modify(|count| *count += 1)
                .or_insert(1);
        }
    }
    let target_of: HashMap<(u32, u32), u32> = between
        .iter()
        .filter_map(|(flanks, middles)| {
            // `max_by_key` breaks a tie by returning the last-seen element,
            // and `middles` is a `HashMap` with randomized iteration order;
            // breaking a tied vote count by the (otherwise arbitrary) gene
            // id keeps the choice fixed by the data instead of by run.
            let (&gene, _) = middles
                .iter()
                .max_by(|a, b| a.1.cmp(b.1).then(b.0.cmp(a.0)))?;
            let count = middles[&gene];
            if count < 2 {
                return None;
            }
            Some((*flanks, gene))
        })
        .collect();
    if target_of.is_empty() {
        return 0;
    }
    let reps = representatives(calls, pools);
    let context_of: HashMap<u32, u32> =
        calls.iter().map(|call| (call.cluster, call.context)).collect();
    let mut present: HashMap<u32, HashSet<u32>> = HashMap::new();
    for call in calls.iter() {
        present.entry(call.sample).or_default().insert(call.cluster);
    }
    // Every sample's own calls are contiguous (sorted by sample first), and
    // a gap can only ever be filled with a gene that sample's own earlier
    // annotations and already-refound genes brought into scope - never by
    // another sample's. So each sample's scan is independent of every other
    // sample's and can run on its own thread; only interning the results
    // into the shared `Pools` below has to stay sequential.
    let mut sample_ranges: Vec<(u32, usize, usize)> = Vec::new();
    let mut index = 0;
    while index < calls.len() {
        let sample = calls[index].sample;
        let mut end = index + 1;
        while end < calls.len() && calls[end].sample == sample {
            end += 1;
        }
        sample_ranges.push((sample, index, end));
        index = end;
    }
    let found: Vec<FoundGene> = sample_ranges
        .par_iter()
        .flat_map(|&(sample, start, end)| {
            let local_present = present.get(&sample).cloned().unwrap_or_default();
            refind_in_sample(&calls[start..end], sample, &target_of, &reps, contigs, local_present)
        })
        .collect();
    let n = found.len();
    for gene in found {
        calls.push(Call {
            sample: gene.sample,
            contig: gene.contig,
            start: gene.start,
            end: gene.end,
            strand: gene.strand,
            dna: pools.intern_dna(gene.dna),
            protein: pools.intern_protein(gene.protein),
            cluster: gene.cluster,
            context: context_of[&gene.cluster],
            locus: None,
            source: None,
            fusion_piece: false,
            symbol: String::new(),
            product: String::new(),
        });
    }
    n
}

/// A gene found in an intergenic gap, before its DNA/protein are interned
/// into the shared `Pools`. Kept separate from `Call` so the parallel,
/// per-sample search phase never has to touch `Pools` (which needs
/// exclusive access to assign ids) and the sequential interning pass that
/// follows is as small as possible.
struct FoundGene {
    sample: u32,
    contig: u32,
    start: u32,
    end: u32,
    strand: i8,
    dna: String,
    protein: Vec<u8>,
    cluster: u32,
}

/// Scan one sample's calls (already contiguous, sorted by contig then
/// start) for intergenic gaps that match a known flank pair, same logic as
/// the old single-threaded scan but scoped to one sample and its own local
/// `present` set, so a gene found via one gap cannot be refound again via a
/// later gap in the same sample - the only ordering guarantee the original
/// sequential version relied on.
fn refind_in_sample(
    sample_calls: &[Call],
    sample: u32,
    target_of: &HashMap<(u32, u32), u32>,
    reps: &HashMap<u32, Vec<u8>>,
    contigs: &[Vec<String>],
    mut present: HashSet<u32>,
) -> Vec<FoundGene> {
    let mut found = Vec::new();
    let mut index = 0;
    while index < sample_calls.len() {
        let contig = sample_calls[index].contig;
        let mut end = index + 1;
        while end < sample_calls.len() && sample_calls[end].contig == contig {
            end += 1;
        }
        let window = &sample_calls[index..end];
        for pair in window.windows(2) {
            let gap = pair[1].start as i32 - pair[0].end as i32 - 1;
            if !(60..=12_000).contains(&gap) {
                continue;
            }
            let Some(&target) = target_of.get(&(pair[0].cluster, pair[1].cluster)) else {
                continue;
            };
            if present.contains(&target) {
                continue;
            }
            let Some(target_protein) = reps.get(&target) else {
                continue;
            };
            let sequence = contig_sequence(contigs, sample, contig);
            let from = pair[0].end as usize;
            let to = (pair[1].start as usize).saturating_sub(1).min(sequence.len());
            if from >= to || to > sequence.len() {
                continue;
            }
            let gap_seq = &sequence[from..to];
            let Some(orf) = best_orf(gap_seq, target_protein) else {
                continue;
            };
            let protein = translate(&orf.dna);
            found.push(FoundGene {
                sample,
                contig,
                start: pair[0].end + orf.offset + 1,
                end: pair[0].end + orf.offset + orf.dna.len() as u32,
                strand: orf.strand,
                dna: orf.dna,
                protein,
                cluster: target,
            });
            present.insert(target);
        }
        index = end;
    }
    found
}

struct FoundOrf {
    dna: String,
    offset: u32,
    strand: i8,
}

fn best_orf(gap: &str, target: &[u8]) -> Option<FoundOrf> {
    let reverse = crate::genetic_code::reverse_complement(gap);
    // Bounds only. The DNA string is built for an ORF that clears 0.80
    // identity, and equal scores keep the earlier strand, frame, and codon.
    let mut spans: Vec<(i8, usize, usize, u32)> = Vec::new();
    for (sequence, strand) in [(gap, 1i8), (reverse.as_str(), -1)] {
        for frame in 0..3 {
            if sequence.len() <= frame {
                continue;
            }
            let frame_seq = &sequence[frame..];
            let mut codon = 0;
            while codon + 1 < frame_seq.len() / 3 {
                let start = codon * 3;
                let codon_bytes = &frame_seq.as_bytes()[start..start + 3];
                if !is_start(codon_bytes) {
                    codon += 1;
                    continue;
                }
                let mut stop_at = None;
                let mut next = codon + 1;
                while next < frame_seq.len() / 3 {
                    let at = next * 3;
                    if is_stop(&frame_seq.as_bytes()[at..at + 3]) {
                        stop_at = Some(next);
                        break;
                    }
                    next += 1;
                }
                let Some(stop_at) = stop_at else {
                    codon += 1;
                    continue;
                };
                let end = (stop_at + 1) * 3;
                let offset = if strand > 0 {
                    (frame + start) as u32
                } else {
                    // Reverse-complement coordinate: the ORF ends `frame+end`
                    // bases before the end of the forward gap.
                    let span = frame + end;
                    gap.len().saturating_sub(span) as u32
                };
                spans.push((strand, frame + start, frame + end, offset));
                codon = stop_at + 1;
            }
        }
    }
    let score = |&(strand, start, end, offset): &(i8, usize, usize, u32)| -> Option<(f64, FoundOrf)> {
        let sequence = if strand > 0 { gap } else { reverse.as_str() };
        let protein = translate_bytes(&sequence.as_bytes()[start..end]);
        let identity = if protein.len().abs_diff(target.len()) * 5 <= target.len() {
            global_identity(&protein, target).unwrap_or(0.0)
        } else if protein.len() <= target.len() {
            fit_short_in_long(&protein, target)
                .map(|fit| fit.identity)
                .unwrap_or(0.0)
        } else {
            fit_short_in_long(target, &protein)
                .map(|fit| fit.identity)
                .unwrap_or(0.0)
        };
        if identity < 0.8 {
            return None;
        }
        Some((
            identity,
            FoundOrf {
                dna: sequence[start..end].to_string(),
                offset,
                strand,
            },
        ))
    };
    let scored: Vec<Option<(f64, FoundOrf)>> = if spans.len() >= 4 {
        spans.par_iter().map(score).collect()
    } else {
        spans.iter().map(score).collect()
    };
    let mut best: Option<(f64, FoundOrf)> = None;
    for found in scored.into_iter().flatten() {
        if best.as_ref().map(|(identity, _)| found.0 > *identity).unwrap_or(true) {
            best = Some(found);
        }
    }
    best.map(|(_, orf)| orf)
}

fn slice_coding(dna: &str, aa_start: usize, aa_end: usize) -> Option<String> {
    if aa_end <= aa_start {
        return None;
    }
    let bytes = dna.as_bytes();
    let start = aa_start * 3;
    let mut end = aa_end * 3;
    if start >= bytes.len() {
        return None;
    }
    end = end.min(bytes.len() - bytes.len() % 3);
    if end <= start {
        return None;
    }
    let mut sequence = String::from_utf8_lossy(&bytes[start..end]).into_owned();
    if sequence.len() >= 3 && !is_stop(&sequence.as_bytes()[sequence.len() - 3..]) {
        sequence.push_str("TAA");
    }
    Some(sequence)
}

/// 4-mer sets for every representative protein, keyed by cluster id and
/// built once up front. `split_fusions` and `merge_fragments` each test this
/// for every pair of clusters; rebuilding the longer protein's k-mer set
/// from scratch inside that O(clusters^2) loop (the old `shares_kmer`) meant
/// redoing the same O(length) work up to once per pair, which dominates the
/// whole pass once a pangenome has more than a few hundred families.
fn kmer_sets(reps: &HashMap<u32, Vec<u8>>) -> HashMap<u32, HashSet<&[u8]>> {
    let ids: Vec<u32> = reps
        .iter()
        .filter(|(_, protein)| protein.len() >= 4)
        .map(|(&id, _)| id)
        .collect();
    let built: Vec<(u32, HashSet<&[u8]>)> = ids
        .par_iter()
        .map(|&id| {
            let protein = &reps[&id];
            (id, protein.windows(4).collect())
        })
        .collect();
    built.into_iter().collect()
}

fn shares_kmer_cached(
    short_id: u32,
    long_id: u32,
    short_len: usize,
    long_len: usize,
    kmer_sets: &HashMap<u32, HashSet<&[u8]>>,
) -> bool {
    if short_len < 4 || long_len < 4 {
        return true;
    }
    match (kmer_sets.get(&short_id), kmer_sets.get(&long_id)) {
        (Some(a), Some(b)) => {
            if a.len() <= b.len() {
                a.iter().any(|window| b.contains(window))
            } else {
                b.iter().any(|window| a.contains(window))
            }
        }
        _ => false,
    }
}

fn representatives(calls: &[Call], pools: &Pools) -> HashMap<u32, Vec<u8>> {
    representative_ids(calls, pools).into_iter()
        .map(|(cluster, protein)| (cluster, pools.protein[protein as usize].clone()))
        .collect()
}

fn representative_ids(calls: &[Call], pools: &Pools) -> HashMap<u32, u32> {
    let mut best = HashMap::new();
    for call in calls {
        let slot = best.entry(call.cluster).or_insert(call.protein);
        if pools.protein[call.protein as usize].len() > pools.protein[*slot as usize].len() {
            *slot = call.protein;
        }
    }
    best
}

fn sample_counts(calls: &[Call]) -> HashMap<u32, u32> {
    sample_sets(calls)
        .into_iter()
        .map(|(cluster, samples)| (cluster, samples.len() as u32))
        .collect()
}

fn sample_sets(calls: &[Call]) -> HashMap<u32, HashSet<u32>> {
    let mut sets: HashMap<u32, HashSet<u32>> = HashMap::new();
    for call in calls {
        sets.entry(call.cluster).or_default().insert(call.sample);
    }
    sets
}

fn cooccurrence(left: Option<&HashSet<u32>>, right: Option<&HashSet<u32>>) -> f64 {
    let (Some(left), Some(right)) = (left, right) else {
        return 0.0;
    };
    let union = left.union(right).count();
    if union == 0 {
        0.0
    } else {
        left.intersection(right).count() as f64 / union as f64
    }
}

#[cfg(test)]
fn modal_flanks(calls: &[Call], pools: &Pools) -> HashMap<u32, (Option<u32>, Option<u32>)> {
    let loci = physical_loci(calls, pools);
    let contexts = locus_contexts(calls, &loci, |call| call.cluster);
    modal_flanks_of(calls, &loci, &contexts)
}

fn modal_flanks_of(
    calls: &[Call], loci: &[Locus], contexts: &[(Option<u32>, Option<u32>)],
) -> HashMap<u32, (Option<u32>, Option<u32>)> {
    let mut votes: HashMap<u32, HashMap<(Option<u32>, Option<u32>), u32>> = HashMap::new();
    for (locus, &context) in loci.iter().zip(contexts) {
        *votes
            .entry(calls[locus.first].cluster)
            .or_default()
            .entry(context)
            .or_default() += 1;
    }
    votes
        .into_iter()
        .map(|(cluster, contexts)| {
            let flank = contexts
                .into_iter()
                .max_by(|a, b| a.1.cmp(&b.1).then(b.0.cmp(&a.0)))
                .map(|(context, _)| context)
                .unwrap_or((None, None));
            (cluster, flank)
        })
        .collect()
}

struct MergeIndex {
    left: HashMap<u32, Vec<u32>>,
    right: HashMap<u32, Vec<u32>>,
    adjacent: HashMap<u32, Vec<u32>>,
}

impl MergeIndex {
    fn new(
        flanks: &HashMap<u32, (Option<u32>, Option<u32>)>,
        adjacency: &HashSet<(u32, u32)>,
    ) -> Self {
        let mut index = Self {
            left: HashMap::new(),
            right: HashMap::new(),
            adjacent: HashMap::new(),
        };
        for (&id, &(left, right)) in flanks {
            if let Some(left) = left {
                index.left.entry(left).or_default().push(id);
            }
            if let Some(right) = right {
                index.right.entry(right).or_default().push(id);
            }
        }
        for &(id, neighbor) in adjacency {
            index.adjacent.entry(id).or_default().push(neighbor);
        }
        index
    }

    fn candidates(
        &self,
        id: u32,
        flanks: &HashMap<u32, (Option<u32>, Option<u32>)>,
    ) -> Vec<u32> {
        let mut candidates = self.adjacent.get(&id).cloned().unwrap_or_default();
        if let Some(&(left, right)) = flanks.get(&id) {
            if let Some(left) = left.and_then(|flank| self.left.get(&flank)) {
                candidates.extend_from_slice(left);
            }
            if let Some(right) = right.and_then(|flank| self.right.get(&flank)) {
                candidates.extend_from_slice(right);
            }
        }
        candidates.sort_unstable();
        candidates.dedup();
        candidates.retain(|&candidate| candidate != id);
        candidates
    }
}

/// Every `(cluster_a, cluster_b)` pair observed abutting in at least one
/// sample, in both directions. `successive(calls)` doesn't depend on which
/// pair is being tested, so building this once instead of re-deriving it
/// inside `merge_fragments`'s `ids x ids` loop turns an O(calls.len()) rescan
/// per candidate pair into a single O(1) lookup per pair.
fn adjacency_pairs(calls: &[Call]) -> HashSet<(u32, u32)> {
    let mut pairs = HashSet::new();
    for pair in successive(calls) {
        pairs.insert((pair.left, pair.right));
        pairs.insert((pair.right, pair.left));
    }
    pairs
}

struct Pair {
    left: u32,
    right: u32,
    gap: i32,
}

fn successive(calls: &[Call]) -> Vec<Pair> {
    let mut pairs = Vec::new();
    for window in calls.windows(2) {
        if window[0].sample != window[1].sample || window[0].contig != window[1].contig {
            continue;
        }
        let gap = window[1].start as i32 - window[0].end as i32 - 1;
        if !(-15..=NEIGHBOR_GAP).contains(&gap) {
            continue;
        }
        pairs.push(Pair {
            left: window[0].cluster,
            right: window[1].cluster,
            gap,
        });
    }
    pairs
}

struct Locus {
    first: usize,
    last: usize,
    end: u32,
}

fn same_physical_locus(
    left: usize,
    right: usize,
    calls: &[Call],
    pools: &Pools,
    template: &[u8],
    cache: &mut HashMap<(u32, u32), bool>,
) -> bool {
    let a = &calls[left];
    let b = &calls[right];
    if right != left + 1 || a.sample != b.sample || a.contig != b.contig || a.strand != b.strand {
        return false;
    }
    if a.locus.is_some() && a.locus == b.locus {
        return true;
    }
    let overlap = i64::from(a.end.min(b.end)) - i64::from(a.start.max(b.start)) + 1;
    let shorter = (i64::from(a.end) - i64::from(a.start) + 1)
        .min(i64::from(b.end) - i64::from(b.start) + 1);
    if shorter > 0 && overlap * 100 >= shorter * 95 {
        return true;
    }
    let gap = i64::from(b.start) - i64::from(a.end) - 1;
    if gap.abs() > i64::from(ABUT_GAP) {
        return false;
    }
    let a_protein = &pools.protein[a.protein as usize];
    let b_protein = &pools.protein[b.protein as usize];
    let (head, tail) = if a.strand >= 0 {
        (a_protein.as_slice(), b_protein.as_slice())
    } else {
        (b_protein.as_slice(), a_protein.as_slice())
    };
    // Repetitive fragments can resemble each other. A longer observed allele
    // can prove their join, whereas two full tandem copies exceed its length.
    if template.len() > head.len().max(tail.len())
        && (head.len() + tail.len()) as f64 >= 0.95 * template.len() as f64
        && terminal_identity_parts(head, tail, template)
            .is_some_and(|identity| identity + f64::EPSILON >= TERMINAL_IDENTITY)
    {
        return true;
    }
    if a.protein == b.protein {
        return false;
    }
    let pair = (a.protein.min(b.protein), a.protein.max(b.protein));
    *cache.entry(pair).or_insert_with(|| {
        global_identity(&pools.protein[a.protein as usize], &pools.protein[b.protein as usize])
            .is_some_and(|identity| identity < 0.4)
    })
}

fn physical_loci(calls: &[Call], pools: &Pools) -> Vec<Locus> {
    let representatives = representative_ids(calls, pools);
    let mut loci: Vec<Locus> = Vec::new();
    let mut cache = HashMap::new();
    for (index, call) in calls.iter().enumerate() {
        if let Some(previous) = loci.last_mut() {
            if calls[previous.first].cluster == call.cluster
                && same_physical_locus(
                    previous.last, index, calls, pools,
                    &pools.protein[representatives[&call.cluster] as usize], &mut cache,
                )
            {
                previous.last = index;
                previous.end = previous.end.max(call.end);
                continue;
            }
        }
        loci.push(Locus { first: index, last: index, end: call.end });
    }
    loci
}

fn locus_contexts(
    calls: &[Call], loci: &[Locus], family: impl Fn(&Call) -> u32,
) -> Vec<(Option<u32>, Option<u32>)> {
    let mut contexts = vec![(None, None); loci.len()];
    let mut index = 0;
    while index < loci.len() {
        let mut end = index + 1;
        let first = &calls[loci[index].first];
        while end < loci.len()
            && calls[loci[end].first].sample == first.sample
            && calls[loci[end].first].contig == first.contig
        {
            end += 1;
        }
        for cursor in index..end {
            let prev = if cursor > index {
                let gap = i64::from(calls[loci[cursor].first].start) - i64::from(loci[cursor - 1].end) - 1;
                if (-15..=i64::from(NEIGHBOR_GAP)).contains(&gap) {
                    Some(family(&calls[loci[cursor - 1].first]))
                } else {
                    None
                }
            } else {
                None
            };
            let next = if cursor + 1 < end {
                let gap = i64::from(calls[loci[cursor + 1].first].start) - i64::from(loci[cursor].end) - 1;
                if (-15..=i64::from(NEIGHBOR_GAP)).contains(&gap) {
                    Some(family(&calls[loci[cursor + 1].first]))
                } else {
                    None
                }
            } else {
                None
            };
            contexts[cursor] = (prev, next);
        }
        index = end;
    }
    contexts
}


fn contig_sequence(contigs: &[Vec<String>], sample: u32, contig: u32) -> &str {
    contigs
        .get(sample as usize)
        .and_then(|sample| sample.get(contig as usize))
        .map(String::as_str)
        .unwrap_or("")
}

fn find(parent: &mut HashMap<u32, u32>, id: u32) -> u32 {
    let mut root = id;
    while let Some(&next) = parent.get(&root) {
        root = next;
    }
    let mut current = id;
    while let Some(&next) = parent.get(&current) {
        parent.insert(current, root);
        current = next;
    }
    root
}

fn sort_calls(calls: &mut [Call]) {
    calls.sort_by(|a, b| {
        a.sample
            .cmp(&b.sample)
            .then(a.contig.cmp(&b.contig))
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
    });
}

pub fn compress_clusters(calls: &mut [Call]) {
    sort_calls(calls);
    let mut map = HashMap::new();
    let mut next = 0u32;
    for call in calls.iter_mut() {
        let id = *map.entry(call.cluster).or_insert_with(|| {
            let id = next;
            next += 1;
            id
        });
        call.cluster = id;
    }
}

/// Drop families that are not usable genes.
///
/// A call with a broken frame, an internal stop, or fewer than 90 bp may sit
/// inside a family that also has an intact allele. A family whose every call
/// is one of those is an annotation artifact and is removed.
///
/// Separately, walk inward from both ends of every contig. A family present
/// in fewer than `max(2, ceil(0.05 * samples))` genomes is removed while every
/// one of its copies is still on that walk: it is not between two neighboring
/// genes, and it does not itself clear the sample threshold. The walk repeats,
/// because deleting the outer gene exposes the next one. A private contig of
/// such families disappears with them. A rare gene that sits between two
/// neighbors is left in place. When no family in the collection clears the
/// threshold, nothing is judged a tail: there is no conserved backbone to
/// walk toward.
pub fn screen_unsupported_families(calls: &mut Vec<Call>, pools: &Pools, n_samples: usize) {
    let invalid = families_without_intact_allele(calls, pools);
    if !invalid.is_empty() {
        record_removed(calls, pools, &invalid, "no_intact_allele");
        calls.retain(|call| !invalid.contains(&call.cluster));
        progress!(
            "rustroo: removed {} families with no intact allele",
            invalid.len()
        );
    }
    let mut removed = 0usize;
    loop {
        if calls.is_empty() {
            break;
        }
        sort_calls(calls);
        let loci = physical_loci(calls, pools);
        let tails = contig_end_tail_families(calls, &loci, n_samples);
        if tails.is_empty() {
            break;
        }
        record_removed(calls, pools, &tails, "contig_end");
        removed += tails.len();
        calls.retain(|call| !tails.contains(&call.cluster));
    }
    if removed > 0 {
        progress!("rustroo: removed {removed} unsupported contig-end families");
    }
}

fn record_removed(calls: &[Call], pools: &Pools, clusters: &HashSet<u32>, reason: &str) {
    if !log::enabled() || clusters.is_empty() {
        return;
    }
    struct Row {
        symbol: String,
        product: String,
        samples: HashSet<u32>,
        copies: u32,
        longest: usize,
    }
    let mut rows: HashMap<u32, Row> = HashMap::new();
    for call in calls.iter().filter(|call| clusters.contains(&call.cluster)) {
        let row = rows.entry(call.cluster).or_insert_with(|| Row {
            symbol: String::new(),
            product: String::new(),
            samples: HashSet::new(),
            copies: 0,
            longest: 0,
        });
        if row.symbol.is_empty() {
            row.symbol.clone_from(&call.symbol);
        }
        if row.product.is_empty() {
            row.product.clone_from(&call.product);
        }
        row.samples.insert(call.sample);
        row.copies += 1;
        row.longest = row.longest.max(pools.dna[call.dna as usize].len());
    }
    log::record("removed\treason\tcluster\tsymbol\tproduct\tsamples\tcopies\tlongest_bp");
    let mut ordered: Vec<_> = rows.into_iter().collect();
    ordered.sort_by(|left, right| left.1.symbol.cmp(&right.1.symbol).then(left.0.cmp(&right.0)));
    for (cluster, row) in ordered {
        log::record(&format!(
            "removed\t{reason}\t{cluster}\t{}\t{}\t{}\t{}\t{}",
            log_field(&row.symbol),
            log_field(&row.product),
            row.samples.len(),
            row.copies,
            row.longest,
        ));
    }
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

fn families_without_intact_allele(calls: &[Call], pools: &Pools) -> HashSet<u32> {
    let mut intact = HashSet::new();
    let mut seen = HashSet::new();
    for call in calls {
        seen.insert(call.cluster);
        if !sequence_cannot_found_a_family(&pools.dna[call.dna as usize]) {
            intact.insert(call.cluster);
        }
    }
    seen.difference(&intact).copied().collect()
}

/// Out of frame, an internal stop, or shorter than 90 bp. A truncation that
/// stays in frame and has no internal stop can still found a family when
/// synteny did not adopt it.
fn sequence_cannot_found_a_family(dna: &str) -> bool {
    let status = orf_status(dna);
    !status.in_frame || status.internal_stops > 0 || dna.len() < 90
}

fn support_threshold(n_samples: usize) -> u32 {
    let fraction = (n_samples as f64 * 0.05).ceil() as u32;
    fraction.max(2)
}

/// Prodigal's empirical balance point. Below this length a positive coding
/// score in a single genome is about as likely to be false as real, so length
/// itself is not evidence. The constant is a label boundary, not a filter.
const LENGTH_PRIOR_BP: usize = 250;

pub const LENGTH_SUFFICIENT: &str = "length_sufficient";
pub const CONSERVED_SHORT: &str = "conserved_short";
pub const UNCONFIRMED_SHORT: &str = "unconfirmed_short";

/// Classify the length uncertainty of one finished family.
///
/// A reference of at least 250 bp needs nothing else from the collection.
/// A shorter gene is conserved when enough genomes carry it, or when any copy
/// sits between two neighboring genes. Otherwise the call is the short,
/// single-sample case the length prior cannot settle.
pub fn length_evidence(
    reference_bp: usize,
    samples_present: u32,
    n_samples: usize,
    between_neighbors: bool,
) -> &'static str {
    if reference_bp >= LENGTH_PRIOR_BP {
        return LENGTH_SUFFICIENT;
    }
    if samples_present >= support_threshold(n_samples) || between_neighbors {
        CONSERVED_SHORT
    } else {
        UNCONFIRMED_SHORT
    }
}

pub fn families_bounded_by_neighbors(calls: &mut [Call], pools: &Pools) -> HashSet<u32> {
    sort_calls(calls);
    let loci = physical_loci(calls, pools);
    families_between_neighbors(calls, &loci)
}

fn weak_flank_threshold(n_samples: usize) -> usize {
    let fraction = (n_samples as f64 * 0.01).ceil() as usize;
    fraction.max(2)
}

fn contig_end_tail_families(calls: &[Call], loci: &[Locus], n_samples: usize) -> HashSet<u32> {
    let support = sample_counts(calls);
    let threshold = support_threshold(n_samples);
    if !support.values().any(|count| *count >= threshold) {
        return HashSet::new();
    }
    let internal = families_between_neighbors(calls, loci);
    let mut consumed = HashSet::new();
    let mut index = 0;
    while index < loci.len() {
        let mut end = index + 1;
        let sample = calls[loci[index].first].sample;
        let contig = calls[loci[index].first].contig;
        while end < loci.len()
            && calls[loci[end].first].sample == sample
            && calls[loci[end].first].contig == contig
        {
            end += 1;
        }
        consume_tail(&mut consumed, calls, loci, index..end, true, &support, &internal, threshold);
        consume_tail(&mut consumed, calls, loci, index..end, false, &support, &internal, threshold);
        index = end;
    }
    // A family stays if any copy was not on a tail walk. Dropping the family
    // would also delete that internal or gap-isolated copy.
    let mut total: HashMap<u32, usize> = HashMap::new();
    let mut hit: HashMap<u32, usize> = HashMap::new();
    for (locus_index, locus) in loci.iter().enumerate() {
        let family = calls[locus.first].cluster;
        *total.entry(family).or_default() += 1;
        if consumed.contains(&locus_index) {
            *hit.entry(family).or_default() += 1;
        }
    }
    total
        .into_iter()
        .filter(|(family, count)| hit.get(family).copied().unwrap_or(0) == *count)
        .map(|(family, _)| family)
        .collect()
}

fn consume_tail(
    consumed: &mut HashSet<usize>,
    calls: &[Call],
    loci: &[Locus],
    range: std::ops::Range<usize>,
    forward: bool,
    support: &HashMap<u32, u32>,
    internal: &HashSet<u32>,
    threshold: u32,
) {
    let cursors: Vec<usize> = if forward {
        range.collect()
    } else {
        range.rev().collect()
    };
    let mut previous: Option<usize> = None;
    for cursor in cursors {
        if let Some(previous) = previous {
            let adjacent = if forward {
                neighbor_gap(calls, loci, previous, cursor)
            } else {
                neighbor_gap(calls, loci, cursor, previous)
            };
            if !adjacent {
                break;
            }
        }
        let family = calls[loci[cursor].first].cluster;
        if support.get(&family).copied().unwrap_or(0) >= threshold || internal.contains(&family) {
            break;
        }
        consumed.insert(cursor);
        previous = Some(cursor);
    }
}

/// A family is between two genes when any copy has a neighbor on both sides
/// within the same gap used for flank evidence. One such copy protects every
/// terminal copy of that family.
fn families_between_neighbors(calls: &[Call], loci: &[Locus]) -> HashSet<u32> {
    let mut internal = HashSet::new();
    let mut index = 0;
    while index < loci.len() {
        let mut end = index + 1;
        let sample = calls[loci[index].first].sample;
        let contig = calls[loci[index].first].contig;
        while end < loci.len()
            && calls[loci[end].first].sample == sample
            && calls[loci[end].first].contig == contig
        {
            end += 1;
        }
        for cursor in index..end {
            let prev = cursor > index && neighbor_gap(calls, loci, cursor - 1, cursor);
            let next = cursor + 1 < end && neighbor_gap(calls, loci, cursor, cursor + 1);
            if prev && next {
                internal.insert(calls[loci[cursor].first].cluster);
            }
        }
        index = end;
    }
    internal
}

fn neighbor_gap(calls: &[Call], loci: &[Locus], left: usize, right: usize) -> bool {
    let gap = i64::from(calls[loci[right].first].start) - i64::from(loci[left].end) - 1;
    (-15..=i64::from(NEIGHBOR_GAP)).contains(&gap)
}

/// A neighbor seen in fewer than `max(2, 1% of samples)` genomes is not a
/// different context when the other flank already identifies a single
/// better-supported neighbor. The rare edge then cannot split a paralog or
/// veto adoption of a fragment. Two equally supported alternatives are left
/// alone: there is no basis for choosing one.
fn stabilize_contexts(
    calls: &[Call],
    loci: &[Locus],
    contexts: &[(Option<u32>, Option<u32>)],
    family_of: impl Fn(&Call) -> u32,
) -> Vec<(Option<u32>, Option<u32>)> {
    let n_samples = calls.iter().map(|call| call.sample as usize + 1).max().unwrap_or(0);
    let threshold = weak_flank_threshold(n_samples);
    let mut support: HashMap<(u32, Option<u32>, Option<u32>), HashSet<u32>> = HashMap::new();
    for (locus, context) in loci.iter().zip(contexts.iter()) {
        let family = family_of(&calls[locus.first]);
        support
            .entry((family, context.0, context.1))
            .or_default()
            .insert(calls[locus.first].sample);
    }
    let count = |family, prev, next| {
        support
            .get(&(family, prev, next))
            .map(|samples| samples.len())
            .unwrap_or(0)
    };
    contexts
        .iter()
        .enumerate()
        .map(|(index, context)| {
            let family = family_of(&calls[loci[index].first]);
            if count(family, context.0, context.1) >= threshold {
                return *context;
            }
            let next = supported_alternative(family, context.0, context.1, true, threshold, &support);
            let prev = supported_alternative(family, context.0, context.1, false, threshold, &support);
            // Replacing both sides at once builds a neighborhood this gene was
            // never in. A missing side is also left missing: an absent flank
            // at a contig break is evidence, not a value to invent.
            if prev != context.0 && next != context.1 {
                return *context;
            }
            (prev, next)
        })
        .collect()
}

/// `fix_next` searches for one strong replacement of the downstream flank
/// among contexts that keep the upstream flank. The opposite search keeps the
/// downstream flank. The observed pair is allowed to compete; a lone strong
/// alternative replaces a weak one.
fn supported_alternative(
    family: u32,
    prev: Option<u32>,
    next: Option<u32>,
    fix_next: bool,
    threshold: usize,
    support: &HashMap<(u32, Option<u32>, Option<u32>), HashSet<u32>>,
) -> Option<u32> {
    let observed = if fix_next { next } else { prev };
    if observed.is_none() {
        return None;
    }
    let mut alternatives: Vec<(usize, u32)> = support
        .iter()
        .filter(|((id, pair_prev, pair_next), _)| {
            *id == family && if fix_next { *pair_prev == prev } else { *pair_next == next }
        })
        .filter_map(|((_, pair_prev, pair_next), samples)| {
            let neighbor = if fix_next { *pair_next } else { *pair_prev };
            let neighbor = neighbor?;
            let count = samples.len();
            (count >= threshold).then_some((count, neighbor))
        })
        .collect();
    alternatives.sort_unstable_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
    alternatives.dedup_by(|a, b| a.1 == b.1);
    let observed = if fix_next { next } else { prev };
    match alternatives.as_slice() {
        [(count, neighbor)] if observed != Some(*neighbor) && *count >= threshold => Some(*neighbor),
        _ => observed,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genetic_code::back_translate;

    fn protein(length: usize, seed: u64) -> Vec<u8> {
        let alphabet = b"ACDEFGHIKLMNPQRSTVWY";
        let mut state = seed;
        (0..length)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                alphabet[(state % 20) as usize]
            })
            .collect()
    }

    fn merge(calls: &mut [Call], pools: &Pools) -> bool {
        let mut lengths = Vec::new();
        for call in calls.iter() {
            lengths.resize_with(lengths.len().max(call.sample as usize + 1), Vec::new);
            let sample = &mut lengths[call.sample as usize];
            sample.resize(sample.len().max(call.contig as usize + 1), 0);
            sample[call.contig as usize] = sample[call.contig as usize].max(call.end as usize + 300);
        }
        merge_fragments(calls, pools, &lengths)
    }

    fn target(
        reps: &HashMap<u32, Vec<u8>>,
        members: &HashMap<u32, HashSet<u32>>,
        flanks: &HashMap<u32, (Option<u32>, Option<u32>)>,
    ) -> Option<MergeTarget> {
        let index = MergeIndex::new(flanks, &HashSet::new());
        best_merge_target(
            0, reps, members, flanks, &kmer_sets(reps), &index.candidates(0, flanks),
        )
    }

    #[test]
    fn indexed_candidates_match_exhaustive_context_and_adjacency_filter() {
        let flanks: HashMap<_, _> = (0..64)
            .map(|id| {
                let left = (id % 3 != 0).then_some(id % 7);
                let right = (id % 5 != 0).then_some(id % 11);
                (id, (left, right))
            })
            .collect();
        let adjacency: HashSet<_> = (0..64).map(|id| (id, (id + 17) % 64)).collect();
        let index = MergeIndex::new(&flanks, &adjacency);
        for id in 0..64 {
            let (left, right) = flanks[&id];
            let expected: Vec<_> = (0..64)
                .filter(|&other| {
                    let (other_left, other_right) = flanks[&other];
                    other != id
                        && ((left.is_some() && left == other_left)
                            || (right.is_some() && right == other_right)
                            || adjacency.contains(&(id, other)))
                })
                .collect();
            assert_eq!(index.candidates(id, &flanks), expected);
        }
    }

    #[test]
    fn moderately_shortened_terminal_alleles_merge_at_length_boundaries() {
        let full = protein(100, 42);
        let members = HashMap::from([(0, HashSet::from([0])), (1, HashSet::from([1]))]);
        let flanks = HashMap::from([
            (0, (Some(10), Some(11))), (1, (Some(10), Some(11))),
        ]);
        for length in [81, 90, 94, 95, 99] {
            for short in [&full[..length], &full[100 - length..]] {
                let reps = HashMap::from([(0, short.to_vec()), (1, full.clone())]);
                let chosen = target(&reps, &members, &flanks).expect("terminal allele");
                assert_eq!(chosen.cluster, 1);
                assert!(chosen.terminal, "length {length}");
            }
        }
    }

    #[test]
    fn terminal_merges_require_both_known_flanks_and_disjoint_samples() {
        let full = protein(100, 42);
        let reps = HashMap::from([(0, full[..90].to_vec()), (1, full)]);
        let mut members = HashMap::from([
            (0, HashSet::from([0])), (1, HashSet::from([1, 2, 3, 4, 5])),
        ]);
        for context in [(Some(10), None), (None, Some(11)), (Some(10), Some(12))] {
            let flanks = HashMap::from([
                (0, context), (1, (Some(10), Some(11))),
            ]);
            assert!(target(&reps, &members, &flanks).is_none());
        }
        let flanks = HashMap::from([
            (0, (Some(10), Some(11))), (1, (Some(10), Some(11))),
        ]);
        members.get_mut(&1).unwrap().insert(0);
        assert!(target(&reps, &members, &flanks).is_none());
    }

    #[test]
    fn interior_domains_and_internal_deletions_are_not_terminal_alleles() {
        let full = protein(100, 42);
        let members = HashMap::from([(0, HashSet::from([0])), (1, HashSet::from([1]))]);
        let flanks = HashMap::from([
            (0, (Some(10), Some(11))), (1, (Some(10), Some(11))),
        ]);
        let mut deleted = full[..45].to_vec();
        deleted.extend_from_slice(&full[55..]);
        for short in [full[5..95].to_vec(), deleted] {
            let reps = HashMap::from([(0, short), (1, full.clone())]);
            assert!(target(&reps, &members, &flanks).is_none());
        }
    }

    #[test]
    fn terminal_identity_requires_full_coverage_and_98_percent_identity() {
        let full = protein(100, 42);
        let members = HashMap::from([(0, HashSet::from([0])), (1, HashSet::from([1]))]);
        let flanks = HashMap::from([
            (0, (Some(10), Some(11))), (1, (Some(10), Some(11))),
        ]);
        let mut short = full[..90].to_vec();
        for changes in 1..=2 {
            let at = changes * 20;
            short[at] = if short[at] == b'A' { b'C' } else { b'A' };
            let reps = HashMap::from([(0, short.clone()), (1, full.clone())]);
            assert_eq!(target(&reps, &members, &flanks).is_some(), changes == 1);
        }
        assert_eq!(terminal_identity(&full[2..92], &full), Some(1.0));
        assert_eq!(terminal_identity(&full[8..98], &full), Some(1.0));
        assert_eq!(terminal_identity(&full, &full), None);
        assert_eq!(terminal_identity(b"", &full), None);
    }

    #[test]
    fn competing_terminal_targets_are_not_merged_arbitrarily() {
        let full = protein(100, 42);
        let reps = HashMap::from([
            (0, full[..90].to_vec()), (1, full.clone()), (2, full[..95].to_vec()),
        ]);
        let members = HashMap::from([
            (0, HashSet::from([0])), (1, HashSet::from([1])), (2, HashSet::from([2])),
        ]);
        let flanks = HashMap::from([
            (0, (Some(10), Some(11))), (1, (Some(10), Some(11))),
            (2, (Some(10), Some(11))),
        ]);
        assert!(target(&reps, &members, &flanks).is_none());
    }

    #[test]
    fn competing_legacy_fragment_targets_are_not_chosen_by_sample_support() {
        let full = protein(120, 42);
        let reps = HashMap::from([
            (0, full[..60].to_vec()), (1, full.clone()), (2, full[..110].to_vec()),
        ]);
        let members = HashMap::from([
            (0, HashSet::from([0])), (1, HashSet::from([1, 2, 3])), (2, HashSet::from([4])),
        ]);
        let flanks = HashMap::from([
            (0, (Some(10), None)), (1, (Some(10), Some(11))), (2, (Some(10), Some(12))),
        ]);
        assert!(target(&reps, &members, &flanks).is_none());
    }

    fn allele_fixture(
        short: &[u8], long: &[u8], side: AlleleSide, reverse: bool,
    ) -> (Pools, Vec<Call>, Vec<Vec<usize>>) {
        let mut pools = Pools::default();
        let mut calls = Vec::new();
        let (left, right) = if reverse { (11, 10) } else { (10, 11) };
        let missing_left = (side == AlleleSide::MissingPrefix && !reverse)
            || (side == AlleleSide::MissingSuffix && reverse);
        let missing_right = (side == AlleleSide::MissingSuffix && !reverse)
            || (side == AlleleSide::MissingPrefix && reverse);
        if !missing_left {
            calls.push(call(0, left, 100, b"MKTLLYVWAACD", &mut pools));
        }
        let mut middle = call(0, 0, if missing_left { 19 } else { 400 }, short, &mut pools);
        middle.strand = if reverse { -1 } else { 1 };
        let short_end = middle.end;
        calls.push(middle);
        if !missing_right {
            calls.push(call(0, right, short_end + 80, b"ACDEFGHIKLMNP", &mut pools));
        }
        for sample in 1..=2 {
            let reverse = sample == 2;
            let (left, right) = if reverse { (11, 10) } else { (10, 11) };
            calls.push(call(sample, left, 100, b"MKTLLYVWAACD", &mut pools));
            let mut middle = call(sample, 1, 400, long, &mut pools);
            middle.strand = if reverse { -1 } else { 1 };
            calls.push(middle);
            calls.push(call(sample, right, 1200, b"ACDEFGHIKLMNP", &mut pools));
        }
        let lengths = vec![vec![if missing_right { short_end as usize + 78 } else { 4000 }],
                           vec![4000], vec![4000]];
        (pools, calls, lengths)
    }

    fn allele_target(calls: &mut [Call], pools: &Pools, lengths: &[Vec<usize>]) -> AlleleDecision {
        sort_calls(calls);
        let loci = physical_loci(calls, pools);
        let contexts = locus_contexts(calls, &loci, |call| call.cluster);
        AlleleIndex::new(calls, pools, &loci, &contexts, lengths).target(0, pools)
    }

    fn intact_protein(length: usize, seed: u64) -> Vec<u8> {
        let mut full = protein(length, seed);
        full[0] = b'M';
        full
    }

    #[test]
    fn strand_normalized_terminal_rescue_handles_both_contig_ends_and_settles() {
        let full = intact_protein(120, 42);
        for reverse in [false, true] {
            for length in [60, 90, 102, 114] {
                for side in [AlleleSide::MissingPrefix, AlleleSide::MissingSuffix] {
                    let short = if side == AlleleSide::MissingPrefix {
                        &full[full.len() - length..]
                    } else {
                        &full[..length]
                    };
                    let (pools, mut calls, lengths) = allele_fixture(short, &full, side, reverse);
                    assert_eq!(allele_target(&mut calls, &pools, &lengths).target, Some(1));
                    freeze_contexts(&mut calls, &pools);
                    assert!(merge_fragments(&mut calls, &pools, &lengths));
                    let id = calls.iter().find(|call| call.sample == 0
                        && call.protein == pools.protein.iter().position(|sequence| sequence == short).unwrap() as u32)
                        .unwrap().cluster;
                    assert!(calls.iter().filter(|call| call.protein == pools.protein.iter()
                        .position(|sequence| sequence == &full).unwrap() as u32).all(|call| call.cluster == id));
                    assert_eq!(split_paralogs(&mut calls, &pools), 0);
                    assert!(!merge_fragments(&mut calls, &pools, &lengths));
                }
            }
        }
    }

    #[test]
    fn both_flanks_can_rescue_shortened_alleles_in_reversed_neighborhoods() {
        let full = intact_protein(120, 42);
        for reverse in [false, true] {
            for short in [&full[..90], &full[30..]] {
                let (pools, mut calls, lengths) = allele_fixture(short, &full, AlleleSide::Either, reverse);
                assert_eq!(allele_target(&mut calls, &pools, &lengths).target, Some(1));
            }
        }
    }

    #[test]
    fn one_flank_requires_a_real_nearby_contig_edge_on_the_matching_side() {
        let full = intact_protein(120, 42);
        let (pools, mut calls, mut lengths) =
            allele_fixture(&full[..90], &full, AlleleSide::MissingSuffix, false);
        lengths[0][0] += 13;
        assert!(allele_target(&mut calls, &pools, &lengths).target.is_none(), "91 bases is not near the edge");
        let (pools, mut calls, lengths) =
            allele_fixture(&full[..90], &full, AlleleSide::MissingPrefix, false);
        assert!(allele_target(&mut calls, &pools, &lengths).target.is_none(), "N-edge evidence cannot justify a C-terminal deletion");
        let (pools, mut calls, lengths) =
            allele_fixture(&full[30..], &full, AlleleSide::MissingPrefix, false);
        calls.iter_mut().find(|call| call.cluster == 0).unwrap().start = 91;
        assert!(allele_target(&mut calls, &pools, &lengths).target.is_none());
        let (mut pools, mut calls, lengths) =
            allele_fixture(&full[30..], &full, AlleleSide::MissingPrefix, false);
        calls.push(call(0, 99, 1, &intact_protein(80, 123), &mut pools));
        assert!(allele_target(&mut calls, &pools, &lengths).target.is_none(),
            "an overlapping predecessor is not a lost contig flank");
    }

    #[test]
    fn rescue_requires_two_intact_full_length_context_witnesses() {
        let full = intact_protein(120, 42);
        let (pools, mut calls, lengths) =
            allele_fixture(&full[30..], &full, AlleleSide::MissingPrefix, false);
        calls.retain(|call| call.sample != 2);
        assert!(allele_target(&mut calls, &pools, &lengths).target.is_none());
        let (mut pools, mut calls, lengths) =
            allele_fixture(&full[30..], &full, AlleleSide::MissingPrefix, false);
        let short_dna = crate::genetic_code::back_translate(&full[..60]);
        let shorter = pools.intern_protein(full[..60].to_vec());
        let dna = pools.intern_dna(short_dna);
        let witness = calls.iter_mut().find(|call| call.sample == 2 && call.cluster == 1).unwrap();
        witness.protein = shorter;
        witness.dna = dna;
        assert!(allele_target(&mut calls, &pools, &lengths).target.is_none(),
            "partial alleles cannot count as full-length witnesses");
        let mut broken = full.clone();
        broken[0] = b'A';
        let (pools, mut calls, lengths) =
            allele_fixture(&full[30..], &broken, AlleleSide::MissingPrefix, false);
        let decision = allele_target(&mut calls, &pools, &lengths);
        assert!(decision.target.is_none());
        assert!(decision.ambiguous);
        assert!(!merge_fragments(&mut calls, &pools, &lengths),
            "a nonintact host must not bypass rescue rejection through the legacy route");
    }

    #[test]
    fn a_one_off_flank_does_not_veto_adoption_when_the_other_flank_is_supported() {
        let full = intact_protein(120, 42);
        let short = &full[..90];
        let mut pools = Pools::default();
        let mut calls = Vec::new();
        let mut lengths = Vec::new();
        for sample in 0..2 {
            lengths.push(vec![4000]);
            bracketed(&mut calls, &mut pools, sample, 1, (10, 11), &[(0, short)]);
        }
        for sample in 2..4 {
            lengths.push(vec![4000]);
            bracketed(&mut calls, &mut pools, sample, 1, (10, 11), &[(1, full.as_slice())]);
        }
        lengths.push(vec![4000]);
        bracketed(&mut calls, &mut pools, 4, 1, (10, 12), &[(0, short)]);
        assert!(allele_target(&mut calls, &pools, &lengths).target.is_none());
        sort_calls(&mut calls);
        let loci = physical_loci(&calls, &pools);
        let raw = locus_contexts(&calls, &loci, |call| call.cluster);
        let stabilized = stabilize_contexts(&calls, &loci, &raw, |call| call.cluster);
        let decision = AlleleIndex::new(&calls, &pools, &loci, &stabilized, &lengths).target(0, &pools);
        assert_eq!(decision.target, Some(1));
        assert!(merge_fragments(&mut calls, &pools, &lengths));
        let short_protein = pools.protein.iter().position(|sequence| sequence == short).unwrap() as u32;
        let full_protein = pools.protein.iter().position(|sequence| sequence == &full).unwrap() as u32;
        let short_cluster = calls.iter().find(|call| call.protein == short_protein).unwrap().cluster;
        let full_cluster = calls.iter().find(|call| call.protein == full_protein).unwrap().cluster;
        assert_eq!(short_cluster, full_cluster);
    }

    #[test]
    fn length_evidence_uses_250_bp_and_the_collection() {
        use super::{length_evidence, CONSERVED_SHORT, LENGTH_SUFFICIENT, UNCONFIRMED_SHORT};
        assert_eq!(length_evidence(250, 1, 10, false), LENGTH_SUFFICIENT);
        assert_eq!(length_evidence(249, 1, 10, false), UNCONFIRMED_SHORT);
        assert_eq!(length_evidence(249, 1, 10, true), CONSERVED_SHORT);
        assert_eq!(length_evidence(249, 2, 10, false), CONSERVED_SHORT);
        assert_eq!(length_evidence(90, 1, 1, false), UNCONFIRMED_SHORT);
    }

    #[test]
    fn rescue_checks_every_short_allele_and_every_observed_context() {
        let full = intact_protein(120, 42);
        for change_context in [false, true] {
            let (mut pools, mut calls, mut lengths) =
                allele_fixture(&full[..90], &full, AlleleSide::Either, false);
            lengths.push(vec![4000]);
            let mut variant = full[..90].to_vec();
            if !change_context {
                for index in [20, 40] {
                    variant[index] = if variant[index] == b'A' { b'C' } else { b'A' };
                }
            }
            bracketed(&mut calls, &mut pools, 3, 1,
                if change_context { (10, 12) } else { (10, 11) }, &[(0, &variant)]);
            assert!(allele_target(&mut calls, &pools, &lengths).target.is_none());
        }
    }

    #[test]
    fn unresolved_multicopy_targets_are_not_allele_hosts() {
        let full = intact_protein(120, 42);
        let (mut pools, mut calls, lengths) =
            allele_fixture(&full[..90], &full, AlleleSide::Either, false);
        calls.push(call(1, 1, 2200, &full, &mut pools));
        assert!(allele_target(&mut calls, &pools, &lengths).target.is_none());
    }

    #[test]
    fn a_cooccurring_sequence_host_does_not_make_a_competing_host_look_unique() {
        let full = intact_protein(120, 42);
        let (mut pools, mut calls, mut lengths) =
            allele_fixture(&full[..90], &full, AlleleSide::Either, false);
        for sample in 3..=4 {
            lengths.push(vec![5000]);
            let model: Vec<_> = calls.iter().filter(|call| call.sample == 1).cloned().collect();
            for mut copy in model {
                copy.sample = sample;
                if copy.cluster == 1 { copy.cluster = 2; }
                calls.push(copy);
            }
        }
        bracketed(&mut calls, &mut pools, 3, 3001, (10, 11), &[(0, &full[..90])]);
        let decision = allele_target(&mut calls, &pools, &lengths);
        assert!(decision.ambiguous);
        assert!(decision.target.is_none());
        assert!(!merge_fragments(&mut calls, &pools, &lengths));
        assert_eq!(calls.iter().filter(|call| call.protein == pools.protein.iter()
            .position(|sequence| sequence == &full[..90]).unwrap() as u32)
            .map(|call| call.cluster).collect::<HashSet<_>>(), HashSet::from([0]));
        calls.retain(|call| call.sample != 4);
        assert!(allele_target(&mut calls, &pools, &lengths).ambiguous,
            "a singleton alternative host is still ambiguity, not permission to choose the frequent host");
    }

    #[test]
    fn damaged_and_bundled_hosts_still_count_as_ambiguity_evidence() {
        let full = intact_protein(120, 71);
        for bundled in [false, true] {
            let (mut pools, mut calls, mut lengths) =
                allele_fixture(&full[30..], &full, AlleleSide::MissingPrefix, false);
            let upstream = pools.protein[calls.iter().find(|call| call.cluster == 10).unwrap().protein as usize].clone();
            let downstream = pools.protein[calls.iter().find(|call| call.cluster == 11).unwrap().protein as usize].clone();
            for sample in 3..=4 {
                calls.extend([
                    call(sample, 10, 200, &upstream, &mut pools),
                    call(sample, 2, 400, &full, &mut pools),
                    call(sample, 11, 850, &downstream, &mut pools),
                ]);
                lengths.push(vec![1200]);
            }
            if bundled {
                calls.push(call(1, 1, 450, &full[15..65], &mut pools));
                calls.push(call(2, 1, 450, &full[15..65], &mut pools));
            } else {
                let mut damaged = back_translate(&full);
                damaged.replace_range(0..3, "GCT");
                assert!(!is_intact(&damaged));
                let protein = pools.intern_protein(translate(&damaged));
                let dna = pools.intern_dna(damaged);
                for call in &mut calls {
                    if call.cluster == 1 {
                        call.dna = dna;
                        call.protein = protein;
                    }
                }
            }
            let decision = allele_target(&mut calls, &pools, &lengths);
            assert!(decision.ambiguous);
            assert!(decision.target.is_none());
            assert!(!merge_fragments(&mut calls, &pools, &lengths));
        }
    }

    #[test]
    fn a_unique_cooccurring_host_is_blocked_during_union_application() {
        let full = intact_protein(120, 42);
        let (pools, mut calls, lengths) =
            allele_fixture(&full[..90], &full, AlleleSide::Either, false);
        let model: Vec<_> = calls.iter().filter(|call| call.sample == 1).cloned().collect();
        for mut copy in model {
            copy.sample = 0;
            copy.start += 2000;
            copy.end += 2000;
            calls.push(copy);
        }
        assert_eq!(allele_target(&mut calls, &pools, &lengths).target, Some(1));
        assert!(!merge_fragments(&mut calls, &pools, &lengths));
        let short = calls.iter().find(|call| call.protein == pools.protein.iter()
            .position(|sequence| sequence == &full[..90]).unwrap() as u32).unwrap();
        assert_eq!(short.cluster, 0);
        // The source group is not combined even though its only sequence host is clear.
        assert!(calls.iter().filter(|call| call.protein == pools.protein.iter()
            .position(|sequence| sequence == &full).unwrap() as u32).all(|call| call.cluster == 1));
    }

    #[test]
    fn shortened_matching_keeps_full_query_identity_and_rejects_domains_and_unknowns() {
        let full = intact_protein(120, 42);
        assert_eq!(shortened_identity(&full[..60], &full, AlleleSide::Either), Some(1.0));
        assert!(shortened_identity(&full[..59], &full, AlleleSide::Either).is_none());
        assert!(shortened_identity(&full, &full, AlleleSide::Either).is_none());
        assert!(shortened_identity(&full[..29], &full[..40], AlleleSide::Either).is_none());
        let mut unknown = full[..90].to_vec();
        unknown[10] = b'X';
        assert!(shortened_identity(&unknown, &full, AlleleSide::Either).is_none());
        unknown[10] = b'*';
        assert!(shortened_identity(&unknown, &full, AlleleSide::Either).is_none());
        let interior = &full[10..100];
        assert!(shortened_identity(interior, &full, AlleleSide::Either).unwrap() < TERMINAL_IDENTITY);
        let mut deleted = full[..45].to_vec();
        deleted.extend_from_slice(&full[75..]);
        assert!(shortened_identity(&deleted, &full, AlleleSide::Either).unwrap() < TERMINAL_IDENTITY);
        let mut variant = full[..100].to_vec();
        for index in [10, 30] {
            variant[index] = if variant[index] == b'A' { b'C' } else { b'A' };
        }
        assert_eq!(shortened_identity(&variant, &full, AlleleSide::Either), Some(0.98));
        variant[50] = if variant[50] == b'A' { b'C' } else { b'A' };
        assert!(shortened_identity(&variant, &full, AlleleSide::Either).unwrap() < TERMINAL_IDENTITY);
    }

    #[test]
    fn legacy_fragment_merges_keep_one_flank_and_containment_rules() {
        let full = protein(100, 42);
        let members = HashMap::from([(0, HashSet::from([0])), (1, HashSet::from([1]))]);
        let flanks = HashMap::from([
            (0, (Some(10), None)), (1, (Some(10), Some(11))),
        ]);
        for short in [full[..80].to_vec(), full[20..80].to_vec()] {
            let reps = HashMap::from([(0, short), (1, full.clone())]);
            let chosen = target(&reps, &members, &flanks).expect("legacy fragment");
            assert_eq!(chosen.cluster, 1);
            assert!(!chosen.terminal);
        }
    }

    fn call(sample: u32, cluster: u32, position: u32, sequence: &[u8], pools: &mut Pools) -> Call {
        let dna = crate::genetic_code::back_translate(sequence);
        Call {
            sample, contig: 0, start: position,
            end: position + dna.len() as u32 - 1, strand: 1,
            dna: pools.intern_dna(dna),
            protein: pools.intern_protein(sequence.to_vec()), cluster, context: cluster, locus: None,
            source: None, fusion_piece: false,
            symbol: String::new(), product: String::new(),
        }
    }

    #[test]
    fn transitive_terminal_merges_cannot_join_cooccurring_copies() {
        for legacy in [false, true] {
            let mut pools = Pools::default();
            let mut calls = vec![
                call(0, 0, 1, b"ACDEFG", &mut pools),
                call(1, 1, 1, b"ACDEFGH", &mut pools),
                call(0, 2, 100, b"ACDEFGHI", &mut pools),
            ];
            let members = sample_sets(&calls);
            let first = [(0, 1)];
            let all = [(0, 1), (1, 2)];
            let (merged, terminal, blocked) = if legacy {
                apply_adoptions(&mut calls, &pools, members, &first, &[(1, 2)])
            } else {
                apply_adoptions(&mut calls, &pools, members, &[], &all)
            };
            assert_eq!(merged, 1);
            assert_eq!(terminal, usize::from(!legacy));
            assert_eq!(blocked, 1);
            assert_eq!(calls[0].cluster, calls[1].cluster);
            assert_ne!(calls[0].cluster, calls[2].cluster);
        }
    }

    #[test]
    fn moderate_terminal_merge_is_applied_and_settles() {
        let mut pools = Pools::default();
        let full = intact_protein(100, 42);
        let left = protein(30, 123);
        let right = protein(31, 456);
        let mut calls = Vec::new();
        for (sample, cluster, middle) in [(0, 0, &full[..90]), (1, 1, full.as_slice())] {
            calls.push(call(sample, 10, 1, &left, &mut pools));
            calls.push(call(sample, cluster, 200, middle, &mut pools));
            calls.push(call(sample, 11, 600, &right, &mut pools));
        }
        let decision = allele_target(&mut calls, &pools, &[vec![1000], vec![1000]]);
        assert!(decision.target.is_none());
        assert!(!decision.ambiguous,
            "one intact witness must preserve eligibility for the original terminal rule");
        assert!(merge(&mut calls, &pools));
        let middle: Vec<_> = calls.iter().filter(|call| call.start == 200).collect();
        assert_eq!(middle[0].cluster, middle[1].cluster);
        assert!(!merge(&mut calls, &pools));
    }

    fn planned_partition(
        records: &[(u32, (Option<u32>, Option<u32>))],
    ) -> (Vec<Vec<usize>>, usize) {
        let mut pools = Pools::default();
        let calls: Vec<_> = records.iter().enumerate().map(|(index, &(sample, _))| {
            call(sample, 0, 1000 * index as u32 + 1, b"ACDEFGHIKL", &mut pools)
        }).collect();
        let loci: Vec<_> = calls.iter().enumerate().map(|(index, call)| {
            Locus { first: index, last: index, end: call.end }
        }).collect();
        let contexts: Vec<_> = records.iter().map(|&(_, context)| context).collect();
        let indexes: Vec<_> = (0..records.len()).collect();
        let plan = paralog_split_groups(&indexes, &calls, &loci, &contexts);
        let assigned: HashSet<_> = plan.groups.iter().flatten().copied().collect();
        let retained = indexes.into_iter().filter(|index| !assigned.contains(index)).collect();
        let mut groups = plan.groups;
        groups.push(retained);
        groups.sort();
        (groups, plan.unresolved_loci)
    }

    #[test]
    fn flank_variants_of_two_paralogs_produce_two_groups_not_four() {
        let records = [
            (0, (Some(10), Some(11))), (0, (Some(20), Some(21))),
            (1, (Some(10), Some(12))), (1, (Some(20), Some(22))),
            (2, (Some(10), None)),
        ];
        let (groups, unresolved) = planned_partition(&records);
        assert_eq!(groups, vec![vec![0, 2, 4], vec![1, 3]]);
        assert_eq!(unresolved, 0);
    }

    #[test]
    fn ambiguous_single_flank_is_not_assigned_to_an_arbitrary_paralog() {
        let records = [
            (0, (Some(10), Some(11))), (0, (Some(10), Some(12))),
            (1, (Some(10), None)), (2, (None, None)), (3, (None, None)),
        ];
        let (groups, unresolved) = planned_partition(&records);
        assert_eq!(groups, vec![vec![0], vec![1], vec![2, 3, 4]]);
        assert_eq!(unresolved, 3);
    }

    #[test]
    fn context_bridge_cannot_combine_cooccurring_paralogs() {
        let records = [
            (0, (Some(10), Some(11))), (1, (Some(10), Some(12))),
            (0, (Some(20), Some(12))),
        ];
        let (groups, _) = planned_partition(&records);
        assert!(!groups.iter().any(|group| group.contains(&0) && group.contains(&2)));
    }

    #[test]
    fn unrelated_single_copy_contexts_are_not_split() {
        let records = [
            (0, (Some(10), Some(11))), (1, (Some(20), Some(21))),
            (2, (None, None)),
        ];
        assert_eq!(planned_partition(&records), (vec![vec![0, 1, 2]], 0));
    }

    #[test]
    fn identical_multicopy_contexts_are_reported_as_unresolved() {
        let records = [
            (0, (Some(10), Some(11))), (0, (Some(10), Some(11))),
        ];
        assert_eq!(planned_partition(&records), (vec![vec![0, 1]], 1));
    }

    #[test]
    fn fragment_loci_have_external_flanks_and_one_vote() {
        let mut pools = Pools::default();
        let full = protein(80, 42);
        let mut calls = vec![
            call(0, 10, 1, b"MKTLLYVWAACD", &mut pools),
            call(0, 0, 100, &full[..40], &mut pools),
        ];
        let start = calls[1].end + 1;
        calls.push(call(0, 0, start, &full[40..], &mut pools));
        calls.push(call(0, 11, 600, b"ACDEFGHIKLMNP", &mut pools));
        let loci = physical_loci(&calls, &pools);
        assert_eq!(loci.len(), 3);
        assert_eq!((loci[1].first, loci[1].last), (1, 2));
        assert_eq!(locus_contexts(&calls, &loci, |call| call.cluster)[1], (Some(10), Some(11)));
        assert_eq!(modal_flanks(&calls, &pools)[&0], (Some(10), Some(11)));
        assert_eq!(split_paralogs(&mut calls, &pools), 0);
    }

    #[test]
    fn touching_full_copies_and_opposite_strands_are_distinct_loci() {
        let mut pools = Pools::default();
        let full = protein(80, 42);
        let first = call(0, 0, 1, &full, &mut pools);
        let second = call(0, 0, first.end + 1, &full, &mut pools);
        assert_eq!(physical_loci(&[first, second], &pools).len(), 2);
        let first = call(0, 0, 1, &full[..40], &mut pools);
        let mut second = call(0, 0, first.end + 1, &full[40..], &mut pools);
        second.strand = -1;
        assert_eq!(physical_loci(&[first, second], &pools).len(), 2);
    }

    #[test]
    fn strongly_overlapping_annotations_are_one_locus() {
        let mut pools = Pools::default();
        let full = protein(100, 42);
        let calls = vec![
            call(0, 0, 1, &full, &mut pools),
            call(0, 0, 1, &full[..90], &mut pools),
        ];
        assert_eq!(physical_loci(&calls, &pools).len(), 1);
    }

    #[test]
    fn repetitive_fragments_use_a_longer_observed_template_not_a_tandem_copy() {
        let mut pools = Pools::default();
        let part = protein(40, 42);
        let mut full = part.clone();
        full.extend_from_slice(&part);
        let first = call(0, 0, 1, &part, &mut pools);
        let second = call(0, 0, first.end + 1, &part, &mut pools);
        let observed = call(1, 0, 1, &full, &mut pools);
        let mut calls = vec![first, second, observed];
        assert_eq!(physical_loci(&calls, &pools).len(), 2);
        assert_eq!(split_paralogs(&mut calls, &pools), 0);
        let first = call(0, 0, 1, &full, &mut pools);
        let second = call(0, 0, first.end + 1, &full, &mut pools);
        assert_eq!(physical_loci(&[first, second], &pools).len(), 2);
    }

    #[test]
    fn repetitive_fragment_union_remains_one_locus_after_assignment() {
        let mut pools = Pools::default();
        let part = protein(40, 42);
        let mut full = part.clone();
        full.extend_from_slice(&part);
        let first = call(0, 0, 1, &part, &mut pools);
        let second = call(0, 1, first.end + 1, &part, &mut pools);
        let observed = call(1, 2, 1, &full, &mut pools);
        let mut calls = vec![first, second, observed];
        let members = sample_sets(&calls);
        assert_eq!(
            apply_adoptions(&mut calls, &pools, members, &[(0, 2), (1, 2)], &[]),
            (2, 0, 0),
        );
        assert_eq!(physical_loci(&calls, &pools).len(), 2);
        assert_eq!(split_paralogs(&mut calls, &pools), 0);
    }

    #[test]
    fn legacy_merges_cannot_recombine_distant_or_other_contig_copies() {
        for other_contig in [false, true] {
            let mut pools = Pools::default();
            let full = protein(100, 42);
            let mut calls = vec![
                call(0, 0, 1, &full[..40], &mut pools),
                call(0, 1, 1000, &full, &mut pools),
            ];
            if other_contig {
                calls[1].contig = 1;
            }
            freeze_contexts(&mut calls, &pools);
            let members = sample_sets(&calls);
            assert_eq!(
                apply_adoptions(&mut calls, &pools, members, &[(0, 1)], &[]),
                (0, 0, 1),
            );
            assert_ne!(calls[0].cluster, calls[1].cluster);
            assert_ne!(calls[0].context, calls[1].context);
            assert_ne!(calls[0].locus, calls[1].locus);
        }
    }

    #[test]
    fn legacy_unions_check_accumulated_sample_conflicts() {
        let mut pools = Pools::default();
        let full = protein(100, 42);
        let mut calls = vec![
            call(0, 0, 1, &full[..40], &mut pools),
            call(0, 2, 1000, &full, &mut pools),
            call(1, 1, 1, &full[..60], &mut pools),
        ];
        let members = sample_sets(&calls);
        assert_eq!(
            apply_adoptions(&mut calls, &pools, members, &[(0, 1), (1, 2)], &[]),
            (1, 0, 1),
        );
        assert_eq!(calls[0].cluster, calls[2].cluster);
        assert_ne!(calls[0].cluster, calls[1].cluster);
    }

    fn bracketed(
        calls: &mut Vec<Call>, pools: &mut Pools, sample: u32, position: u32,
        flanks: (u32, u32), parts: &[(u32, &[u8])],
    ) {
        calls.push(call(sample, flanks.0, position, b"MKTLLYVWAACD", pools));
        let mut start = position + 100;
        for &(cluster, sequence) in parts {
            let part = call(sample, cluster, start, sequence, pools);
            start = part.end + 1;
            calls.push(part);
        }
        calls.push(call(sample, flanks.1, position + 600, b"ACDEFGHIKLMNP", pools));
    }

    #[test]
    fn splitting_neighbors_does_not_invent_new_paralog_contexts() {
        let mut pools = Pools::default();
        let full = protein(80, 42);
        let parts: &[(u32, &[u8])] = &[(0, &full)];
        let mut calls = Vec::new();
        for sample in 0..2 {
            bracketed(&mut calls, &mut pools, sample, 1, (10, 11), parts);
            bracketed(&mut calls, &mut pools, sample, 1001, (10, 11), parts);
        }
        freeze_contexts(&mut calls, &pools);
        assert_eq!(split_paralogs(&mut calls, &pools), 2, "only the flanks split");
        let loci = physical_loci(&calls, &pools);
        let current = locus_contexts(&calls, &loci, |call| call.cluster);
        let stable = locus_contexts(&calls, &loci, |call| call.context);
        let core_contexts = |contexts: &[(Option<u32>, Option<u32>)]| {
            loci.iter().zip(contexts).filter_map(|(locus, &context)| {
                (pools.protein[calls[locus.first].protein as usize] == full).then_some(context)
            }).collect::<HashSet<_>>()
        };
        assert_eq!(core_contexts(&current).len(), 2);
        assert_eq!(core_contexts(&stable).len(), 1);
        for _ in 0..4 {
            assert_eq!(split_paralogs(&mut calls, &pools), 0);
            compress_clusters(&mut calls);
        }
        let groups: HashSet<_> = calls.iter().filter_map(|call| {
            (pools.protein[call.protein as usize] == full).then_some(call.cluster)
        }).collect();
        assert_eq!(groups.len(), 1, "identical contexts must remain explicitly unresolved");
    }

    #[test]
    fn repetitive_fragment_proof_survives_removing_the_full_paralog() {
        let mut pools = Pools::default();
        let half = protein(40, 42);
        let full = [half.as_slice(), half.as_slice()].concat();
        let fragments: &[(u32, &[u8])] = &[(0, &half), (0, &half)];
        let complete: &[(u32, &[u8])] = &[(0, &full)];
        let mut calls = Vec::new();
        for sample in 0..2 {
            bracketed(&mut calls, &mut pools, sample, 1, (10, 11), fragments);
            bracketed(&mut calls, &mut pools, sample, 1001, (20, 21), complete);
        }
        freeze_contexts(&mut calls, &pools);
        assert_eq!(split_paralogs(&mut calls, &pools), 1);
        for _ in 0..4 {
            let loci = physical_loci(&calls, &pools);
            assert_eq!(loci.iter().filter(|locus| {
                pools.protein[calls[locus.first].protein as usize].len() >= 40
            }).count(), 4, "each fragmented copy is still one locus");
            assert_eq!(split_paralogs(&mut calls, &pools), 0);
        }
        let a = calls.iter().find(|call| call.sample == 0 && call.start == 101).unwrap();
        let b = calls.iter().find(|call| call.sample == 0 && call.start == 1101).unwrap();
        assert_ne!(a.cluster, b.cluster);
    }

    #[test]
    fn proven_fragment_unions_coarsen_contexts_across_existing_descendants() {
        let mut pools = Pools::default();
        let full = protein(100, 42);
        let mut calls = vec![
            call(0, 10, 1, &full[..90], &mut pools),
            call(1, 11, 1, &full, &mut pools),
            call(2, 10, 1, &full[..90], &mut pools),
        ];
        freeze_contexts(&mut calls, &pools);
        calls[2].cluster = 12;
        let members = sample_sets(&calls);
        assert_eq!(apply_adoptions(&mut calls, &pools, members, &[], &[(10, 11)]), (1, 1, 0));
        assert_eq!(calls[0].cluster, calls[1].cluster);
        assert_ne!(calls[0].cluster, calls[2].cluster, "anchors are not ortholog assignments");
        assert!(calls.iter().all(|call| call.context == 10));
        compress_clusters(&mut calls);
        assert!(calls.iter().all(|call| call.context == 10));
    }

    #[test]
    fn exact_whole_alleles_bridge_both_changed_flanks_without_combining_copies() {
        let mut pools = Pools::default();
        let a = protein(80, 42);
        let b = protein(80, 123);
        let mut calls = Vec::new();
        for (sample, offset) in [(0, 0), (1, 20)] {
            bracketed(&mut calls, &mut pools, sample, 1, (10 + offset, 11 + offset), &[(0, &a)]);
            bracketed(&mut calls, &mut pools, sample, 1001, (12 + offset, 13 + offset), &[(0, &b)]);
        }
        freeze_contexts(&mut calls, &pools);
        assert_eq!(split_paralogs(&mut calls, &pools), 1);
        let cluster_at = |sample, position| {
            calls.iter().find(|call| call.sample == sample && call.start == position).unwrap().cluster
        };
        assert_eq!(cluster_at(0, 101), cluster_at(1, 101));
        assert_eq!(cluster_at(0, 1101), cluster_at(1, 1101));
        assert_ne!(cluster_at(0, 101), cluster_at(0, 1101));
        assert_eq!(split_paralogs(&mut calls, &pools), 0);
    }

    #[test]
    fn identical_allele_rescue_does_not_choose_between_cooccurring_hosts() {
        let mut pools = Pools::default();
        let full = protein(80, 42);
        let mut calls = Vec::new();
        for (sample, offset) in [(0, 0), (1, 20)] {
            bracketed(&mut calls, &mut pools, sample, 1, (10 + offset, 11 + offset), &[(0, &full)]);
            bracketed(&mut calls, &mut pools, sample, 1001, (12 + offset, 13 + offset), &[(0, &full)]);
        }
        freeze_contexts(&mut calls, &pools);
        assert_eq!(split_paralogs(&mut calls, &pools), 3);
        let groups: HashSet<_> = calls.iter().filter_map(|call| {
            (pools.protein[call.protein as usize] == full).then_some(call.cluster)
        }).collect();
        assert_eq!(groups.len(), 4, "ambiguous cross-sample copies must not be assigned by order");
        assert_eq!(split_paralogs(&mut calls, &pools), 0);
    }

    #[test]
    fn allele_rescue_requires_exact_proteins_not_near_identical_ones() {
        let mut pools = Pools::default();
        let a = protein(80, 42);
        let mut variant = a.clone();
        variant[10] = if variant[10] == b'A' { b'C' } else { b'A' };
        let b = protein(80, 123);
        let mut calls = Vec::new();
        for (sample, offset, allele) in [(0, 0, a.as_slice()), (1, 20, variant.as_slice())] {
            bracketed(&mut calls, &mut pools, sample, 1, (10 + offset, 11 + offset), &[(0, allele)]);
            bracketed(&mut calls, &mut pools, sample, 1001, (12 + offset, 13 + offset), &[(0, &b)]);
        }
        freeze_contexts(&mut calls, &pools);
        assert_eq!(split_paralogs(&mut calls, &pools), 2);
        let cluster_at = |sample, position| {
            calls.iter().find(|call| call.sample == sample && call.start == position).unwrap().cluster
        };
        assert_ne!(cluster_at(0, 101), cluster_at(1, 101));
        assert_eq!(cluster_at(0, 1101), cluster_at(1, 1101));
    }

    #[test]
    fn a_fragment_part_is_not_whole_allele_rescue_evidence() {
        let mut pools = Pools::default();
        let head = b"AAAAAAAAAAAAAAAAAAAA";
        let tail = b"CCCCCCCCCCCCCCCCCCCC";
        let full = protein(40, 42);
        let mut calls = Vec::new();
        bracketed(&mut calls, &mut pools, 0, 1, (10, 11), &[(0, head), (0, tail)]);
        bracketed(&mut calls, &mut pools, 0, 1001, (20, 21), &[(0, &full)]);
        bracketed(&mut calls, &mut pools, 1, 1, (30, 31), &[(0, head)]);
        freeze_contexts(&mut calls, &pools);
        assert_eq!(split_paralogs(&mut calls, &pools), 2);
        let a = calls.iter().find(|call| call.sample == 0 && call.start == 101).unwrap();
        let b = calls.iter().find(|call| call.sample == 1 && call.start == 101).unwrap();
        assert_ne!(a.cluster, b.cluster);
    }

    #[test]
    fn context_and_allele_partitioning_is_idempotent_and_preserves_sample_constraints() {
        for seed in 1u64..=64 {
            let mut state = seed;
            let mut next = || {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                state
            };
            let mut pools = Pools::default();
            let mut contexts = Vec::new();
            let calls: Vec<_> = (0..30).map(|index| {
                let sample = (next() % 6) as u32;
                let left = (next() % 5) as u32;
                let right = (next() % 5) as u32;
                contexts.push(((left != 0).then_some(left), (right != 0).then_some(right)));
                let allele = if next() % 2 == 0 { b"ACDEFGHIKL" } else { b"MNPQRSTVWY" };
                call(sample, 0, index * 1000 + 1, allele, &mut pools)
            }).collect();
            let loci: Vec<_> = calls.iter().enumerate().map(|(index, call)| {
                Locus { first: index, last: index, end: call.end }
            }).collect();
            let indexes: Vec<_> = (0..calls.len()).collect();
            let plan = paralog_split_groups(&indexes, &calls, &loci, &contexts);
            let assigned: HashSet<_> = plan.groups.iter().flatten().copied().collect();
            let retained: Vec<_> = indexes.into_iter().filter(|index| !assigned.contains(index)).collect();
            let groups = plan.groups.into_iter().chain(std::iter::once(retained));
            for group in groups {
                assert!(paralog_split_groups(&group, &calls, &loci, &contexts).groups.is_empty(), "seed {seed}");
                for (position, &a) in group.iter().enumerate() {
                    for &b in &group[position + 1..] {
                        if calls[a].sample == calls[b].sample {
                            assert_eq!(contexts[a], contexts[b], "co-occurring contexts joined for seed {seed}");
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn fragments_stay_merged_across_paralog_rounds_without_collapsing_copies() {
        let mut pools = Pools::default();
        let full = protein(80, 42);
        let fragments: &[(u32, &[u8])] = &[(0, &full[..40]), (0, &full[40..])];
        let unmerged: &[(u32, &[u8])] = &[(1, &full[..40]), (2, &full[40..])];
        let complete: &[(u32, &[u8])] = &[(0, &full)];
        let mut calls = Vec::new();
        bracketed(&mut calls, &mut pools, 0, 1, (10, 11), fragments);
        bracketed(&mut calls, &mut pools, 0, 1001, (20, 21), complete);
        bracketed(&mut calls, &mut pools, 1, 1, (10, 11), complete);
        bracketed(&mut calls, &mut pools, 1, 1001, (20, 21), complete);
        bracketed(&mut calls, &mut pools, 2, 1, (10, 11), unmerged);
        bracketed(&mut calls, &mut pools, 2, 1001, (20, 21), complete);
        freeze_contexts(&mut calls, &pools);
        assert_eq!(split_paralogs(&mut calls, &pools), 1);
        assert!(merge(&mut calls, &pools));
        compress_clusters(&mut calls);
        let at_a = calls.iter().find(|call| call.sample == 1 && call.start == 101).unwrap().cluster;
        let at_b = calls.iter().find(|call| call.sample == 1 && call.start == 1101).unwrap().cluster;
        assert_ne!(at_a, at_b);
        for call in calls.iter().filter(|call| call.sample == 2 && (101..500).contains(&call.start)) {
            assert_eq!(call.cluster, at_a);
        }
        for _ in 0..3 {
            assert_eq!(split_paralogs(&mut calls, &pools), 0);
            assert!(!merge(&mut calls, &pools));
        }
    }

    #[test]
    fn paralog_assignments_do_not_depend_on_input_iteration_order() {
        let mut pools = Pools::default();
        let full = protein(80, 42);
        let parts: &[(u32, &[u8])] = &[(0, &full)];
        let mut calls = Vec::new();
        for (sample, right_a, right_b) in [(0, 11, 21), (1, 12, 22)] {
            bracketed(&mut calls, &mut pools, sample, 1, (10, right_a), parts);
            bracketed(&mut calls, &mut pools, sample, 1001, (20, right_b), parts);
        }
        let mut reversed = calls.clone();
        reversed.reverse();
        assert_eq!(split_paralogs(&mut calls, &pools), 1);
        assert_eq!(split_paralogs(&mut reversed, &pools), 1);
        let labels = |calls: &[Call]| {
            calls.iter().map(|call| (call.sample, call.start, call.cluster)).collect::<Vec<_>>()
        };
        assert_eq!(labels(&calls), labels(&reversed));
    }

    #[test]
    fn union_find_resolves_and_compresses_chains_longer_than_64() {
        let mut parent: HashMap<u32, u32> = (0..200).map(|id| (id, id + 1)).collect();
        assert_eq!(find(&mut parent, 0), 200);
        assert_eq!(parent[&0], 200);
        assert_eq!(find(&mut parent, 1), 200);
    }
}
