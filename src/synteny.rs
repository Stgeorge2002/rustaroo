//! Ortholog groups are homology plus gene order.
//!
//! Sequence identity alone keeps a fragment out of the full gene and keeps a
//! fused CDS glued to whichever piece is longer. Gene order repairs both.
//! Paralogs that landed in one cluster are split by their flanking genes.
//! Fusions are given back to the two genes they contain. Fragments are adopted
//! by the full gene they are a piece of.

use crate::fxhash::{HashMap, HashSet};

use rayon::prelude::*;

use crate::align::{
    fit_short_in_long, fit_short_in_long_indexed, global_identity, terminal_identity_parts, LongIndex,
};
use crate::packed::PackedSeq;
use crate::genetic_code::{is_intact, is_start, is_stop, orf_status, translate, translate_bytes};
use crate::log::{self, progress};
use crate::model::{Call, Label, Pools};

const CONTAIN_IDENTITY: f64 = 0.90;
const CONTAIN_COVERAGE: f64 = 0.85;
const MAX_COOCCUR: f64 = 0.25;
const NEIGHBOR_GAP: i32 = 15_000;
const NEIGHBOR_OVERLAP: i64 = 30;
const ABUT_GAP: i32 = 3;
const TERMINAL_IDENTITY: f64 = 0.98;
pub(crate) const ALLELE_EDGE_BP: usize = 90;
const ALLELE_MIN_AA: usize = 30;

/// Split fusion CDS into the two genes they contain. Returns the number split.
pub fn split_fusions(
    calls: &mut Vec<Call>,
    pools: &mut Pools,
    contigs: &[Vec<PackedSeq>],
) -> usize {
    let (reps, counts) = rayon::join(|| representatives(calls, pools), || sample_counts(calls));
    let mut cluster_ids: Vec<u32> = reps.keys().copied().collect();
    // Sorted so the `best` tie-breaks below (strict `>` on a support score)
    // resolve the same way on every run, instead of depending on this
    // HashMap's randomized iteration order.
    cluster_ids.sort_unstable();
    // Finding the best (left, right) split for one `long_id` only reads
    // `reps`/`counts`/`kmers` - nothing here is mutated until the groups are
    // actually applied below - so every family's search is independent of
    // every other family's and can run on its own thread.
    // PERF PROTOTYPE: raw 4-window postings, ids ascending. Candidate lists
    // drawn from it are exactly the ids `shares_kmer_cached` / `fit_short_in_long`
    // could accept, visited in the same ascending order.
    let window_index = FusionIndex::new(&reps, &cluster_ids);
    let compositions: HashMap<u32, Composition> =
        reps.iter().map(|(&id, protein)| (id, composition(protein))).collect();
    // Each representative's 4-mer table, built once instead of once per fit.
    let long_indexes: HashMap<u32, LongIndex> = cluster_ids
        .par_iter()
        .map(|&id| (id, LongIndex::new(&reps[&id])))
        .collect();
    let planned: Vec<(u32, u32, u32, usize)> = cluster_ids
        .par_iter()
        .filter_map(|&long_id| {
            best_fusion_split_indexed(
                long_id,
                &reps,
                &counts,
                &window_index,
                &compositions,
                &long_indexes,
            )
                .map(|(left_id, right_id, end)| (long_id, left_id, right_id, end))
        })
        .collect();
    let plan_of: HashMap<u32, (u32, u32, u32, usize)> =
        planned.iter().map(|item| (item.0, *item)).collect();
    if planned.is_empty() {
        return 0;
    }
    let mut twins = Vec::new();
    let mut split = 0usize;
    let mut not_contiguous = 0usize;
    for call in calls.iter_mut() {
        let Some((_, left_id, right_id, end)) = plan_of.get(&call.cluster)
        else {
            continue;
        };
        let protein_len = pools.protein[call.protein as usize].len();
        let cut = (*end).min(protein_len);
        if cut >= protein_len {
            // Nothing is left for the right-hand gene: the whole call is the
            // left gene, unchanged.
            call.cluster = *left_id;
            continue;
        }
        let dna = pools.dna[call.dna as usize].clone();
        // FIX(fusion): the pieces are the genome's own bases. No stop codon is
        // invented for the prefix, the suffix keeps the real terminal codon,
        // and both intervals are derived from the CDS's 5' end so they tile
        // the original span exactly (no overlap, no off-by-one).
        let Some(pieces) = fusion_pieces(call, &dna, cut, contigs) else {
            // A CDS whose bases are not one contiguous genome interval (joined
            // parts, phase offset, trimmed tail) cannot be cut into exact
            // intervals. It is left whole rather than given made-up coordinates.
            not_contiguous += 1;
            continue;
        };
        call.fusion_piece = true;
        let mut twin = call.clone();
        twin.cluster = *right_id;
        twin.symbol = Label::default();
        twin.product = Label::default();
        twin.protein = pools.intern_protein(translate(&pieces.suffix));
        twin.dna = pools.intern_dna(pieces.suffix);
        (twin.start, twin.end) = pieces.suffix_span;
        twins.push(twin);
        call.protein = pools.intern_protein(translate(&pieces.prefix));
        call.dna = pools.intern_dna(pieces.prefix);
        (call.start, call.end) = pieces.prefix_span;
        call.cluster = *left_id;
        split += 1;
    }
    if not_contiguous > 0 {
        progress!(
            "rustroo: left {not_contiguous} fused annotations whole: their CDS is not one contiguous genome interval"
        );
    }
    calls.extend(twins);
    split
}

/// The best `(left_id, right_id, cut_point)` fusion split for `long_id`, if
/// any: `left` has to explain the start of `long`, and the leftover suffix
/// has to be explained by `right`. Pulled out of `split_fusions` so the
/// search for each `long_id` can run as an independent unit of parallel work.
#[allow(dead_code)]
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

pub(crate) struct FusionIndex {
    postings: HashMap<[u8; 4], Vec<u32>>,
    /// Clusters whose representative is shorter than 4 residues: the old
    /// `shares_kmer_cached` treated those as always sharing.
    tiny: Vec<u32>,
}

impl FusionIndex {
    pub(crate) fn new(reps: &HashMap<u32, Vec<u8>>, cluster_ids: &[u32]) -> Self {
        let mut postings: HashMap<[u8; 4], Vec<u32>> = HashMap::default();
        let mut tiny = Vec::new();
        for &id in cluster_ids {
            let protein = &reps[&id];
            if protein.len() < 4 {
                tiny.push(id);
                continue;
            }
            for w in protein.windows(4) {
                let list = postings.entry([w[0], w[1], w[2], w[3]]).or_default();
                if list.last() != Some(&id) {
                    list.push(id);
                }
            }
        }
        Self { postings, tiny }
    }

    /// Ascending ids sharing at least one raw 4-window with `sequence`.
    pub(crate) fn sharing(&self, sequence: &[u8], include_tiny: bool) -> Vec<u32> {
        let mut out = Vec::new();
        if sequence.len() >= 4 {
            for w in sequence.windows(4) {
                if let Some(list) = self.postings.get(&[w[0], w[1], w[2], w[3]]) {
                    out.extend_from_slice(list);
                }
            }
        }
        if include_tiny {
            out.extend_from_slice(&self.tiny);
        }
        out.sort_unstable();
        out.dedup();
        out
    }
}

/// PERF PROTOTYPE of `best_fusion_split`: identical result. Lefts that share
/// no 4-window with `long` were skipped by `shares_kmer_cached`; rights that
/// share no 4-window with the suffix cannot get a `fit_short_in_long` result.
fn best_fusion_split_indexed(
    long_id: u32,
    reps: &HashMap<u32, Vec<u8>>,
    counts: &HashMap<u32, u32>,
    index: &FusionIndex,
    compositions: &HashMap<u32, Composition>,
    long_indexes: &HashMap<u32, LongIndex>,
) -> Option<(u32, u32, usize)> {
    let long = reps.get(&long_id)?;
    let long_index = &long_indexes[&long_id];
    let long_composition = composition(long);
    let long_samples = counts.get(&long_id).copied().unwrap_or(0);
    let mut best: Option<(u32, u32, usize, u32)> = None;
    let lefts = if long.len() < 4 {
        // shares_kmer_cached returns true for every left.
        let mut all: Vec<u32> = reps.keys().copied().collect();
        all.sort_unstable();
        all
    } else {
        index.sharing(long, true)
    };
    let mut suffix_cache: HashMap<usize, (Vec<u32>, Composition, LongIndex)> = HashMap::default();
    for left_id in lefts {
        if left_id == long_id || counts.get(&left_id).copied().unwrap_or(0) <= long_samples {
            continue;
        }
        let Some(left) = reps.get(&left_id) else {
            continue;
        };
        // Exact skip: a fit's identity is diagonal matches / left.len(), and
        // matches can never exceed the shared residue composition.
        if left.len() <= long.len()
            && !left.is_empty()
            && (composition_bound(&compositions[&left_id], &long_composition) as f64
                / left.len() as f64)
                < CONTAIN_IDENTITY
        {
            continue;
        }
        let Some(fit) = fit_short_in_long_indexed(left, long, long_index) else {
            continue;
        };
        if fit.identity < CONTAIN_IDENTITY
            || fit.coverage < CONTAIN_COVERAGE
            || fit.long_start > 2
            || fit.internal_gap > 5
        {
            continue;
        }
        let cut = fit.long_end.min(long.len());
        let suffix = &long[cut..];
        let (rights, suffix_composition, suffix_index) = suffix_cache.entry(cut).or_insert_with(|| {
            (
                if suffix.len() >= 4 { index.sharing(suffix, false) } else { Vec::new() },
                composition(suffix),
                LongIndex::new(suffix),
            )
        });
        for &right_id in rights.iter() {
            if right_id == long_id || right_id == left_id {
                continue;
            }
            if counts.get(&right_id).copied().unwrap_or(0) <= long_samples {
                continue;
            }
            let Some(right) = reps.get(&right_id) else {
                continue;
            };
            // Same exact composition skip for the right piece (identity is
            // over the shorter of right / suffix).
            let shorter = right.len().min(suffix.len());
            if shorter > 0
                && (composition_bound(&compositions[&right_id], suffix_composition) as f64
                    / shorter as f64)
                    < 0.85
            {
                continue;
            }
            let matched = if right.len() <= suffix.len() {
                fit_short_in_long_indexed(right, suffix, suffix_index)
            } else if suffix.len() >= 4 {
                fit_short_in_long_indexed(suffix, right, &long_indexes[&right_id])
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

/// Residue counts bucketed by `byte & 31`. Merging residues into one bucket
/// can only raise `sum min(count)`, so the bound stays a true upper bound on
/// identical aligned pairs.
pub(crate) type Composition = [u32; 32];

pub(crate) fn composition(sequence: &[u8]) -> Composition {
    let mut counts = [0u32; 32];
    for &residue in sequence {
        counts[(residue & 31) as usize] += 1;
    }
    counts
}

pub(crate) fn composition_bound(a: &Composition, b: &Composition) -> usize {
    a.iter().zip(b).map(|(x, y)| (*x).min(*y) as usize).sum()
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
    let mut by_cluster: HashMap<u32, Vec<usize>> = HashMap::default();
    for (index, locus) in loci.iter().enumerate() {
        by_cluster
            .entry(calls[locus.first].cluster)
            .or_default()
            .push(index);
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
    let mut per_sample: HashMap<u32, usize> = HashMap::default();
    for &index in indexes {
        *per_sample
            .entry(calls[loci[index].first].sample)
            .or_default() += 1;
    }
    if per_sample.values().all(|&count| count <= 1) {
        return ParalogPlan::default();
    }
    let mut by_context: HashMap<(Option<u32>, Option<u32>), Vec<usize>> = HashMap::default();
    for &index in indexes {
        by_context.entry(contexts[index]).or_default().push(index);
    }
    let mut ordered: Vec<_> = by_context.into_iter().collect();
    ordered.sort_by_key(|&(context, _)| context);
    let mut unresolved_loci = 0;
    let buckets: Vec<ContextBucket> = ordered
        .into_iter()
        .map(|(context, indexes)| {
            let mut samples = HashSet::default();
            for &index in &indexes {
                if !samples.insert(calls[loci[index].first].sample) {
                    // Equal contexts cannot resolve two copies in one sample.
                    unresolved_loci += 1;
                }
            }
            let anchored = samples.iter().any(|sample| per_sample[sample] > 1);
            ContextBucket {
                context,
                loci: indexes,
                samples,
                anchored,
            }
        })
        .collect();
    let mut left: HashMap<u32, Vec<usize>> = HashMap::default();
    let mut right: HashMap<u32, Vec<usize>> = HashMap::default();
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
    let mut parents: HashMap<u32, u32> = HashMap::default();
    let mut samples: Vec<HashSet<u32>> = buckets
        .iter()
        .map(|bucket| bucket.samples.clone())
        .collect();
    reconcile_disjoint_pairs(pairs, &mut parents, &mut samples);
    // Exact whole alleles can bridge rearranged contexts. Use roots after
    // flank reconciliation so weaker sequence evidence cannot undo it.
    let mut by_protein: HashMap<u32, Vec<usize>> = HashMap::default();
    for (id, bucket) in buckets.iter().enumerate() {
        let root = find(&mut parents, id as u32) as usize;
        for &index in &bucket.loci {
            let locus = &loci[index];
            if locus.first == locus.last {
                by_protein
                    .entry(calls[locus.first].protein)
                    .or_default()
                    .push(root);
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
    let mut grouped: HashMap<u32, (Vec<usize>, bool)> = HashMap::default();
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
    ParalogPlan {
        groups: groups.into_iter().skip(1).collect(),
        unresolved_loci,
        allele_unions,
    }
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
            let mut seen = HashSet::default();
            hosts
                .iter()
                .any(|&host| samples[host].iter().any(|&sample| !seen.insert(sample)))
        })
        .collect();
    pairs.sort_by(|&(a, b), &(c, d)| {
        let support = |left: usize, right: usize| samples[left].len() + samples[right].len();
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
        || {
            rayon::join(
                || modal_flanks_of(calls, &loci, &contexts),
                || kmer_sets(&reps),
            )
        },
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
    legacy.extend(abutting_pieces(calls, pools, &reps, &members));
    terminal.extend(alleles.same_locus_adoptions(calls, pools, &members));
    let contained = contained_adoptions(calls, pools, &reps, &members, &kmers);
    let (merged, terminal_merged, blocked) =
        apply_adoptions(calls, pools, members, &legacy, &terminal, &contained);
    if merged > 0 {
        progress!(
            "rustroo: merged {merged} groups, including {terminal_merged} guarded shortened-allele groups"
        );
    }
    if blocked > 0 {
        progress!(
            "rustroo: blocked {blocked} merges that would combine distinct co-occurring loci"
        );
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
    fn same_locus_adoptions(
        &self,
        calls: &[Call],
        pools: &Pools,
        members: &HashMap<u32, HashSet<u32>>,
    ) -> Vec<(u32, u32)> {
        let mut counts: HashMap<u32, HashMap<u32, usize>> = HashMap::default();
        let mut intact = HashMap::default();
        for call in calls {
            if !call.refound
                && *intact.entry(call.dna)
                    .or_insert_with(|| is_intact(&pools.dna[call.dna as usize]))
            {
                *counts.entry(call.cluster).or_default().entry(call.protein).or_default() += 1;
            }
        }
        let representatives: HashMap<u32, Vec<u32>> = counts
            .into_iter()
            .map(|(cluster, counts)| {
                let mut ranked: Vec<_> = counts.into_iter().collect();
                ranked.sort_by(|a, b| {
                    b.1.cmp(&a.1).then_with(|| {
                        pools.protein[a.0 as usize].cmp(&pools.protein[b.0 as usize])
                    })
                });
                ranked.truncate(3);
                (cluster, ranked.into_iter().map(|(protein, _)| protein).collect())
            })
            .collect();
        let mut ids: Vec<_> = self.occurrences.keys().copied().collect();
        ids.sort_unstable();
        ids.par_iter()
            .filter_map(|&child| {
                let occurrences = &self.occurrences[&child];
                let mut candidates = HashSet::default();
                let mut placed = 0;
                for occurrence in occurrences {
                    if let Some(hosts) = self.hosts(occurrence) {
                        candidates.extend(hosts.iter().copied());
                        placed += 1;
                    }
                }
                if placed < 2 {
                    return None;
                }
                let mut candidates: Vec<_> = candidates.into_iter().collect();
                candidates.sort_unstable();
                let mut target = None;
                let mut cache = HashMap::default();
                for parent in candidates {
                    if parent == child
                        || !self.occurrences.contains_key(&parent)
                    {
                        continue;
                    }
                    let (Some(witnesses), Some(reps)) =
                        (self.replicated.get(&parent), representatives.get(&parent))
                    else {
                        continue;
                    };
                    let context_matches = occurrences.iter().all(|occurrence| {
                        match occurrence.flanks {
                            (Some(a), Some(b)) => witnesses.contains(&(a, b)),
                            (None, Some(b)) if occurrence.upstream_edge => {
                                witnesses.iter().any(|&(_, flank)| flank == b)
                            }
                            (Some(a), None) if occurrence.downstream_edge => {
                                witnesses.iter().any(|&(flank, _)| flank == a)
                            }
                            // Missing placement is not evidence for a different locus.
                            (None, None) => true,
                            _ => false,
                        }
                    });
                    if !context_matches {
                        log::record(&format!(
                            "reconcile\trejected_context\tchild\t{child}\tparent\t{parent}"
                        ));
                        continue;
                    }
                    let sequence_matches = occurrences.iter().all(|occurrence| {
                        reps.iter().any(|&rep| {
                            *cache.entry((occurrence.protein, rep)).or_insert_with(|| {
                                Self::same_locus_sequence(
                                    &pools.protein[occurrence.protein as usize],
                                    &pools.protein[rep as usize],
                                )
                            })
                        })
                    });
                    if !sequence_matches {
                        log::record(&format!(
                            "reconcile\trejected_sequence\tchild\t{child}\tparent\t{parent}"
                        ));
                        continue;
                    }
                    if target.is_some() {
                        log::record(&format!(
                            "reconcile\tambiguous\tchild\t{child}\tadditional_target\t{parent}"
                        ));
                        return None;
                    }
                    target = Some(parent);
                }
                target.filter(|&parent| {
                    // Count competing hosts before overlap/direction rejection.
                    members[&child].is_disjoint(&members[&parent])
                        && (members[&parent].len(), std::cmp::Reverse(parent))
                            > (members[&child].len(), std::cmp::Reverse(child))
                }).map(|parent| {
                    log::record(&format!(
                        "reconcile\tproposed\tchild\t{child}\tparent\t{parent}\tplaced\t{placed}"
                    ));
                    (child, parent)
                })
            })
            .collect()
    }

    fn new(
        calls: &[Call],
        pools: &Pools,
        loci: &[Locus],
        contexts: &[(Option<u32>, Option<u32>)],
        contig_lengths: &[Vec<usize>],
    ) -> Self {
        let representatives = representative_ids(calls, pools);
        let mut occurrences: HashMap<u32, Vec<AlleleOccurrence>> = HashMap::default();
        let mut samples: HashMap<u32, HashSet<u32>> = HashMap::default();
        let mut ineligible = HashSet::default();
        let mut observed: HashMap<(u32, u32), HashSet<u32>> = HashMap::default();
        let mut votes: HashMap<(u32, u32), HashMap<u32, usize>> = HashMap::default();
        let mut intact_representatives = HashSet::default();
        let mut intact_cache = HashMap::default();
        for (position, (locus, &context)) in loci.iter().zip(contexts).enumerate() {
            let call = &calls[locus.first];
            let cluster = call.cluster;
            if locus.first != locus.last || !samples.entry(cluster).or_default().insert(call.sample)
            {
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
            // `context` is already (upstream, downstream); the contig-edge
            // tests need the coordinate-left and coordinate-right sides.
            let (left_flank, right_flank) = if call.strand >= 0 {
                (context.0, context.1)
            } else {
                (context.1, context.0)
            };
            let start_edge =
                first_on_contig
                    && left_flank.is_none()
                    && (call.start as usize).saturating_sub(1) <= ALLELE_EDGE_BP;
            let end_edge = last_on_contig
                && right_flank.is_none()
                && locus.end as usize <= length
                && length - locus.end as usize <= ALLELE_EDGE_BP;
            let (flanks, upstream_edge, downstream_edge) = if call.strand >= 0 {
                (context, start_edge, end_edge)
            } else {
                (context, end_edge, start_edge)
            };
            occurrences
                .entry(cluster)
                .or_default()
                .push(AlleleOccurrence {
                    protein: call.protein,
                    flanks,
                    upstream_edge,
                    downstream_edge,
                });
            let intact = *intact_cache
                .entry(call.dna)
                .or_insert_with(|| is_intact(&pools.dna[call.dna as usize]));
            let representative_length = pools.protein[representatives[&cluster] as usize].len();
            let full_length = (locus.first..=locus.last).any(|index| {
                pools.protein[calls[index].protein as usize].len() * 100
                    >= representative_length * 95
            });
            if full_length {
                if let (Some(upstream), Some(downstream)) = flanks {
                    observed
                        .entry((upstream, downstream))
                        .or_default()
                        .insert(cluster);
                    if locus.first == locus.last && intact {
                        *votes
                            .entry((upstream, downstream))
                            .or_default()
                            .entry(cluster)
                            .or_default() += 1;
                    }
                }
            }
            if call.protein == representatives[&cluster] && intact {
                intact_representatives.insert(cluster);
            }
        }
        occurrences.retain(|cluster, _| !ineligible.contains(cluster));
        let mut pairs: HashMap<(u32, u32), HashSet<u32>> = HashMap::default();
        let mut upstream: HashMap<u32, HashSet<u32>> = HashMap::default();
        let mut downstream: HashMap<u32, HashSet<u32>> = HashMap::default();
        let mut replicated: HashMap<u32, HashSet<(u32, u32)>> = HashMap::default();
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
        Self {
            occurrences,
            representatives,
            intact_representatives,
            pairs,
            upstream,
            downstream,
            replicated,
        }
    }

    fn same_locus_sequence(a: &[u8], b: &[u8]) -> bool {
        let (short, long) = if a.len() <= b.len() { (a, b) } else { (b, a) };
        if short.len() < ALLELE_MIN_AA {
            return false;
        }
        if short.len() * 5 >= long.len() * 4 {
            global_identity(short, long).is_some_and(|identity| identity >= CONTAIN_IDENTITY)
        } else {
            // Small shared domains are not enough; a fragment must reach a terminus.
            terminal_identity(short, long)
                .is_some_and(|identity| identity >= TERMINAL_IDENTITY)
        }
    }

    fn hosts(&self, occurrence: &AlleleOccurrence) -> Option<&HashSet<u32>> {
        match occurrence.flanks {
            (Some(upstream), Some(downstream)) => self.pairs.get(&(upstream, downstream)),
            (None, Some(downstream)) if occurrence.upstream_edge => {
                self.downstream.get(&downstream)
            }
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
        let mut cache: HashMap<(u32, u32, AlleleSide), bool> = HashMap::default();
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
                *cache
                    .entry((occurrence.protein, long_id, side))
                    .or_insert_with(|| {
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
                return AlleleDecision {
                    target: None,
                    ambiguous: true,
                };
            }
            matched = Some(candidate);
        }
        let Some(candidate) = matched else {
            return AlleleDecision::default();
        };
        if !self.occurrences.contains_key(&candidate)
            || !self.intact_representatives.contains(&candidate)
        {
            return AlleleDecision {
                target: None,
                ambiguous: true,
            };
        }
        let supported = self.replicated.get(&candidate).is_some_and(|witnesses| {
            occurrences
                .iter()
                .all(|occurrence| match occurrence.flanks {
                    (Some(upstream), Some(downstream)) => {
                        witnesses.contains(&(upstream, downstream))
                    }
                    (None, Some(downstream)) => {
                        witnesses.iter().any(|&(_, flank)| flank == downstream)
                    }
                    (Some(upstream), None) => witnesses.iter().any(|&(flank, _)| flank == upstream),
                    _ => false,
                })
        });
        AlleleDecision {
            target: supported.then_some(candidate),
            ambiguous: false,
        }
    }
}

fn shortened_identity(short: &[u8], long: &[u8], side: AlleleSide) -> Option<f64> {
    if short.len() < ALLELE_MIN_AA
        || short.len() >= long.len()
        || short.len() * 2 < long.len()
        || short
            .iter()
            .chain(long)
            .any(|&residue| residue == b'X' || residue == b'*')
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

#[derive(Clone, Copy, PartialEq, Eq)]
enum Adoption {
    /// Flank or abutment evidence. Overlap is allowed only for one physical locus.
    Legacy,
    /// A shortened allele. Any shared sample blocks the merge.
    Terminal,
    /// The short protein is part of the long one. A shared sample is allowed
    /// only when the short call sits inside the long call.
    Contained,
}

fn apply_adoptions(
    calls: &mut [Call],
    pools: &Pools,
    mut members: HashMap<u32, HashSet<u32>>,
    legacy: &[(u32, u32)],
    terminal: &[(u32, u32)],
    contained: &[(u32, u32)],
) -> (usize, usize, usize) {
    let mut parent: HashMap<u32, u32> = HashMap::default();
    let mut locations: HashMap<u32, Vec<usize>> = HashMap::default();
    for (index, call) in calls.iter().enumerate() {
        locations.entry(call.cluster).or_default().push(index);
    }
    let mut pair_cache = HashMap::default();
    let mut merged = 0;
    let mut terminal_merged = 0;
    let mut contained_merged = 0;
    let mut blocked = 0;
    for (child, parent_id, kind) in legacy
        .iter()
        .map(|&(child, parent)| (child, parent, Adoption::Legacy))
        .chain(
            terminal
                .iter()
                .map(|&(child, parent)| (child, parent, Adoption::Terminal)),
        )
        .chain(
            contained
                .iter()
                .map(|&(child, parent)| (child, parent, Adoption::Contained)),
        )
    {
        let root = find(&mut parent, parent_id);
        let child_root = find(&mut parent, child);
        if child_root == root {
            continue;
        }
        // Pairwise-disjoint proposals can still form a chain that joins two
        // copies in one sample. Check the accumulated groups, not just inputs.
        if !members[&child_root].is_disjoint(&members[&root]) {
            let allowed = match kind {
                Adoption::Terminal => false,
                Adoption::Legacy => cooccurring_loci_compatible(
                    &locations[&child_root],
                    &locations[&root],
                    &members[&child_root],
                    &members[&root],
                    calls,
                    pools,
                    &mut pair_cache,
                ),
                Adoption::Contained => intervals_inside(
                    &locations[&child_root],
                    &locations[&root],
                    Some(&members[&child_root]),
                    Some(&members[&root]),
                    calls,
                ),
            };
            if !allowed {
                log::record(&format!(
                    "reconcile\trejected_cooccurrence\tchild\t{child_root}\tparent\t{root}"
                ));
                blocked += 1;
                continue;
            }
        }
        let child_samples = members.remove(&child_root).expect("child root exists");
        members
            .get_mut(&root)
            .expect("parent root exists")
            .extend(child_samples);
        parent.insert(child_root, root);
        let child_locations = locations
            .remove(&child_root)
            .expect("child locations exist");
        locations
            .get_mut(&root)
            .expect("parent locations exist")
            .extend(child_locations);
        merged += 1;
        let method = match kind {
            Adoption::Legacy => "fragment",
            Adoption::Terminal => "same_locus_or_terminal",
            Adoption::Contained => "contained",
        };
        log::record(&format!(
            "reconcile\tapplied\tchild\t{child_root}\tparent\t{root}\tmethod\t{method}"
        ));
        terminal_merged += usize::from(kind == Adoption::Terminal);
        contained_merged += usize::from(kind == Adoption::Contained);
    }
    if contained_merged > 0 {
        progress!(
            "rustroo: adopted {contained_merged} length-contained fragments into the longer gene"
        );
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
    let mut first_context = HashMap::default();
    let mut parents = HashMap::default();
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
    let template_id = left
        .iter()
        .chain(right)
        .map(|&index| calls[index].protein)
        .max_by(|a, b| {
            pools.protein[*a as usize]
                .len()
                .cmp(&pools.protein[*b as usize].len())
                .then(b.cmp(a))
        })
        .expect("merge groups contain calls");
    let template = &pools.protein[template_id as usize];
    for &sample in left_samples.intersection(right_samples) {
        let mut indexes: Vec<usize> = left
            .iter()
            .chain(right)
            .copied()
            .filter(|&index| calls[index].sample == sample)
            .collect();
        indexes.sort_unstable();
        if !indexes
            .windows(2)
            .all(|pair| same_physical_locus(pair[0], pair[1], calls, pools, template, cache))
        {
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
        if fit.identity < CONTAIN_IDENTITY
            || fit.coverage < CONTAIN_COVERAGE
            || fit.internal_gap > 5
        {
            continue;
        }
        legacy_target = Some(long_id);
        legacy_matches += 1;
    }
    if terminal_matches == 1 {
        terminal_target.map(|cluster| MergeTarget {
            cluster,
            terminal: true,
        })
    } else if terminal_matches == 0 && legacy_matches == 1 {
        legacy_target.map(|cluster| MergeTarget {
            cluster,
            terminal: false,
        })
    } else {
        None
    }
}

/// A Bakta call that is only a piece of a longer gene.
///
/// Homology will not join them: the length ratio stays strict so a real
/// paralogue a few percent shorter is not swallowed. Adoption is sequence
/// containment at 90% identity over 90% of the short protein, and the placed
/// span must cover less than 90% of the long protein. A sample that has both
/// calls joins only when every short interval lies inside a long call on the
/// same contig and strand. Calls that never share a genome join as alleles.
/// Two eligible parents, or a short call sitting beside the long gene, stay
/// separate.
fn contained_adoptions(
    calls: &[Call],
    pools: &Pools,
    reps: &HashMap<u32, Vec<u8>>,
    members: &HashMap<u32, HashSet<u32>>,
    kmers: &HashMap<u32, HashSet<&[u8]>>,
) -> Vec<(u32, u32)> {
    let postings = containment_postings(kmers);
    let mut locations: HashMap<u32, Vec<usize>> = HashMap::default();
    for (index, call) in calls.iter().enumerate() {
        locations.entry(call.cluster).or_default().push(index);
    }
    let intact_hosts = intact_hosts(calls, pools);
    let mut ids: Vec<u32> = kmers.keys().copied().collect();
    ids.sort_unstable();
    let mut pairs: Vec<(u32, u32)> = ids
        .par_iter()
        .filter_map(|&short_id| {
            let short = reps.get(&short_id)?;
            let short_kmers = kmers.get(&short_id)?;
            let needed = short_kmers.len() as u32;
            if needed == 0 {
                return None;
            }
            let mut scores: HashMap<u32, u32> = HashMap::default();
            for window in short_kmers {
                let Some(hits) = postings.get(window) else {
                    continue;
                };
                for &long_id in hits {
                    if long_id != short_id {
                        *scores.entry(long_id).or_default() += 1;
                    }
                }
            }
            let mut ranked: Vec<(u32, u32)> = scores
                .into_iter()
                .filter(|&(long_id, score)| {
                    reps.get(&long_id)
                        .is_some_and(|long| long.len() > short.len())
                        && score * 2 >= needed
                })
                .collect();
            ranked.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
            ranked.truncate(8);
            let mut fits = Vec::new();
            for (long_id, _) in ranked {
                let long = &reps[&long_id];
                let Some(fit) = fit_short_in_long(short, long) else {
                    continue;
                };
                let span = fit.long_end.saturating_sub(fit.long_start);
                if fit.identity + f64::EPSILON < 0.90
                    || fit.coverage + f64::EPSILON < 0.90
                    || fit.internal_gap > 5
                    || span * 10 >= long.len() * 9
                {
                    continue;
                }
                fits.push(long_id);
            }
            // A second sequence match is ambiguity, including a damaged host
            // and one the coordinates later reject.
            if fits.len() != 1 {
                return None;
            }
            let long_id = fits[0];
            if !intact_hosts.contains(&long_id) {
                return None;
            }
            intervals_inside(
                locations.get(&short_id).map_or(&[], Vec::as_slice),
                locations.get(&long_id).map_or(&[], Vec::as_slice),
                members.get(&short_id),
                members.get(&long_id),
                calls,
            )
            .then_some((short_id, long_id))
        })
        .collect();
    pairs.sort_unstable();
    pairs
}

/// Families with at least one intact call carrying the family's own
/// representative protein, in one pass over the calls rather than a full
/// scan per adoption candidate. `pools.protein` is deduplicated by content,
/// so a pool id picks out the representative exactly as comparing its bytes
/// did, and `pools.dna` is shared by every call with the same sequence, so
/// each distinct allele is tested for an intact ORF once.
fn intact_hosts(calls: &[Call], pools: &Pools) -> HashSet<u32> {
    let representatives = representative_ids(calls, pools);
    let mut dna_intact: HashMap<u32, bool> = HashMap::default();
    let mut hosts = HashSet::default();
    for call in calls {
        if representatives.get(&call.cluster) != Some(&call.protein) {
            continue;
        }
        if *dna_intact
            .entry(call.dna)
            .or_insert_with(|| is_intact(&pools.dna[call.dna as usize]))
        {
            hosts.insert(call.cluster);
        }
    }
    hosts
}

fn containment_postings<'a>(
    kmers: &'a HashMap<u32, HashSet<&'a [u8]>>,
) -> HashMap<&'a [u8], Vec<u32>> {
    let mut postings: HashMap<&[u8], Vec<u32>> = HashMap::default();
    for (&id, windows) in kmers {
        for window in windows {
            postings.entry(window).or_default().push(id);
        }
    }
    // A 4-mer shared by dozens of families cannot point at one parent.
    postings.retain(|_, ids| ids.len() <= 64);
    postings
}

/// Every short call in a shared genome lies inside a long call.
/// No shared genome is containment, not a conflict.
fn intervals_inside(
    short_indexes: &[usize],
    long_indexes: &[usize],
    short_samples: Option<&HashSet<u32>>,
    long_samples: Option<&HashSet<u32>>,
    calls: &[Call],
) -> bool {
    let (Some(short_samples), Some(long_samples)) = (short_samples, long_samples) else {
        return true;
    };
    for &sample in short_samples.intersection(long_samples) {
        let mut pending = 0usize;
        for &index in short_indexes {
            let short = &calls[index];
            if short.sample != sample {
                continue;
            }
            pending += 1;
            let inside = long_indexes.iter().any(|&long_index| {
                let long = &calls[long_index];
                long.sample == sample
                    && long.contig == short.contig
                    && long.strand == short.strand
                    && short.start >= long.start
                    && short.end <= long.end
            });
            if !inside {
                return false;
            }
        }
        if pending == 0 {
            return false;
        }
    }
    true
}

fn terminal_identity(short: &[u8], long: &[u8]) -> Option<f64> {
    if short.is_empty() || short.len() >= long.len() {
        return None;
    }
    terminal_identity_parts(short, &[], long)
}

fn abutting_pieces(
    calls: &[Call],
    pools: &Pools,
    reps: &HashMap<u32, Vec<u8>>,
    members: &HashMap<u32, HashSet<u32>>,
) -> Vec<(u32, u32)> {
    let mut gaps: HashMap<(u32, u32), Vec<i32>> = HashMap::default();
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
    // One intact allele means the cluster is already a gene. Two of those
    // are an operon (metN then metP, rpsC then rplP), not the two halves of
    // a broken CDS. A fragment has no intact allele and may still complete
    // the one gene it sits against.
    // PERF PROTOTYPE: one ORF scan per distinct DNA id instead of per call.
    let mut dna_intact: HashMap<u32, bool> = HashMap::default();
    let intact: HashSet<u32> = calls
        .iter()
        .filter(|call| *dna_intact.entry(call.dna).or_insert_with(|| is_intact(&pools.dna[call.dna as usize])))
        .map(|call| call.cluster)
        .collect();
    let pairs: Vec<(u32, u32)> = keys
        .par_iter()
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
            if intact.contains(&left) && intact.contains(&right) {
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
        .collect();
    // A scrap between two genes would otherwise be a bridge: each edge is
    // legal (fragment plus one gene) and the union would still glue the
    // two genes into one reference.
    separate_finished_genes(pairs, &intact)
}

/// Keep fragment joins. Drop an edge that would put two intact genes in one
/// component, including through a chain of fragments.
fn separate_finished_genes(pairs: Vec<(u32, u32)>, intact: &HashSet<u32>) -> Vec<(u32, u32)> {
    let mut parent: HashMap<u32, u32> = HashMap::default();
    let mut finished: HashMap<u32, bool> = HashMap::default();
    let mut kept = Vec::with_capacity(pairs.len());
    for (child, parent_id) in pairs {
        let left = find(&mut parent, child);
        let right = find(&mut parent, parent_id);
        if left == right {
            continue;
        }
        let left_finished = *finished.entry(left).or_insert(intact.contains(&child));
        let right_finished = *finished.entry(right).or_insert(intact.contains(&parent_id));
        if left_finished && right_finished {
            continue;
        }
        parent.insert(left, right);
        finished.insert(right, left_finished || right_finished);
        finished.remove(&left);
        kept.push((child, parent_id));
    }
    kept
}

pub fn refind(
    calls: &mut Vec<Call>,
    pools: &mut Pools,
    contigs: &[Vec<PackedSeq>],
    min_cds_bp: usize,
) -> usize {
    sort_calls(calls);
    let loci = physical_loci(calls, pools);
    let contexts = locus_contexts(calls, &loci, |call| call.cluster);
    let mut support: HashMap<(u32, u32, u32), (usize, u32)> = HashMap::default();
    for (locus, context) in loci.iter().zip(contexts) {
        if let (Some(left), Some(right)) = context {
            // FIX(orientation): the gap scan below does not know the strand of
            // the missing gene, so the flank pair is unordered.
            let (left, right) = flank_key(left, right);
            record_context_sample(
                &mut support,
                (left, right, calls[locus.first].cluster),
                calls[locus.first].sample,
            );
        }
    }
    let mut target_support: HashMap<(u32, u32), Vec<(u32, usize)>> = HashMap::default();
    for ((left, right, middle), (samples, _)) in support {
        insert_context_target(
            target_support.entry((left, right)).or_default(),
            middle,
            samples,
        );
    }
    let target_of: HashMap<(u32, u32), Vec<u32>> = target_support
        .into_iter()
        .map(|(flanks, targets)| {
            (
                flanks,
                targets.into_iter().map(|(family, _)| family).collect(),
            )
        })
        .collect();
    if target_of.is_empty() {
        return 0;
    }
    let graph_targets = target_of.values().map(Vec::len).sum::<usize>();
    progress!(
        "rustroo: context graph has {} observed flank pairs and {graph_targets} candidate families",
        target_of.len()
    );
    let reps = representatives(calls, pools);
    let context_of: HashMap<u32, u32> = calls
        .iter()
        .map(|call| (call.cluster, call.context))
        .collect();
    let mut present: HashMap<u32, HashSet<u32>> = HashMap::default();
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
            refind_in_sample(
                &calls[start..end],
                sample,
                &target_of,
                &reps,
                contigs,
                local_present,
                min_cds_bp,
            )
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
            refound: true,
            fusion_piece: false,
            symbol: Label::default(),
            product: Label::default(),
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
/// start) for intergenic gaps that match a supported flank pair. The graph
/// stores a bounded set of middle-family alternatives; sequence evidence
/// chooses among them, while tied assignments to one ORF are rejected.
fn refind_in_sample(
    sample_calls: &[Call],
    sample: u32,
    target_of: &HashMap<(u32, u32), Vec<u32>>,
    reps: &HashMap<u32, Vec<u8>>,
    contigs: &[Vec<PackedSeq>],
    present: HashSet<u32>,
    min_cds_bp: usize,
) -> Vec<FoundGene> {
    // FIX(orientation): every gap is scored first, then each missing gene is
    // placed in its single best gap. The old first-gap-along-the-contig rule
    // picked a different gap when a contig was reverse-complemented. The
    // choice is by identity, then contig, then distance from the nearer contig
    // end, then the ORF bases: all unchanged by flipping a contig.
    struct Candidate {
        identity: f64,
        contig: u32,
        edge_distance: usize,
        gene: FoundGene,
    }
    let mut best: HashMap<u32, Candidate> = HashMap::default();
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
            let Some(targets) = target_of.get(&flank_key(pair[0].cluster, pair[1].cluster)) else {
                continue;
            };
            // A missing contig reads as "" in the original, which always
            // fails the range check below.
            let Some(sequence) = contig_sequence(contigs, sample, contig) else {
                continue;
            };
            let from = pair[0].end as usize;
            let to = (pair[1].start as usize)
                .saturating_sub(1)
                .min(sequence.len());
            if from >= to || to > sequence.len() {
                continue;
            }
            let gap_seq = sequence.slice(from, to);
            let mut gap_best: Option<Candidate> = None;
            let mut ambiguous = false;
            for &target in targets {
                if present.contains(&target) {
                    continue;
                }
                let Some(target_protein) = reps.get(&target) else {
                    continue;
                };
                let Some((identity, orf)) = best_orf(&gap_seq, target_protein, min_cds_bp) else {
                    continue;
                };
                let protein = translate(&orf.dna);
                let start = pair[0].end + orf.offset + 1;
                let end_bp = pair[0].end + orf.offset + orf.dna.len() as u32;
                let edge_distance = (start as usize - 1).min(sequence.len() - end_bp as usize);
                let candidate = Candidate {
                    identity,
                    contig,
                    edge_distance,
                    gene: FoundGene {
                        sample,
                        contig,
                        start,
                        end: end_bp,
                        strand: orf.strand,
                        dna: orf.dna,
                        protein,
                        cluster: target,
                    },
                };
                match &gap_best {
                    None => gap_best = Some(candidate),
                    Some(current) => match candidate.identity.total_cmp(&current.identity) {
                        std::cmp::Ordering::Greater => {
                            gap_best = Some(candidate);
                            ambiguous = false;
                        }
                        std::cmp::Ordering::Equal
                            if candidate.gene.cluster != current.gene.cluster =>
                        {
                            ambiguous = true;
                        }
                        _ => {}
                    },
                };
            }
            if ambiguous {
                continue;
            }
            let Some(candidate) = gap_best else {
                continue;
            };
            let target = candidate.gene.cluster;
            let better = match best.get(&target) {
                None => true,
                Some(current) => candidate
                    .identity
                    .total_cmp(&current.identity)
                    .then(current.contig.cmp(&candidate.contig))
                    .then(current.edge_distance.cmp(&candidate.edge_distance))
                    .then(current.gene.dna.cmp(&candidate.gene.dna))
                    .is_gt(),
            };
            if better {
                best.insert(target, candidate);
            }
        }
        index = end;
    }
    let mut found: Vec<FoundGene> = best.into_values().map(|candidate| candidate.gene).collect();
    // One gene per target: ordering by target keeps pool interning order
    // independent of contig orientation.
    found.sort_by_key(|gene| gene.cluster);
    found
}

const MIN_CONTEXT_SAMPLES: usize = 1;
const MAX_CONTEXT_TARGETS: usize = 8;

fn record_context_sample(
    support: &mut HashMap<(u32, u32, u32), (usize, u32)>,
    key: (u32, u32, u32),
    sample: u32,
) {
    // Calls are sorted by sample, so one last-sample value replaces a per-edge
    // sample set while still counting each genome only once.
    let entry = support.entry(key).or_insert((1, sample));
    if entry.1 != sample {
        entry.0 += 1;
        entry.1 = sample;
    }
}

fn insert_context_target(targets: &mut Vec<(u32, usize)>, family: u32, support: usize) {
    if support < MIN_CONTEXT_SAMPLES {
        return;
    }
    let position = targets
        .iter()
        .position(|(existing, count)| {
            support > *count || (support == *count && family < *existing)
        })
        .unwrap_or(targets.len());
    if position >= MAX_CONTEXT_TARGETS {
        return;
    }
    targets.insert(position, (family, support));
    targets.truncate(MAX_CONTEXT_TARGETS);
}

/// Unordered flank pair: a gene between A and B is between B and A on the
/// reverse-complemented contig.
fn flank_key(a: u32, b: u32) -> (u32, u32) {
    (a.min(b), a.max(b))
}

struct FoundOrf {
    dna: String,
    offset: u32,
    strand: i8,
}

fn best_orf(gap: &str, target: &[u8], min_cds_bp: usize) -> Option<(f64, FoundOrf)> {
    let reverse = crate::genetic_code::reverse_complement(gap);
    // Short matches cannot stand in for a gene. Rank accepted ORFs by
    // identity weighted by shorter/longer coverage, not identity alone.
    let mut spans: Vec<(i8, usize, usize, u32)> = Vec::new();
    // First stop codon at or after each codon of the current frame, swept
    // once from the back instead of re-scanned from every start codon.
    // `codons` stands for "no stop after here", which is what the forward
    // scan reported by running off the end.
    let mut next_stop: Vec<usize> = Vec::new();
    for (sequence, strand) in [(gap, 1i8), (reverse.as_str(), -1)] {
        for frame in 0..3 {
            if sequence.len() <= frame {
                continue;
            }
            let frame_seq = &sequence[frame..];
            let codons = frame_seq.len() / 3;
            next_stop.clear();
            next_stop.resize(codons + 1, codons);
            for index in (0..codons).rev() {
                let at = index * 3;
                next_stop[index] = if is_stop(&frame_seq.as_bytes()[at..at + 3]) {
                    index
                } else {
                    next_stop[index + 1]
                };
            }
            let mut codon = 0;
            while codon + 1 < codons {
                let start = codon * 3;
                let codon_bytes = &frame_seq.as_bytes()[start..start + 3];
                if !is_start(codon_bytes) {
                    codon += 1;
                    continue;
                }
                let stop_at = next_stop[codon + 1];
                if stop_at == codons {
                    codon += 1;
                    continue;
                }
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
    let score =
        |&(strand, start, end, offset): &(i8, usize, usize, u32)| -> Option<(f64, FoundOrf)> {
            let sequence = if strand > 0 { gap } else { reverse.as_str() };
            if end - start < min_cds_bp {
                return None;
            }
            let protein = translate_bytes(&sequence.as_bytes()[start..end]);
            let ratio = crate::align::length_ratio(&protein, target);
            if ratio < 0.80 {
                return None;
            }
            let identity = global_identity(&protein, target)?;
            if identity < 0.8 {
                return None;
            }
            Some((
                identity * ratio,
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
    // FIX(orientation): an identity tie goes to the smaller ORF sequence and
    // then to the ORF nearer a gap end, not to whichever strand was scanned
    // first, so a reverse-complemented gap gives the same ORF.
    let gap_len = gap.len() as u32;
    let edge = |orf: &FoundOrf| orf.offset.min(gap_len.saturating_sub(orf.offset + orf.dna.len() as u32));
    let mut best: Option<(f64, FoundOrf)> = None;
    for found in scored.into_iter().flatten() {
        let better = match &best {
            None => true,
            Some((identity, current)) => found
                .0
                .total_cmp(identity)
                .then(current.dna.cmp(&found.1.dna))
                .then(edge(current).cmp(&edge(&found.1)))
                .is_gt(),
        };
        if better {
            best = Some(found);
        }
    }
    best
}

/// The two genome-exact pieces of a fused CDS cut after `cut_aa` codons.
pub(crate) struct FusionPieces {
    pub prefix: String,
    pub prefix_span: (u32, u32),
    pub suffix: String,
    pub suffix_span: (u32, u32),
}

/// Cut `dna` (the call's 5'->3' coding sequence) after `cut_aa` codons.
///
/// Both pieces are plain substrings of the CDS: the prefix ends where the
/// fused reading frame continued (no stop is appended), the suffix runs to the
/// CDS's own last base and so keeps its real stop codon. Coordinates are
/// counted from the 5' end (`start` on +, `end` on -), so the two intervals
/// tile `start..=end` with no gap and no overlap. `None` unless the CDS is
/// exactly the genome interval it is annotated on, base for base.
pub(crate) fn fusion_pieces(
    call: &Call,
    dna: &str,
    cut_aa: usize,
    contigs: &[Vec<PackedSeq>],
) -> Option<FusionPieces> {
    let cut_bp = cut_aa.checked_mul(3)?;
    if cut_bp == 0 || cut_bp >= dna.len() || call.end < call.start {
        return None;
    }
    let span = (call.end - call.start) as usize + 1;
    if span != dna.len() {
        return None;
    }
    let sequence = contig_sequence(contigs, call.sample, call.contig)?;
    if call.start == 0 || call.end as usize > sequence.len() {
        return None;
    }
    let genome = sequence.slice(call.start as usize - 1, call.end as usize);
    let genome = if call.strand < 0 {
        crate::genetic_code::reverse_complement(&genome)
    } else {
        genome
    };
    if genome != dna {
        return None;
    }
    let cut = cut_bp as u32;
    let (prefix_span, suffix_span) = if call.strand >= 0 {
        ((call.start, call.start + cut - 1), (call.start + cut, call.end))
    } else {
        ((call.end - cut + 1, call.end), (call.start, call.end - cut))
    };
    Some(FusionPieces {
        prefix: dna[..cut_bp].to_string(),
        prefix_span,
        suffix: dna[cut_bp..].to_string(),
        suffix_span,
    })
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
    representative_ids(calls, pools)
        .into_iter()
        .map(|(cluster, protein)| (cluster, pools.protein[protein as usize].clone()))
        .collect()
}

fn representative_ids(calls: &[Call], pools: &Pools) -> HashMap<u32, u32> {
    let mut best = HashMap::default();
    for call in calls {
        let slot = best.entry(call.cluster).or_insert(call.protein);
        // FIX(orientation): equal lengths go to the smaller sequence, not to
        // whichever copy came first in coordinate order.
        let (candidate, current) = (&pools.protein[call.protein as usize], &pools.protein[*slot as usize]);
        if candidate.len() > current.len() || (candidate.len() == current.len() && candidate < current) {
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
    let mut sets: HashMap<u32, HashSet<u32>> = HashMap::default();
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
    calls: &[Call],
    loci: &[Locus],
    contexts: &[(Option<u32>, Option<u32>)],
) -> HashMap<u32, (Option<u32>, Option<u32>)> {
    let mut votes: HashMap<u32, HashMap<(Option<u32>, Option<u32>), u32>> = HashMap::default();
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
            left: HashMap::default(),
            right: HashMap::default(),
            adjacent: HashMap::default(),
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

    fn candidates(&self, id: u32, flanks: &HashMap<u32, (Option<u32>, Option<u32>)>) -> Vec<u32> {
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
    let mut pairs = HashSet::default();
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
    let shorter =
        (i64::from(a.end) - i64::from(a.start) + 1).min(i64::from(b.end) - i64::from(b.start) + 1);
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
        global_identity(
            &pools.protein[a.protein as usize],
            &pools.protein[b.protein as usize],
        )
        .is_some_and(|identity| identity < 0.4)
    })
}

fn physical_loci(calls: &[Call], pools: &Pools) -> Vec<Locus> {
    let representatives = representative_ids(calls, pools);
    let mut loci: Vec<Locus> = Vec::new();
    let mut cache = HashMap::default();
    for (index, call) in calls.iter().enumerate() {
        if let Some(previous) = loci.last_mut() {
            if calls[previous.first].cluster == call.cluster
                && same_physical_locus(
                    previous.last,
                    index,
                    calls,
                    pools,
                    &pools.protein[representatives[&call.cluster] as usize],
                    &mut cache,
                )
            {
                previous.last = index;
                previous.end = previous.end.max(call.end);
                continue;
            }
        }
        loci.push(Locus {
            first: index,
            last: index,
            end: call.end,
        });
    }
    loci
}

/// Whether a finished family was annotated, and whether its copies share one
/// neighborhood. Calls must already be ordered by sample, contig, and start.
#[derive(Clone, Copy)]
pub struct ClusterFacts {
    pub annotated: bool,
    pub neighbors_agree: bool,
}

pub fn cluster_facts(calls: &[Call], pools: &Pools) -> HashMap<u32, ClusterFacts> {
    debug_assert!(calls.windows(2).all(|pair| {
        pair[0].sample < pair[1].sample
            || (pair[0].sample == pair[1].sample && pair[0].contig < pair[1].contig)
            || (pair[0].sample == pair[1].sample
                && pair[0].contig == pair[1].contig
                && pair[0].start <= pair[1].start)
    }));
    let mut facts: HashMap<u32, ClusterFacts> = HashMap::default();
    for call in calls {
        let fact = facts.entry(call.cluster).or_insert(ClusterFacts {
            annotated: false,
            neighbors_agree: false,
        });
        if !call.refound {
            fact.annotated = true;
        }
    }
    let loci = physical_loci(calls, pools);
    let contexts = locus_contexts(calls, &loci, |call| call.cluster);
    let mut seen: HashMap<u32, Vec<(Option<u32>, Option<u32>)>> = HashMap::default();
    for (locus, context) in loci.iter().zip(contexts) {
        seen.entry(calls[locus.first].cluster)
            .or_default()
            .push(context);
    }
    for (cluster, contexts) in seen {
        if let Some(fact) = facts.get_mut(&cluster) {
            fact.neighbors_agree = same_neighborhood(&contexts);
        }
    }
    facts
}

/// One neighborhood covers more than half of the copies that have a neighbor.
/// A single copy between two genes agrees with itself. Copies with no neighbor
/// on either side do not count as support.
fn same_neighborhood(contexts: &[(Option<u32>, Option<u32>)]) -> bool {
    let mut counts: HashMap<(Option<u32>, Option<u32>), usize> = HashMap::default();
    let mut known = 0usize;
    for &context in contexts {
        if context.0.is_none() && context.1.is_none() {
            continue;
        }
        known += 1;
        *counts.entry(context).or_default() += 1;
    }
    let best = counts.into_values().max().unwrap_or(0);
    known > 0 && best * 2 > known
}

/// `(upstream, downstream)` neighbor families of every locus, in the
/// locus's own reading direction (5' side first).
fn locus_contexts(
    calls: &[Call],
    loci: &[Locus],
    family: impl Fn(&Call) -> u32,
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
                let gap = i64::from(calls[loci[cursor].first].start)
                    - i64::from(loci[cursor - 1].end)
                    - 1;
                if (-NEIGHBOR_OVERLAP..=i64::from(NEIGHBOR_GAP)).contains(&gap) {
                    Some(family(&calls[loci[cursor - 1].first]))
                } else {
                    None
                }
            } else {
                None
            };
            let next = if cursor + 1 < end {
                let gap = i64::from(calls[loci[cursor + 1].first].start)
                    - i64::from(loci[cursor].end)
                    - 1;
                if (-NEIGHBOR_OVERLAP..=i64::from(NEIGHBOR_GAP)).contains(&gap) {
                    Some(family(&calls[loci[cursor + 1].first]))
                } else {
                    None
                }
            } else {
                None
            };
            // FIX(orientation): flanks are read in the gene's own direction,
            // (upstream, downstream). A contig assembled in the opposite
            // orientation flips both coordinate order and strand, so the
            // same locus gets the same context either way.
            contexts[cursor] = if calls[loci[cursor].first].strand < 0 {
                (next, prev)
            } else {
                (prev, next)
            };
        }
        index = end;
    }
    contexts
}

fn contig_sequence(contigs: &[Vec<PackedSeq>], sample: u32, contig: u32) -> Option<&PackedSeq> {
    contigs
        .get(sample as usize)
        .and_then(|sample| sample.get(contig as usize))
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

/// Renumber clusters densely from 0.
///
/// FIX(orientation): ids used to follow first appearance in coordinate order,
/// so reverse-complementing one contig of the first sample renumbered
/// families, and every later tie broken "by the smaller id" (and the order
/// names are handed out in) changed with it. Families are now ordered by
/// where they first occur by sample and contig, then by the smallest member
/// protein by content, then by size; position along a contig is only the
/// last resort for families that are otherwise indistinguishable.
pub fn compress_clusters(calls: &mut [Call], pools: &Pools) {
    sort_calls(calls);
    let rank = protein_ranks(pools);
    // cluster -> (sample, contig, protein rank, -calls, first index)
    let mut keys: HashMap<u32, (u32, u32, u32, i64, usize)> = HashMap::default();
    for (index, call) in calls.iter().enumerate() {
        let protein = rank[call.protein as usize];
        let entry = keys
            .entry(call.cluster)
            .or_insert((call.sample, call.contig, protein, 0, index));
        if (call.sample, call.contig) == (entry.0, entry.1) && protein < entry.2 {
            entry.2 = protein;
        }
        entry.3 -= 1;
    }
    let mut ordered: Vec<(u32, (u32, u32, u32, i64, usize))> = keys.into_iter().collect();
    ordered.sort_unstable_by(|a, b| a.1.cmp(&b.1));
    let map: HashMap<u32, u32> = ordered
        .iter()
        .enumerate()
        .map(|(id, &(cluster, _))| (cluster, id as u32))
        .collect();
    for call in calls.iter_mut() {
        call.cluster = map[&call.cluster];
    }
}

/// Rank of every pooled protein by content (length, then bytes), so a
/// tie-break can use the sequence rather than the order it was first seen.
/// The ranking is a pure function of the pool, and `compress_clusters` runs
/// after every synteny round, so the pool keeps the result between the
/// rounds that do not intern anything new.
fn protein_ranks(pools: &Pools) -> std::sync::Arc<[u32]> {
    pools.protein_ranks(|proteins: &[Vec<u8>]| {
        let mut order: Vec<u32> = (0..proteins.len() as u32).collect();
        order.par_sort_unstable_by(|&a, &b| {
            let (pa, pb) = (&proteins[a as usize], &proteins[b as usize]);
            pa.len().cmp(&pb.len()).then_with(|| pa.cmp(pb)).then(a.cmp(&b))
        });
        let mut rank = vec![0u32; order.len()];
        for (position, &protein) in order.iter().enumerate() {
            rank[protein as usize] = position as u32;
        }
        rank
    })
}

/// Protein identity of two alleles over the shorter one: a banded global
/// alignment when their lengths are comparable, the short-in-long fit when
/// one is a fragment of the other. `0` when no alignment is found.
pub fn allele_identity(a: &[u8], b: &[u8]) -> f64 {
    if a == b {
        return 1.0;
    }
    let (short, long) = if a.len() <= b.len() { (a, b) } else { (b, a) };
    if short.is_empty() {
        return 0.0;
    }
    if (long.len() - short.len()) * 5 <= long.len() {
        if let Some(identity) = global_identity(a, b) {
            return identity;
        }
    }
    fit_short_in_long(short, long)
        .map(|fit| fit.identity)
        .unwrap_or(0.0)
}

/// What `split_divergent_members` did.
#[derive(Default, Debug, Clone, PartialEq)]
pub struct DivergentSplit {
    /// Families that lost at least one divergent group.
    pub families: usize,
    /// New families created from divergent groups.
    pub new_families: usize,
    /// Physical loci moved into the new families.
    pub loci: usize,
    /// Divergent groups left in place because they sit in the family's own
    /// neighbourhood and never share a genome with it (divergent alleles of
    /// one locus, not paralogs).
    pub kept_alleles: usize,
    /// Divergent groups left in place because they have no neighbour at all.
    pub kept_unplaced: usize,
}

/// FIX(paralogs): split distant paralogs that homology and adoption put in
/// one family even when no genome carries two copies.
///
/// `split_paralogs` only runs for families with a multi-copy genome, so two
/// paralogs that never co-occur stayed one gene however different they are.
/// Here every locus whose longest protein is below `threshold` identity to the
/// family representative (the longest member protein, as everywhere else) is
/// a divergent candidate. Candidates are grouped greedily, longest first, at
/// the same threshold against each group's own representative (the CD-HIT
/// rule used by homology). A group becomes its own family when synteny agrees
/// it is another locus:
///   * it shares a genome with the rest of the family (two copies), or
///   * none of its gene neighbours is a neighbour of the rest of the family.
/// A divergent group in the family's own neighbourhood that never co-occurs
/// with it is an allele of the same locus and stays. A group with no
/// neighbour at all has no synteny evidence and stays.
///
/// `threshold` defaults to the family identity (0.70): homology already
/// promises that much identity to the representative, so a member below it
/// was only ever joined through adoption or refinding.
pub fn split_divergent_members(calls: &mut [Call], pools: &Pools, threshold: f64) -> DivergentSplit {
    if threshold <= 0.0 || calls.is_empty() {
        return DivergentSplit::default();
    }
    sort_calls(calls);
    let loci = physical_loci(calls, pools);
    let contexts = locus_contexts(calls, &loci, |call| call.cluster);
    let reps = representative_ids(calls, pools);
    let mut by_cluster: HashMap<u32, Vec<usize>> = HashMap::default();
    for (index, locus) in loci.iter().enumerate() {
        by_cluster.entry(calls[locus.first].cluster).or_default().push(index);
    }
    // The longest protein of each locus (a broken locus is judged by its
    // longest piece), smaller sequence on a length tie.
    let locus_proteins: Vec<u32> = loci
        .iter()
        .map(|locus| {
            (locus.first..=locus.last)
                .map(|index| calls[index].protein)
                .max_by(|&a, &b| {
                    let (pa, pb) = (&pools.protein[a as usize], &pools.protein[b as usize]);
                    pa.len().cmp(&pb.len()).then_with(|| pb.cmp(pa))
                })
                .expect("a locus has at least one call")
        })
        .collect();
    let mut clusters: Vec<u32> = by_cluster.keys().copied().collect();
    clusters.sort_unstable();
    struct Plan {
        groups: Vec<Vec<usize>>,
        kept_alleles: usize,
        kept_unplaced: usize,
    }
    let plans: Vec<Plan> = clusters
        .par_iter()
        .map(|cluster| {
            let members = &by_cluster[cluster];
            let rep = &pools.protein[reps[cluster] as usize];
            let mut identity_of: HashMap<u32, f64> = HashMap::default();
            let mut divergent: Vec<usize> = Vec::new();
            for &index in members {
                let protein = locus_proteins[index];
                let identity = *identity_of
                    .entry(protein)
                    .or_insert_with(|| allele_identity(&pools.protein[protein as usize], rep));
                if identity < threshold {
                    divergent.push(index);
                }
            }
            let mut plan = Plan {
                groups: Vec::new(),
                kept_alleles: 0,
                kept_unplaced: 0,
            };
            if divergent.is_empty() {
                return plan;
            }
            // Greedy longest-first grouping of the divergent proteins.
            let mut proteins: Vec<u32> = divergent.iter().map(|&index| locus_proteins[index]).collect();
            proteins.sort_unstable_by(|&a, &b| {
                let (pa, pb) = (&pools.protein[a as usize], &pools.protein[b as usize]);
                pb.len().cmp(&pa.len()).then_with(|| pa.cmp(pb))
            });
            proteins.dedup();
            let mut group_reps: Vec<u32> = Vec::new();
            let mut group_of: HashMap<u32, usize> = HashMap::default();
            for &protein in &proteins {
                let sequence = &pools.protein[protein as usize];
                let joined = group_reps.iter().position(|&head| {
                    allele_identity(sequence, &pools.protein[head as usize]) >= threshold
                });
                let group = joined.unwrap_or_else(|| {
                    group_reps.push(protein);
                    group_reps.len() - 1
                });
                group_of.insert(protein, group);
            }
            let mut grouped: Vec<Vec<usize>> = vec![Vec::new(); group_reps.len()];
            for &index in &divergent {
                grouped[group_of[&locus_proteins[index]]].push(index);
            }
            let moved: HashSet<usize> = divergent.iter().copied().collect();
            let rest: Vec<usize> = members.iter().copied().filter(|index| !moved.contains(index)).collect();
            let flanks_of = |indexes: &[usize]| -> HashSet<u32> {
                indexes
                    .iter()
                    .flat_map(|&index| [contexts[index].0, contexts[index].1])
                    .flatten()
                    .collect()
            };
            let samples_of = |indexes: &[usize]| -> HashSet<u32> {
                indexes.iter().map(|&index| calls[loci[index].first].sample).collect()
            };
            let rest_flanks = flanks_of(&rest);
            let rest_samples = samples_of(&rest);
            for group in grouped {
                let flanks = flanks_of(&group);
                let shares_genome = !samples_of(&group).is_disjoint(&rest_samples);
                let other_neighbourhood = !flanks.is_empty() && flanks.is_disjoint(&rest_flanks);
                if rest.is_empty() {
                    // Every locus is divergent from the representative's own
                    // protein only if the representative is itself one of
                    // them; nothing to split from.
                    continue;
                }
                if shares_genome || other_neighbourhood {
                    plan.groups.push(group);
                } else if flanks.is_empty() {
                    plan.kept_unplaced += 1;
                } else {
                    plan.kept_alleles += 1;
                }
            }
            plan
        })
        .collect();
    let mut next_id = calls.iter().map(|call| call.cluster).max().unwrap_or(0) + 1;
    let mut result = DivergentSplit::default();
    for (cluster, plan) in clusters.iter().zip(plans) {
        result.kept_alleles += plan.kept_alleles;
        result.kept_unplaced += plan.kept_unplaced;
        if plan.groups.is_empty() {
            continue;
        }
        result.families += 1;
        for group in plan.groups {
            let id = next_id;
            next_id += 1;
            result.new_families += 1;
            result.loci += group.len();
            let mut samples = HashSet::default();
            for &index in &group {
                samples.insert(calls[loci[index].first].sample);
                for call in loci[index].first..=loci[index].last {
                    calls[call].cluster = id;
                }
            }
            if log::enabled() {
                let protein = locus_proteins[group[0]];
                let identity = allele_identity(
                    &pools.protein[protein as usize],
                    &pools.protein[reps[cluster] as usize],
                );
                log::record(&format!(
                    "divergent_split\tfamily\t{cluster}\tsymbol\t{}\tloci\t{}\tsamples\t{}\tidentity_to_representative\t{identity:.3}",
                    log_field(&calls[loci[group[0]].first].symbol),
                    group.len(),
                    samples.len(),
                ));
            }
        }
    }
    result
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
    let mut rows: HashMap<u32, Row> = HashMap::default();
    for call in calls.iter().filter(|call| clusters.contains(&call.cluster)) {
        let row = rows.entry(call.cluster).or_insert_with(|| Row {
            symbol: String::new(),
            product: String::new(),
            samples: HashSet::default(),
            copies: 0,
            longest: 0,
        });
        if row.symbol.is_empty() {
            row.symbol = call.symbol.to_string();
        }
        if row.product.is_empty() {
            row.product = call.product.to_string();
        }
        row.samples.insert(call.sample);
        row.copies += 1;
        row.longest = row.longest.max(pools.dna[call.dna as usize].len());
    }
    log::record("removed\treason\tcluster\tsymbol\tproduct\tsamples\tcopies\tlongest_bp");
    let mut ordered: Vec<_> = rows.into_iter().collect();
    ordered.sort_by(|left, right| {
        left.1
            .symbol
            .cmp(&right.1.symbol)
            .then(left.0.cmp(&right.0))
    });
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
    let mut intact = HashSet::default();
    let mut seen = HashSet::default();
    let mut usable: HashMap<u32, bool> = HashMap::default();
    for call in calls {
        seen.insert(call.cluster);
        if *usable.entry(call.dna).or_insert_with(|| !sequence_cannot_found_a_family(&pools.dna[call.dna as usize])) {
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
        return HashSet::default();
    }
    let internal = families_between_neighbors(calls, loci);
    let mut consumed = HashSet::default();
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
        consume_tail(
            &mut consumed,
            calls,
            loci,
            index..end,
            true,
            &support,
            &internal,
            threshold,
        );
        consume_tail(
            &mut consumed,
            calls,
            loci,
            index..end,
            false,
            &support,
            &internal,
            threshold,
        );
        index = end;
    }
    // A family stays if any copy was not on a tail walk. Dropping the family
    // would also delete that internal or gap-isolated copy.
    let mut total: HashMap<u32, usize> = HashMap::default();
    let mut hit: HashMap<u32, usize> = HashMap::default();
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
    let mut internal = HashSet::default();
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
    let n_samples = calls
        .iter()
        .map(|call| call.sample as usize + 1)
        .max()
        .unwrap_or(0);
    let threshold = weak_flank_threshold(n_samples);
    let mut support: HashMap<(u32, Option<u32>, Option<u32>), HashSet<u32>> = HashMap::default();
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
    // PERF PROTOTYPE: the alternatives for one family are the support keys
    // with that family; grouping them once replaces a scan of every key.
    let mut by_family: HashMap<u32, Vec<(Option<u32>, Option<u32>, usize)>> = HashMap::default();
    for (&(family, prev, next), samples) in &support {
        by_family.entry(family).or_default().push((prev, next, samples.len()));
    }
    contexts
        .iter()
        .enumerate()
        .map(|(index, context)| {
            let family = family_of(&calls[loci[index].first]);
            if count(family, context.0, context.1) >= threshold {
                return *context;
            }
            let rows = by_family.get(&family).map(Vec::as_slice).unwrap_or(&[]);
            let next = supported_alternative_indexed(context.0, context.1, true, threshold, rows);
            let prev = supported_alternative_indexed(context.0, context.1, false, threshold, rows);
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

fn supported_alternative_indexed(
    prev: Option<u32>,
    next: Option<u32>,
    fix_next: bool,
    threshold: usize,
    rows: &[(Option<u32>, Option<u32>, usize)],
) -> Option<u32> {
    let observed = if fix_next { next } else { prev };
    observed?;
    let mut alternatives: Vec<(usize, u32)> = rows
        .iter()
        .filter(|(pair_prev, pair_next, _)| if fix_next { *pair_prev == prev } else { *pair_next == next })
        .filter_map(|&(pair_prev, pair_next, count)| {
            let neighbor = if fix_next { pair_next } else { pair_prev }?;
            (count >= threshold).then_some((count, neighbor))
        })
        .collect();
    alternatives.sort_unstable_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
    alternatives.dedup_by(|a, b| a.1 == b.1);
    match alternatives.as_slice() {
        [(count, neighbor)] if observed != Some(*neighbor) && *count >= threshold => Some(*neighbor),
        _ => observed,
    }
}

#[allow(dead_code)]
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
            *id == family
                && if fix_next {
                    *pair_prev == prev
                } else {
                    *pair_next == next
                }
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
        [(count, neighbor)] if observed != Some(*neighbor) && *count >= threshold => {
            Some(*neighbor)
        }
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

    #[test]
    fn refind_rejects_tiny_and_low_coverage_orfs() {
        let mut full = protein(100, 42);
        full[0] = b'M';
        let tiny = back_translate(&full[..5]);
        assert!(best_orf(&tiny, &full, 90).is_none());
        assert!(best_orf(&tiny, &full, 0).is_none());
        let fragment = back_translate(&full[..40]);
        assert!(best_orf(&fragment, &full, 90).is_none());
        let dna = back_translate(&full);
        assert_eq!(best_orf(&dna, &full, 90).unwrap().1.dna, dna);
        assert!(best_orf(&dna, &full, dna.len() + 1).is_none());
    }

    #[test]
    fn refind_prefers_whole_gene_over_perfect_short_match_and_is_orientation_invariant() {
        let mut full = protein(100, 42);
        full[0] = b'M';
        let tiny = back_translate(&full[..5]);
        let dna = back_translate(&full);
        let gap = format!("{tiny}CCC{dna}");
        let found = best_orf(&gap, &full, 90).unwrap().1;
        assert_eq!(found.dna, dna);
        let reverse = crate::genetic_code::reverse_complement(&gap);
        let flipped = best_orf(&reverse, &full, 90).unwrap().1;
        assert_eq!(flipped.dna, dna);
        assert_eq!(found.strand, -flipped.strand);
    }

    #[test]
    fn refind_uses_sequence_to_choose_between_supported_context_families() {
        let mut left = protein(40, 4);
        let mut right = protein(40, 5);
        let mut first_middle = protein(100, 6);
        let mut second_middle = protein(100, 7);
        for sequence in [&mut left, &mut right, &mut first_middle, &mut second_middle] {
            sequence[0] = b'M';
        }
        let left_dna = back_translate(&left);
        let middle_dna = back_translate(&second_middle);
        let right_dna = back_translate(&right);
        let spacer = "N".repeat(60);
        let right_start = (left_dna.len() + spacer.len() + middle_dna.len() + spacer.len() + 1) as u32;
        let contig = format!("{left_dna}{spacer}{middle_dna}{spacer}{right_dna}");

        let mut pools = Pools::default();
        let calls = vec![
            call(4, 10, 1, &left, &mut pools),
            call(4, 11, right_start, &right, &mut pools),
        ];
        let targets = HashMap::from_iter([(
            flank_key(10, 11),
            vec![20, 21],
        )]);
        let reps = HashMap::from_iter([(20, first_middle), (21, second_middle)]);
        let mut contigs: Vec<Vec<PackedSeq>> = (0..5).map(|_| Vec::new()).collect();
        contigs[4].push(PackedSeq::new(contig));

        let found = refind_in_sample(
            &calls,
            4,
            &targets,
            &reps,
            &contigs,
            HashSet::default(),
            90,
        );
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].cluster, 21);
    }

    #[test]
    fn refind_rejects_an_orf_tied_between_context_families() {
        let mut left = protein(40, 14);
        let mut middle = protein(100, 15);
        let mut right = protein(40, 16);
        for sequence in [&mut left, &mut middle, &mut right] {
            sequence[0] = b'M';
        }
        let left_dna = back_translate(&left);
        let middle_dna = back_translate(&middle);
        let right_dna = back_translate(&right);
        let spacer = "N".repeat(60);
        let right_start = (left_dna.len() + spacer.len() + middle_dna.len() + spacer.len() + 1) as u32;
        let contig = format!("{left_dna}{spacer}{middle_dna}{spacer}{right_dna}");

        let mut pools = Pools::default();
        let calls = vec![
            call(0, 10, 1, &left, &mut pools),
            call(0, 11, right_start, &right, &mut pools),
        ];
        let targets = HashMap::from_iter([(
            flank_key(10, 11),
            vec![20, 21],
        )]);
        let reps = HashMap::from_iter([(20, middle.clone()), (21, middle)]);
        let contigs = vec![vec![PackedSeq::new(contig)]];

        let found = refind_in_sample(
            &calls,
            0,
            &targets,
            &reps,
            &contigs,
            HashSet::default(),
            90,
        );
        assert!(found.is_empty());
    }

    fn reconcile_fixture() -> (Pools, Vec<Call>, Vec<Vec<usize>>) {
        let mut pools = Pools::default();
        let mut full = protein(100, 42);
        full[0] = b'M';
        let mut variant = full.clone();
        variant[40] = if variant[40] == b'A' { b'G' } else { b'A' };
        let mut calls = Vec::new();
        for sample in 0..4 {
            calls.push(call(sample, 10, 100, b"MKTLLYVWAACD", &mut pools));
            calls.push(call(sample, u32::from(sample >= 2), 400,
                if sample < 2 { &full } else { &variant }, &mut pools));
            calls.push(call(sample, 11, 1000, b"MACDEFGHIKLNP", &mut pools));
        }
        (pools, calls, vec![vec![4000]; 4])
    }

    #[test]
    fn same_locus_reconciles_equal_length_alleles_and_preserves_cooccurring_copies() {
        let (mut pools, mut calls, lengths) = reconcile_fixture();
        assert!(merge_fragments(&mut calls, &pools, &lengths));
        let middle: HashSet<_> = calls.iter().filter(|call| call.start == 400)
            .map(|call| call.cluster).collect();
        assert_eq!(middle.len(), 1);
        assert!(!merge_fragments(&mut calls, &pools, &lengths));

        let (_, mut calls, _) = reconcile_fixture();
        let mut duplicate = calls.iter().find(|call| call.sample == 2 && call.start == 400)
            .unwrap().clone();
        duplicate.sample = 0;
        duplicate.start = 2000;
        duplicate.end = 2302;
        calls.push(duplicate);
        assert!(!merge_fragments(&mut calls, &pools, &lengths));
        // An unrelated same-size protein must not pass a context-only merge.
        let unrelated = pools.intern_protein(protein(100, 91));
        assert!(!AlleleIndex::same_locus_sequence(&pools.protein[unrelated as usize],
            &pools.protein[calls[1].protein as usize]));
    }

    #[test]
    fn same_locus_reconciliation_rejects_conflicting_context_and_internal_domains() {
        let (pools, mut calls, lengths) = reconcile_fixture();
        for call in &mut calls {
            if call.sample >= 2 && call.cluster == 11 {
                call.cluster = 12;
            }
        }
        assert!(!merge_fragments(&mut calls, &pools, &lengths));
        let full = protein(100, 42);
        assert!(!AlleleIndex::same_locus_sequence(&full[20..60], &full));
        assert!(AlleleIndex::same_locus_sequence(&full[..40], &full));
    }

    #[test]
    fn shallow_overlaps_retain_neighborhood_anchors() {
        let (pools, mut calls, _) = reconcile_fixture();
        for sample in 0..4 {
            let middle = calls.iter().find(|call| call.sample == sample && call.start == 400)
                .unwrap().end;
            let right = calls.iter_mut().find(|call| call.sample == sample && call.cluster == 11)
                .unwrap();
            right.start = middle - 19;
            right.end = right.start + 41;
        }
        sort_calls(&mut calls);
        let loci = physical_loci(&calls, &pools);
        let contexts = locus_contexts(&calls, &loci, |call| call.cluster);
        for (locus, context) in loci.iter().zip(contexts) {
            if calls[locus.first].start == 400 {
                assert_eq!(context, (Some(10), Some(11)));
            }
        }
    }

    fn merge(calls: &mut [Call], pools: &Pools) -> bool {
        let mut lengths = Vec::new();
        for call in calls.iter() {
            lengths.resize_with(lengths.len().max(call.sample as usize + 1), Vec::new);
            let sample = &mut lengths[call.sample as usize];
            sample.resize(sample.len().max(call.contig as usize + 1), 0);
            sample[call.contig as usize] =
                sample[call.contig as usize].max(call.end as usize + 300);
        }
        merge_fragments(calls, pools, &lengths)
    }

    fn target(
        reps: &HashMap<u32, Vec<u8>>,
        members: &HashMap<u32, HashSet<u32>>,
        flanks: &HashMap<u32, (Option<u32>, Option<u32>)>,
    ) -> Option<MergeTarget> {
        let index = MergeIndex::new(flanks, &HashSet::default());
        best_merge_target(
            0,
            reps,
            members,
            flanks,
            &kmer_sets(reps),
            &index.candidates(0, flanks),
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
        let members = HashMap::from_iter([(0, HashSet::from_iter([0])), (1, HashSet::from_iter([1]))]);
        let flanks = HashMap::from_iter([(0, (Some(10), Some(11))), (1, (Some(10), Some(11)))]);
        for length in [81, 90, 94, 95, 99] {
            for short in [&full[..length], &full[100 - length..]] {
                let reps = HashMap::from_iter([(0, short.to_vec()), (1, full.clone())]);
                let chosen = target(&reps, &members, &flanks).expect("terminal allele");
                assert_eq!(chosen.cluster, 1);
                assert!(chosen.terminal, "length {length}");
            }
        }
    }

    #[test]
    fn terminal_merges_require_both_known_flanks_and_disjoint_samples() {
        let full = protein(100, 42);
        let reps = HashMap::from_iter([(0, full[..90].to_vec()), (1, full)]);
        let mut members =
            HashMap::from_iter([(0, HashSet::from_iter([0])), (1, HashSet::from_iter([1, 2, 3, 4, 5]))]);
        for context in [(Some(10), None), (None, Some(11)), (Some(10), Some(12))] {
            let flanks = HashMap::from_iter([(0, context), (1, (Some(10), Some(11)))]);
            assert!(target(&reps, &members, &flanks).is_none());
        }
        let flanks = HashMap::from_iter([(0, (Some(10), Some(11))), (1, (Some(10), Some(11)))]);
        members.get_mut(&1).unwrap().insert(0);
        assert!(target(&reps, &members, &flanks).is_none());
    }

    #[test]
    fn interior_domains_and_internal_deletions_are_not_terminal_alleles() {
        let full = protein(100, 42);
        let members = HashMap::from_iter([(0, HashSet::from_iter([0])), (1, HashSet::from_iter([1]))]);
        let flanks = HashMap::from_iter([(0, (Some(10), Some(11))), (1, (Some(10), Some(11)))]);
        let mut deleted = full[..45].to_vec();
        deleted.extend_from_slice(&full[55..]);
        for short in [full[5..95].to_vec(), deleted] {
            let reps = HashMap::from_iter([(0, short), (1, full.clone())]);
            assert!(target(&reps, &members, &flanks).is_none());
        }
    }

    #[test]
    fn terminal_identity_requires_full_coverage_and_98_percent_identity() {
        let full = protein(100, 42);
        let members = HashMap::from_iter([(0, HashSet::from_iter([0])), (1, HashSet::from_iter([1]))]);
        let flanks = HashMap::from_iter([(0, (Some(10), Some(11))), (1, (Some(10), Some(11)))]);
        let mut short = full[..90].to_vec();
        for changes in 1..=2 {
            let at = changes * 20;
            short[at] = if short[at] == b'A' { b'C' } else { b'A' };
            let reps = HashMap::from_iter([(0, short.clone()), (1, full.clone())]);
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
        let reps = HashMap::from_iter([
            (0, full[..90].to_vec()),
            (1, full.clone()),
            (2, full[..95].to_vec()),
        ]);
        let members = HashMap::from_iter([
            (0, HashSet::from_iter([0])),
            (1, HashSet::from_iter([1])),
            (2, HashSet::from_iter([2])),
        ]);
        let flanks = HashMap::from_iter([
            (0, (Some(10), Some(11))),
            (1, (Some(10), Some(11))),
            (2, (Some(10), Some(11))),
        ]);
        assert!(target(&reps, &members, &flanks).is_none());
    }

    #[test]
    fn competing_legacy_fragment_targets_are_not_chosen_by_sample_support() {
        let full = protein(120, 42);
        let reps = HashMap::from_iter([
            (0, full[..60].to_vec()),
            (1, full.clone()),
            (2, full[..110].to_vec()),
        ]);
        let members = HashMap::from_iter([
            (0, HashSet::from_iter([0])),
            (1, HashSet::from_iter([1, 2, 3])),
            (2, HashSet::from_iter([4])),
        ]);
        let flanks = HashMap::from_iter([
            (0, (Some(10), None)),
            (1, (Some(10), Some(11))),
            (2, (Some(10), Some(12))),
        ]);
        assert!(target(&reps, &members, &flanks).is_none());
    }

    fn allele_fixture(
        short: &[u8],
        long: &[u8],
        side: AlleleSide,
        reverse: bool,
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
        let lengths = vec![
            vec![if missing_right {
                short_end as usize + 78
            } else {
                4000
            }],
            vec![4000],
            vec![4000],
        ];
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
                    let id = calls
                        .iter()
                        .find(|call| {
                            call.sample == 0
                                && call.protein
                                    == pools
                                        .protein
                                        .iter()
                                        .position(|sequence| sequence == short)
                                        .unwrap() as u32
                        })
                        .unwrap()
                        .cluster;
                    assert!(calls
                        .iter()
                        .filter(|call| call.protein
                            == pools
                                .protein
                                .iter()
                                .position(|sequence| sequence == &full)
                                .unwrap() as u32)
                        .all(|call| call.cluster == id));
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
                let (pools, mut calls, lengths) =
                    allele_fixture(short, &full, AlleleSide::Either, reverse);
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
        assert!(
            allele_target(&mut calls, &pools, &lengths).target.is_none(),
            "91 bases is not near the edge"
        );
        let (pools, mut calls, lengths) =
            allele_fixture(&full[..90], &full, AlleleSide::MissingPrefix, false);
        assert!(
            allele_target(&mut calls, &pools, &lengths).target.is_none(),
            "N-edge evidence cannot justify a C-terminal deletion"
        );
        let (pools, mut calls, lengths) =
            allele_fixture(&full[30..], &full, AlleleSide::MissingPrefix, false);
        // FIX(orientation): edges count the bases outside the call on both
        // sides, so "91 bases is not near the edge" is `start = 92` here, the
        // mirror of `lengths += 13` above (it used to be `start = 91`, i.e.
        // 90 bases, which only the left side treated as interior).
        calls
            .iter_mut()
            .find(|call| call.cluster == 0)
            .unwrap()
            .start = 92;
        assert!(allele_target(&mut calls, &pools, &lengths).target.is_none());
        let (mut pools, mut calls, lengths) =
            allele_fixture(&full[30..], &full, AlleleSide::MissingPrefix, false);
        calls.push(call(0, 99, 1, &intact_protein(80, 123), &mut pools));
        assert!(
            allele_target(&mut calls, &pools, &lengths).target.is_none(),
            "an overlapping predecessor is not a lost contig flank"
        );
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
        let witness = calls
            .iter_mut()
            .find(|call| call.sample == 2 && call.cluster == 1)
            .unwrap();
        witness.protein = shorter;
        witness.dna = dna;
        assert!(
            allele_target(&mut calls, &pools, &lengths).target.is_none(),
            "partial alleles cannot count as full-length witnesses"
        );
        let mut broken = full.clone();
        broken[0] = b'A';
        let (pools, mut calls, lengths) =
            allele_fixture(&full[30..], &broken, AlleleSide::MissingPrefix, false);
        let decision = allele_target(&mut calls, &pools, &lengths);
        assert!(decision.target.is_none());
        assert!(decision.ambiguous);
        assert!(
            !merge_fragments(&mut calls, &pools, &lengths),
            "a nonintact host must not bypass rescue rejection through the legacy route"
        );
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
            bracketed(
                &mut calls,
                &mut pools,
                sample,
                1,
                (10, 11),
                &[(1, full.as_slice())],
            );
        }
        lengths.push(vec![4000]);
        bracketed(&mut calls, &mut pools, 4, 1, (10, 12), &[(0, short)]);
        assert!(allele_target(&mut calls, &pools, &lengths).target.is_none());
        sort_calls(&mut calls);
        let loci = physical_loci(&calls, &pools);
        let raw = locus_contexts(&calls, &loci, |call| call.cluster);
        let stabilized = stabilize_contexts(&calls, &loci, &raw, |call| call.cluster);
        let decision =
            AlleleIndex::new(&calls, &pools, &loci, &stabilized, &lengths).target(0, &pools);
        assert_eq!(decision.target, Some(1));
        assert!(merge_fragments(&mut calls, &pools, &lengths));
        let short_protein = pools
            .protein
            .iter()
            .position(|sequence| sequence == short)
            .unwrap() as u32;
        let full_protein = pools
            .protein
            .iter()
            .position(|sequence| sequence == &full)
            .unwrap() as u32;
        let short_cluster = calls
            .iter()
            .find(|call| call.protein == short_protein)
            .unwrap()
            .cluster;
        let full_cluster = calls
            .iter()
            .find(|call| call.protein == full_protein)
            .unwrap()
            .cluster;
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
            bracketed(
                &mut calls,
                &mut pools,
                3,
                1,
                if change_context { (10, 12) } else { (10, 11) },
                &[(0, &variant)],
            );
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
            let model: Vec<_> = calls
                .iter()
                .filter(|call| call.sample == 1)
                .cloned()
                .collect();
            for mut copy in model {
                copy.sample = sample;
                if copy.cluster == 1 {
                    copy.cluster = 2;
                }
                calls.push(copy);
            }
        }
        bracketed(
            &mut calls,
            &mut pools,
            3,
            3001,
            (10, 11),
            &[(0, &full[..90])],
        );
        let decision = allele_target(&mut calls, &pools, &lengths);
        assert!(decision.ambiguous);
        assert!(decision.target.is_none());
        assert!(!merge_fragments(&mut calls, &pools, &lengths));
        assert_eq!(
            calls
                .iter()
                .filter(|call| call.protein
                    == pools
                        .protein
                        .iter()
                        .position(|sequence| sequence == &full[..90])
                        .unwrap() as u32)
                .map(|call| call.cluster)
                .collect::<HashSet<_>>(),
            HashSet::from_iter([0])
        );
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
            let upstream = pools.protein[calls
                .iter()
                .find(|call| call.cluster == 10)
                .unwrap()
                .protein as usize]
                .clone();
            let downstream = pools.protein[calls
                .iter()
                .find(|call| call.cluster == 11)
                .unwrap()
                .protein as usize]
                .clone();
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
        let model: Vec<_> = calls
            .iter()
            .filter(|call| call.sample == 1)
            .cloned()
            .collect();
        for mut copy in model {
            copy.sample = 0;
            copy.start += 2000;
            copy.end += 2000;
            calls.push(copy);
        }
        assert_eq!(allele_target(&mut calls, &pools, &lengths).target, Some(1));
        assert!(!merge_fragments(&mut calls, &pools, &lengths));
        let short = calls
            .iter()
            .find(|call| {
                call.protein
                    == pools
                        .protein
                        .iter()
                        .position(|sequence| sequence == &full[..90])
                        .unwrap() as u32
            })
            .unwrap();
        assert_eq!(short.cluster, 0);
        // The source group is not combined even though its only sequence host is clear.
        assert!(calls
            .iter()
            .filter(|call| call.protein
                == pools
                    .protein
                    .iter()
                    .position(|sequence| sequence == &full)
                    .unwrap() as u32)
            .all(|call| call.cluster == 1));
    }

    #[test]
    fn shortened_matching_keeps_full_query_identity_and_rejects_domains_and_unknowns() {
        let full = intact_protein(120, 42);
        assert_eq!(
            shortened_identity(&full[..60], &full, AlleleSide::Either),
            Some(1.0)
        );
        assert!(shortened_identity(&full[..59], &full, AlleleSide::Either).is_none());
        assert!(shortened_identity(&full, &full, AlleleSide::Either).is_none());
        assert!(shortened_identity(&full[..29], &full[..40], AlleleSide::Either).is_none());
        let mut unknown = full[..90].to_vec();
        unknown[10] = b'X';
        assert!(shortened_identity(&unknown, &full, AlleleSide::Either).is_none());
        unknown[10] = b'*';
        assert!(shortened_identity(&unknown, &full, AlleleSide::Either).is_none());
        let interior = &full[10..100];
        assert!(
            shortened_identity(interior, &full, AlleleSide::Either).unwrap() < TERMINAL_IDENTITY
        );
        let mut deleted = full[..45].to_vec();
        deleted.extend_from_slice(&full[75..]);
        assert!(
            shortened_identity(&deleted, &full, AlleleSide::Either).unwrap() < TERMINAL_IDENTITY
        );
        let mut variant = full[..100].to_vec();
        for index in [10, 30] {
            variant[index] = if variant[index] == b'A' { b'C' } else { b'A' };
        }
        assert_eq!(
            shortened_identity(&variant, &full, AlleleSide::Either),
            Some(0.98)
        );
        variant[50] = if variant[50] == b'A' { b'C' } else { b'A' };
        assert!(
            shortened_identity(&variant, &full, AlleleSide::Either).unwrap() < TERMINAL_IDENTITY
        );
    }

    fn cluster_of(calls: &[Call], pools: &Pools, sequence: &[u8]) -> u32 {
        let protein = pools
            .protein
            .iter()
            .position(|stored| stored == sequence)
            .unwrap() as u32;
        calls
            .iter()
            .find(|call| call.protein == protein)
            .unwrap()
            .cluster
    }

    #[test]
    fn a_short_call_inside_the_longer_gene_joins_it() {
        let full = intact_protein(120, 42);
        let piece = &full[30..70];
        let mut pools = Pools::default();
        let mut calls = Vec::new();
        for sample in 0..2u32 {
            let long = call(sample, 1, 1000, &full, &mut pools);
            let mut inner = call(sample, 0, long.start + 60, piece, &mut pools);
            assert!(inner.end <= long.end);
            inner.strand = long.strand;
            calls.push(long);
            calls.push(inner);
        }
        let lengths = vec![vec![8000], vec![8000]];
        assert!(merge_fragments(&mut calls, &pools, &lengths));
        assert_eq!(
            cluster_of(&calls, &pools, piece),
            cluster_of(&calls, &pools, &full)
        );
    }

    #[test]
    fn a_short_call_beside_the_longer_gene_stays_separate() {
        let full = intact_protein(120, 42);
        let piece = &full[30..70];
        let mut pools = Pools::default();
        let mut calls = Vec::new();
        for sample in 0..2u32 {
            let long = call(sample, 1, 1000, &full, &mut pools);
            let beside = call(sample, 0, long.end + 50, piece, &mut pools);
            calls.push(long);
            calls.push(beside);
        }
        let lengths = vec![vec![8000], vec![8000]];
        assert!(!merge_fragments(&mut calls, &pools, &lengths));
        assert_ne!(
            cluster_of(&calls, &pools, piece),
            cluster_of(&calls, &pools, &full)
        );
    }

    #[test]
    fn a_short_call_that_never_shares_a_genome_joins_the_longer_gene() {
        let full = intact_protein(120, 42);
        let piece = &full[30..70];
        let mut pools = Pools::default();
        let mut calls = vec![
            call(0, 0, 100, piece, &mut pools),
            call(1, 1, 1000, &full, &mut pools),
        ];
        let lengths = vec![vec![8000], vec![8000]];
        assert!(merge_fragments(&mut calls, &pools, &lengths));
        assert_eq!(
            cluster_of(&calls, &pools, piece),
            cluster_of(&calls, &pools, &full)
        );
    }

    #[test]
    fn two_longer_genes_that_both_contain_the_short_call_are_not_chosen() {
        let full = intact_protein(120, 42);
        let piece = &full[30..70];
        let mut pools = Pools::default();
        let mut calls = vec![
            call(0, 0, 100, piece, &mut pools),
            call(1, 1, 1000, &full, &mut pools),
            call(2, 2, 1000, &full, &mut pools),
        ];
        let lengths = vec![vec![8000], vec![8000], vec![8000]];
        assert!(!merge_fragments(&mut calls, &pools, &lengths));
        assert_eq!(cluster_of(&calls, &pools, piece), 0);
    }

    #[test]
    fn legacy_fragment_merges_keep_one_flank_and_containment_rules() {
        let full = protein(100, 42);
        let members = HashMap::from_iter([(0, HashSet::from_iter([0])), (1, HashSet::from_iter([1]))]);
        let flanks = HashMap::from_iter([(0, (Some(10), None)), (1, (Some(10), Some(11)))]);
        for short in [full[..80].to_vec(), full[20..80].to_vec()] {
            let reps = HashMap::from_iter([(0, short), (1, full.clone())]);
            let chosen = target(&reps, &members, &flanks).expect("legacy fragment");
            assert_eq!(chosen.cluster, 1);
            assert!(!chosen.terminal);
        }
    }

    fn call(sample: u32, cluster: u32, position: u32, sequence: &[u8], pools: &mut Pools) -> Call {
        let dna = crate::genetic_code::back_translate(sequence);
        Call {
            sample,
            contig: 0,
            start: position,
            end: position + dna.len() as u32 - 1,
            strand: 1,
            dna: pools.intern_dna(dna),
            protein: pools.intern_protein(sequence.to_vec()),
            cluster,
            context: cluster,
            locus: None,
            source: None,
            refound: false,
            fusion_piece: false,
            symbol: Label::default(),
            product: Label::default(),
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
                apply_adoptions(&mut calls, &pools, members, &first, &[(1, 2)], &[])
            } else {
                apply_adoptions(&mut calls, &pools, members, &[], &all, &[])
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
        assert!(
            !decision.ambiguous,
            "one intact witness must preserve eligibility for the original terminal rule"
        );
        assert!(merge(&mut calls, &pools));
        let middle: Vec<_> = calls.iter().filter(|call| call.start == 200).collect();
        assert_eq!(middle[0].cluster, middle[1].cluster);
        assert!(!merge(&mut calls, &pools));
    }

    fn planned_partition(
        records: &[(u32, (Option<u32>, Option<u32>))],
    ) -> (Vec<Vec<usize>>, usize) {
        let mut pools = Pools::default();
        let calls: Vec<_> = records
            .iter()
            .enumerate()
            .map(|(index, &(sample, _))| {
                call(
                    sample,
                    0,
                    1000 * index as u32 + 1,
                    b"ACDEFGHIKL",
                    &mut pools,
                )
            })
            .collect();
        let loci: Vec<_> = calls
            .iter()
            .enumerate()
            .map(|(index, call)| Locus {
                first: index,
                last: index,
                end: call.end,
            })
            .collect();
        let contexts: Vec<_> = records.iter().map(|&(_, context)| context).collect();
        let indexes: Vec<_> = (0..records.len()).collect();
        let plan = paralog_split_groups(&indexes, &calls, &loci, &contexts);
        let assigned: HashSet<_> = plan.groups.iter().flatten().copied().collect();
        let retained = indexes
            .into_iter()
            .filter(|index| !assigned.contains(index))
            .collect();
        let mut groups = plan.groups;
        groups.push(retained);
        groups.sort();
        (groups, plan.unresolved_loci)
    }

    #[test]
    fn flank_variants_of_two_paralogs_produce_two_groups_not_four() {
        let records = [
            (0, (Some(10), Some(11))),
            (0, (Some(20), Some(21))),
            (1, (Some(10), Some(12))),
            (1, (Some(20), Some(22))),
            (2, (Some(10), None)),
        ];
        let (groups, unresolved) = planned_partition(&records);
        assert_eq!(groups, vec![vec![0, 2, 4], vec![1, 3]]);
        assert_eq!(unresolved, 0);
    }

    #[test]
    fn ambiguous_single_flank_is_not_assigned_to_an_arbitrary_paralog() {
        let records = [
            (0, (Some(10), Some(11))),
            (0, (Some(10), Some(12))),
            (1, (Some(10), None)),
            (2, (None, None)),
            (3, (None, None)),
        ];
        let (groups, unresolved) = planned_partition(&records);
        assert_eq!(groups, vec![vec![0], vec![1], vec![2, 3, 4]]);
        assert_eq!(unresolved, 3);
    }

    #[test]
    fn context_bridge_cannot_combine_cooccurring_paralogs() {
        let records = [
            (0, (Some(10), Some(11))),
            (1, (Some(10), Some(12))),
            (0, (Some(20), Some(12))),
        ];
        let (groups, _) = planned_partition(&records);
        assert!(!groups
            .iter()
            .any(|group| group.contains(&0) && group.contains(&2)));
    }

    #[test]
    fn unrelated_single_copy_contexts_are_not_split() {
        let records = [
            (0, (Some(10), Some(11))),
            (1, (Some(20), Some(21))),
            (2, (None, None)),
        ];
        assert_eq!(planned_partition(&records), (vec![vec![0, 1, 2]], 0));
    }

    #[test]
    fn identical_multicopy_contexts_are_reported_as_unresolved() {
        let records = [(0, (Some(10), Some(11))), (0, (Some(10), Some(11)))];
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
        assert_eq!(
            locus_contexts(&calls, &loci, |call| call.cluster)[1],
            (Some(10), Some(11))
        );
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
            apply_adoptions(&mut calls, &pools, members, &[(0, 2), (1, 2)], &[], &[]),
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
                apply_adoptions(&mut calls, &pools, members, &[(0, 1)], &[], &[]),
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
            apply_adoptions(&mut calls, &pools, members, &[(0, 1), (1, 2)], &[], &[]),
            (1, 0, 1),
        );
        assert_eq!(calls[0].cluster, calls[2].cluster);
        assert_ne!(calls[0].cluster, calls[1].cluster);
    }

    fn bracketed(
        calls: &mut Vec<Call>,
        pools: &mut Pools,
        sample: u32,
        position: u32,
        flanks: (u32, u32),
        parts: &[(u32, &[u8])],
    ) {
        calls.push(call(sample, flanks.0, position, b"MKTLLYVWAACD", pools));
        let mut start = position + 100;
        for &(cluster, sequence) in parts {
            let part = call(sample, cluster, start, sequence, pools);
            start = part.end + 1;
            calls.push(part);
        }
        calls.push(call(
            sample,
            flanks.1,
            position + 600,
            b"ACDEFGHIKLMNP",
            pools,
        ));
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
        assert_eq!(
            split_paralogs(&mut calls, &pools),
            2,
            "only the flanks split"
        );
        let loci = physical_loci(&calls, &pools);
        let current = locus_contexts(&calls, &loci, |call| call.cluster);
        let stable = locus_contexts(&calls, &loci, |call| call.context);
        let core_contexts = |contexts: &[(Option<u32>, Option<u32>)]| {
            loci.iter()
                .zip(contexts)
                .filter_map(|(locus, &context)| {
                    (pools.protein[calls[locus.first].protein as usize] == full).then_some(context)
                })
                .collect::<HashSet<_>>()
        };
        assert_eq!(core_contexts(&current).len(), 2);
        assert_eq!(core_contexts(&stable).len(), 1);
        for _ in 0..4 {
            assert_eq!(split_paralogs(&mut calls, &pools), 0);
            compress_clusters(&mut calls, &pools);
        }
        let groups: HashSet<_> = calls
            .iter()
            .filter_map(|call| {
                (pools.protein[call.protein as usize] == full).then_some(call.cluster)
            })
            .collect();
        assert_eq!(
            groups.len(),
            1,
            "identical contexts must remain explicitly unresolved"
        );
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
            assert_eq!(
                loci.iter()
                    .filter(|locus| {
                        pools.protein[calls[locus.first].protein as usize].len() >= 40
                    })
                    .count(),
                4,
                "each fragmented copy is still one locus"
            );
            assert_eq!(split_paralogs(&mut calls, &pools), 0);
        }
        let a = calls
            .iter()
            .find(|call| call.sample == 0 && call.start == 101)
            .unwrap();
        let b = calls
            .iter()
            .find(|call| call.sample == 0 && call.start == 1101)
            .unwrap();
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
        assert_eq!(
            apply_adoptions(&mut calls, &pools, members, &[], &[(10, 11)], &[]),
            (1, 1, 0)
        );
        assert_eq!(calls[0].cluster, calls[1].cluster);
        assert_ne!(
            calls[0].cluster, calls[2].cluster,
            "anchors are not ortholog assignments"
        );
        assert!(calls.iter().all(|call| call.context == 10));
        compress_clusters(&mut calls, &pools);
        assert!(calls.iter().all(|call| call.context == 10));
    }

    #[test]
    fn exact_whole_alleles_bridge_both_changed_flanks_without_combining_copies() {
        let mut pools = Pools::default();
        let a = protein(80, 42);
        let b = protein(80, 123);
        let mut calls = Vec::new();
        for (sample, offset) in [(0, 0), (1, 20)] {
            bracketed(
                &mut calls,
                &mut pools,
                sample,
                1,
                (10 + offset, 11 + offset),
                &[(0, &a)],
            );
            bracketed(
                &mut calls,
                &mut pools,
                sample,
                1001,
                (12 + offset, 13 + offset),
                &[(0, &b)],
            );
        }
        freeze_contexts(&mut calls, &pools);
        assert_eq!(split_paralogs(&mut calls, &pools), 1);
        let cluster_at = |sample, position| {
            calls
                .iter()
                .find(|call| call.sample == sample && call.start == position)
                .unwrap()
                .cluster
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
            bracketed(
                &mut calls,
                &mut pools,
                sample,
                1,
                (10 + offset, 11 + offset),
                &[(0, &full)],
            );
            bracketed(
                &mut calls,
                &mut pools,
                sample,
                1001,
                (12 + offset, 13 + offset),
                &[(0, &full)],
            );
        }
        freeze_contexts(&mut calls, &pools);
        assert_eq!(split_paralogs(&mut calls, &pools), 3);
        let groups: HashSet<_> = calls
            .iter()
            .filter_map(|call| {
                (pools.protein[call.protein as usize] == full).then_some(call.cluster)
            })
            .collect();
        assert_eq!(
            groups.len(),
            4,
            "ambiguous cross-sample copies must not be assigned by order"
        );
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
            bracketed(
                &mut calls,
                &mut pools,
                sample,
                1,
                (10 + offset, 11 + offset),
                &[(0, allele)],
            );
            bracketed(
                &mut calls,
                &mut pools,
                sample,
                1001,
                (12 + offset, 13 + offset),
                &[(0, &b)],
            );
        }
        freeze_contexts(&mut calls, &pools);
        assert_eq!(split_paralogs(&mut calls, &pools), 2);
        let cluster_at = |sample, position| {
            calls
                .iter()
                .find(|call| call.sample == sample && call.start == position)
                .unwrap()
                .cluster
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
        bracketed(
            &mut calls,
            &mut pools,
            0,
            1,
            (10, 11),
            &[(0, head), (0, tail)],
        );
        bracketed(&mut calls, &mut pools, 0, 1001, (20, 21), &[(0, &full)]);
        bracketed(&mut calls, &mut pools, 1, 1, (30, 31), &[(0, head)]);
        freeze_contexts(&mut calls, &pools);
        assert_eq!(split_paralogs(&mut calls, &pools), 2);
        let a = calls
            .iter()
            .find(|call| call.sample == 0 && call.start == 101)
            .unwrap();
        let b = calls
            .iter()
            .find(|call| call.sample == 1 && call.start == 101)
            .unwrap();
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
            let calls: Vec<_> = (0..30)
                .map(|index| {
                    let sample = (next() % 6) as u32;
                    let left = (next() % 5) as u32;
                    let right = (next() % 5) as u32;
                    contexts.push(((left != 0).then_some(left), (right != 0).then_some(right)));
                    let allele = if next() % 2 == 0 {
                        b"ACDEFGHIKL"
                    } else {
                        b"MNPQRSTVWY"
                    };
                    call(sample, 0, index * 1000 + 1, allele, &mut pools)
                })
                .collect();
            let loci: Vec<_> = calls
                .iter()
                .enumerate()
                .map(|(index, call)| Locus {
                    first: index,
                    last: index,
                    end: call.end,
                })
                .collect();
            let indexes: Vec<_> = (0..calls.len()).collect();
            let plan = paralog_split_groups(&indexes, &calls, &loci, &contexts);
            let assigned: HashSet<_> = plan.groups.iter().flatten().copied().collect();
            let retained: Vec<_> = indexes
                .into_iter()
                .filter(|index| !assigned.contains(index))
                .collect();
            let groups = plan.groups.into_iter().chain(std::iter::once(retained));
            for group in groups {
                assert!(
                    paralog_split_groups(&group, &calls, &loci, &contexts)
                        .groups
                        .is_empty(),
                    "seed {seed}"
                );
                for (position, &a) in group.iter().enumerate() {
                    for &b in &group[position + 1..] {
                        if calls[a].sample == calls[b].sample {
                            assert_eq!(
                                contexts[a], contexts[b],
                                "co-occurring contexts joined for seed {seed}"
                            );
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
        compress_clusters(&mut calls, &pools);
        let at_a = calls
            .iter()
            .find(|call| call.sample == 1 && call.start == 101)
            .unwrap()
            .cluster;
        let at_b = calls
            .iter()
            .find(|call| call.sample == 1 && call.start == 1101)
            .unwrap()
            .cluster;
        assert_ne!(at_a, at_b);
        for call in calls
            .iter()
            .filter(|call| call.sample == 2 && (101..500).contains(&call.start))
        {
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
            calls
                .iter()
                .map(|call| (call.sample, call.start, call.cluster))
                .collect::<Vec<_>>()
        };
        assert_eq!(labels(&calls), labels(&reversed));
    }

    #[test]
    fn adjacent_intact_genes_stay_separate() {
        let left = block(b"ACDEFGHIKL", 80);
        let right = block(b"NPQRSTVWY", 90);
        let identity = crate::align::global_identity(&left, &right).expect("comparable proteins");
        assert!(identity < 0.4, "{identity}");
        let mut pools = Pools::default();
        let mut calls = Vec::new();
        for sample in 0..2 {
            let first = call(sample, 1, 100, &left, &mut pools);
            let second = call(sample, 2, first.end + 1, &right, &mut pools);
            calls.push(first);
            calls.push(second);
        }
        let lengths = vec![vec![5000], vec![5000]];
        merge_fragments(&mut calls, &pools, &lengths);
        let clusters: HashSet<u32> = calls.iter().map(|call| call.cluster).collect();
        assert_eq!(clusters.len(), 2, "{clusters:?}");
    }

    #[test]
    fn a_fragment_still_joins_the_one_gene_it_belongs_to() {
        let gene = block(b"ACDEFGHIKL", 80);
        let fragment = unfinished(b"NPQRSTWY", 40);
        let identity = crate::align::global_identity(&gene, &fragment).expect("comparable");
        assert!(identity < 0.4, "{identity}");
        let mut pools = Pools::default();
        let mut calls = Vec::new();
        for sample in 0..2 {
            let first = call(sample, 1, 100, &gene, &mut pools);
            calls.push(call(sample, 2, first.end + 1, &fragment, &mut pools));
            calls.push(first);
        }
        let lengths = vec![vec![5000], vec![5000]];
        assert!(merge_fragments(&mut calls, &pools, &lengths));
        let clusters: HashSet<u32> = calls.iter().map(|call| call.cluster).collect();
        assert_eq!(clusters.len(), 1, "{clusters:?}");
    }

    #[test]
    fn a_fragment_between_two_genes_does_not_bridge_them() {
        let left = block(b"ACDEFGHIKL", 80);
        let middle = unfinished(b"Q", 30);
        let right = block(b"NPRSTWY", 70);
        for (one, other) in [(&left, &middle), (&middle, &right), (&left, &right)] {
            let identity = crate::align::global_identity(one, other).expect("comparable");
            assert!(identity < 0.4, "{identity}");
        }
        let mut pools = Pools::default();
        let mut calls = Vec::new();
        for sample in 0..2 {
            let first = call(sample, 1, 100, &left, &mut pools);
            let mid = call(sample, 2, first.end + 1, &middle, &mut pools);
            calls.push(call(sample, 3, mid.end + 1, &right, &mut pools));
            calls.push(first);
            calls.push(mid);
        }
        let lengths = vec![vec![5000], vec![5000]];
        merge_fragments(&mut calls, &pools, &lengths);
        let cluster = |sequence: &[u8]| cluster_of(&calls, &pools, sequence);
        assert_ne!(cluster(&left), cluster(&right));
    }

    fn unfinished(residues: &[u8], length: usize) -> Vec<u8> {
        let mut sequence = Vec::with_capacity(length);
        for index in 0..length {
            let residue = residues[index % residues.len()];
            // M and V back-translate to start codons, which would make this
            // piece a finished gene and refuse the join.
            sequence.push(if residue == b'M' || residue == b'V' {
                b'A'
            } else {
                residue
            });
        }
        sequence
    }

    fn block(residues: &[u8], length: usize) -> Vec<u8> {
        let mut sequence = Vec::with_capacity(length);
        sequence.push(b'M');
        for index in 1..length {
            let residue = residues[index % residues.len()];
            sequence.push(if residue == b'M' { b'L' } else { residue });
        }
        sequence
    }

    #[test]
    fn union_find_resolves_and_compresses_chains_longer_than_64() {
        let mut parent: HashMap<u32, u32> = (0..200).map(|id| (id, id + 1)).collect();
        assert_eq!(find(&mut parent, 0), 200);
        assert_eq!(parent[&0], 200);
        assert_eq!(find(&mut parent, 1), 200);
    }

    // ---- FIX tests: orientation, fusion pieces, divergent members ----

    fn reversed(calls: &[Call], length: u32) -> Vec<Call> {
        calls
            .iter()
            .map(|call| {
                let mut flipped = call.clone();
                flipped.start = length - call.end + 1;
                flipped.end = length - call.start + 1;
                flipped.strand = -call.strand;
                flipped
            })
            .collect()
    }

    #[test]
    fn fix_contexts_are_the_same_on_a_reversed_contig() {
        let mut pools = Pools::default();
        let mut calls = vec![
            call(0, 10, 1, &intact_protein(60, 1), &mut pools),
            call(0, 20, 400, &intact_protein(70, 2), &mut pools),
            call(0, 11, 900, &intact_protein(80, 3), &mut pools),
        ];
        calls[1].strand = -1;
        let length = 2000;
        let mut flipped = reversed(&calls, length);
        sort_calls(&mut calls);
        sort_calls(&mut flipped);
        let context_of = |calls: &[Call], pools: &Pools| {
            let loci = physical_loci(calls, pools);
            let contexts = locus_contexts(calls, &loci, |call| call.cluster);
            loci.iter()
                .zip(contexts)
                .map(|(locus, context)| (calls[locus.first].cluster, context))
                .collect::<HashMap<u32, (Option<u32>, Option<u32>)>>()
        };
        let forward = context_of(&calls, &pools);
        let backward = context_of(&flipped, &pools);
        assert_eq!(forward, backward);
        // The minus-strand middle gene reads 11 upstream and 10 downstream.
        assert_eq!(forward[&20], (Some(11), Some(10)));
        assert_eq!(forward[&10], (None, Some(20)));
    }

    #[test]
    fn fix_refind_flank_pairs_are_unordered() {
        assert_eq!(flank_key(3, 7), flank_key(7, 3));
        assert_eq!(flank_key(5, 5), (5, 5));
    }

    #[test]
    fn context_graph_ranks_distinct_sample_support_and_bounds_candidates() {
        let mut targets = Vec::new();
        for (family, support) in [(9, 1), (20, 2), (7, 2), (4, 3), (2, 1)]
            .into_iter()
            .chain((10..=20).map(|family| (family, 2)))
        {
            insert_context_target(&mut targets, family, support);
        }
        assert_eq!(
            targets,
            vec![
                (4, 3),
                (7, 2),
                (10, 2),
                (11, 2),
                (12, 2),
                (13, 2),
                (14, 2),
                (15, 2)
            ]
        );
        let mut singleton_context = Vec::new();
        insert_context_target(&mut singleton_context, 30, 1);
        assert_eq!(singleton_context, vec![(30, 1)]);
    }

    #[test]
    fn context_graph_counts_repeated_loci_once_per_sample() {
        let mut support = HashMap::default();
        let key = (3, 7, 12);
        record_context_sample(&mut support, key, 0);
        record_context_sample(&mut support, key, 0);
        assert_eq!(support[&key].0, 1);
        record_context_sample(&mut support, key, 1);
        assert_eq!(support[&key].0, 2);
    }

    #[test]
    fn fix_cluster_numbering_ignores_contig_orientation() {
        let mut pools = Pools::default();
        let x = intact_protein(60, 11);
        let y = intact_protein(60, 12);
        let mut forward = vec![
            call(0, 5, 1, &x, &mut pools),
            call(0, 9, 500, &y, &mut pools),
            call(1, 9, 1, &y, &mut pools),
        ];
        let mut backward = forward.clone();
        let first_sample = reversed(&backward[..2], 1000);
        backward.splice(0..2, first_sample);
        compress_clusters(&mut forward, &pools);
        compress_clusters(&mut backward, &pools);
        let id_of = |calls: &[Call], protein: &[u8]| {
            calls
                .iter()
                .find(|call| pools.protein[call.protein as usize] == protein)
                .unwrap()
                .cluster
        };
        assert_eq!(id_of(&forward, &x), id_of(&backward, &x));
        assert_eq!(id_of(&forward, &y), id_of(&backward, &y));
    }

    #[test]
    fn fix_fusion_pieces_are_genome_bases_on_both_strands() {
        let mut coding = crate::genetic_code::back_translate(&intact_protein(60, 77));
        // A real TGA stop, so a rewritten TAA would be caught.
        coding.replace_range(coding.len() - 3.., "TGA");
        let filler = "GATTACA".repeat(20);
        let cut = 25usize;
        for strand in [1i8, -1] {
            let placed = if strand > 0 {
                coding.clone()
            } else {
                crate::genetic_code::reverse_complement(&coding)
            };
            let contig = format!("{filler}{placed}{filler}");
            let contigs = vec![vec![PackedSeq::new(contig.clone())]];
            let mut pools = Pools::default();
            let mut fused = call(0, 0, filler.len() as u32 + 1, &intact_protein(60, 77), &mut pools);
            fused.strand = strand;
            fused.dna = pools.intern_dna(coding.clone());
            let pieces = fusion_pieces(&fused, &coding, cut, &contigs).expect("exact CDS splits");
            assert_eq!(pieces.prefix, coding[..cut * 3]);
            assert_eq!(pieces.suffix, coding[cut * 3..]);
            assert!(pieces.suffix.ends_with("TGA"), "the real stop is kept");
            assert_eq!(pieces.prefix.len() + pieces.suffix.len(), coding.len(), "nothing is invented");
            let genome = |(start, end): (u32, u32)| {
                let bases = &contig[start as usize - 1..end as usize];
                if strand > 0 {
                    bases.to_string()
                } else {
                    crate::genetic_code::reverse_complement(bases)
                }
            };
            assert_eq!(genome(pieces.prefix_span), pieces.prefix);
            assert_eq!(genome(pieces.suffix_span), pieces.suffix);
            let (first, second) = if strand > 0 {
                (pieces.prefix_span, pieces.suffix_span)
            } else {
                (pieces.suffix_span, pieces.prefix_span)
            };
            assert_eq!(first.0, fused.start);
            assert_eq!(first.1 + 1, second.0, "pieces abut with no overlap");
            assert_eq!(second.1, fused.end);
            // Bases that are not the genome's, or a span that is not the CDS,
            // are left whole.
            let mut wrong = coding.clone();
            wrong.replace_range(0..3, "ATA");
            assert!(fusion_pieces(&fused, &wrong, cut, &contigs).is_none());
            let mut longer = fused.clone();
            longer.end += 3;
            assert!(fusion_pieces(&longer, &coding, cut, &contigs).is_none());
            assert!(fusion_pieces(&fused, &coding, 0, &contigs).is_none());
        }
    }

    #[test]
    fn fix_allele_identity_is_one_for_equal_and_low_for_unrelated() {
        let a = intact_protein(200, 5);
        assert_eq!(allele_identity(&a, &a), 1.0);
        assert!(allele_identity(&a, &intact_protein(195, 6)) < 0.3);
        let mut close = a.clone();
        for index in (10..200).step_by(10) {
            close[index] = if close[index] == b'W' { b'Y' } else { b'W' };
        }
        let identity = allele_identity(&a, &close);
        assert!(identity > 0.85 && identity < 1.0, "{identity}");
        // A fragment is scored over its own length.
        assert!(allele_identity(&a[50..120], &a) > 0.99);
    }

    /// Family 0: P between flanks 10 and 11 in samples 0-3, a divergent Q
    /// between other flanks (12, 13) in samples 4-5, and a divergent Q2 in the
    /// family's own neighbourhood in sample 6.
    fn divergent_fixture() -> (Pools, Vec<Call>, Vec<u8>, Vec<u8>, Vec<u8>) {
        let mut pools = Pools::default();
        let p = intact_protein(200, 21);
        let q = intact_protein(190, 22);
        let q2 = intact_protein(195, 23);
        let mut calls = Vec::new();
        let mut triplet = |sample: u32, left: u32, middle: &[u8], right: u32, calls: &mut Vec<Call>| {
            calls.push(call(sample, left, 1, &intact_protein(80, u64::from(left)), &mut pools));
            calls.push(call(sample, 0, 400, middle, &mut pools));
            calls.push(call(sample, right, 1200, &intact_protein(90, u64::from(right)), &mut pools));
        };
        for sample in 0..4 {
            triplet(sample, 10, &p, 11, &mut calls);
        }
        for sample in 4..6 {
            triplet(sample, 12, &q, 13, &mut calls);
        }
        triplet(6, 10, &q2, 11, &mut calls);
        (pools, calls, p, q, q2)
    }

    #[test]
    fn fix_divergent_members_in_another_neighbourhood_are_split() {
        let (pools, mut calls, p, q, q2) = divergent_fixture();
        let result = split_divergent_members(&mut calls, &pools, 0.70);
        assert_eq!(result.new_families, 1);
        assert_eq!(result.families, 1);
        assert_eq!(result.loci, 2);
        assert_eq!(result.kept_alleles, 1, "Q2 sits between the family's own flanks");
        let cluster_of = |protein: &[u8]| -> HashSet<u32> {
            calls
                .iter()
                .filter(|call| pools.protein[call.protein as usize] == protein)
                .map(|call| call.cluster)
                .collect()
        };
        let (p_ids, q_ids, q2_ids) = (cluster_of(&p), cluster_of(&q), cluster_of(&q2));
        assert_eq!(q_ids.len(), 1, "both Q copies form one new family");
        assert!(p_ids.is_disjoint(&q_ids));
        assert_eq!(p_ids, q2_ids);
        // Disabled, nothing moves.
        let (pools, mut calls, ..) = divergent_fixture();
        let before: Vec<u32> = calls.iter().map(|call| call.cluster).collect();
        assert_eq!(split_divergent_members(&mut calls, &pools, 0.0), DivergentSplit::default());
        assert_eq!(before, calls.iter().map(|call| call.cluster).collect::<Vec<_>>());
    }

    #[test]
    fn fix_divergent_copy_sharing_a_genome_is_split_even_beside_the_same_flanks() {
        let (mut pools, mut calls, _, _, q2) = divergent_fixture();
        // A second, divergent copy in sample 0, far down the same contig.
        calls.push(call(0, 10, 5000, &intact_protein(80, 10), &mut pools));
        calls.push(call(0, 0, 5400, &q2, &mut pools));
        calls.push(call(0, 11, 6200, &intact_protein(90, 11), &mut pools));
        let result = split_divergent_members(&mut calls, &pools, 0.70);
        // Q (other neighbourhood) and Q2 (now co-occurring with P) both leave.
        assert_eq!(result.new_families, 2);
        assert_eq!(result.kept_alleles, 0);
    }

    #[test]
    fn fix_divergent_members_without_neighbours_stay() {
        let mut pools = Pools::default();
        let mut calls = vec![
            call(0, 0, 400, &intact_protein(200, 31), &mut pools),
            call(1, 0, 400, &intact_protein(200, 31), &mut pools),
            call(2, 0, 400, &intact_protein(190, 32), &mut pools),
        ];
        let result = split_divergent_members(&mut calls, &pools, 0.70);
        assert_eq!(result.new_families, 0);
        assert_eq!(result.kept_unplaced, 1);
    }
}
