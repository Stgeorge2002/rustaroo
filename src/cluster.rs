//! Two-step homology.
//!
//! Identical proteins collapse by hash. What remains is clustered longest-first,
//! the same order CD-HIT uses, so a fragment cannot become the representative
//! of a full gene. The first pass is 98% identity. Families are a second pass
//! on those representatives at the family threshold. Length ratio stays strict
//! here; a half-length fragment is the synteny step's problem, not this one's.

use std::cell::RefCell;
use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

use rayon::prelude::*;

use crate::align::{global_identity, length_ratio};
use crate::log::progress;

// Both caps used to bound the cost of the one thing that was expensive:
// every extra candidate meant another O(n*m) full alignment. Each candidate
// now costs O(n*band) in `global_identity`, typically several times cheaper
// for the near-equal-length pairs that reach this point, so both caps can
// afford to be wider - catching homologs that previously fell outside the
// cutoff in large, highly conserved gene families - for a fraction of the
// old per-candidate cost.
const MAX_CANDIDATES: usize = 160;
const SATURATED_POSTING: usize = 8_000;
// Below this many candidates, spreading the (now cheap) alignment work
// across threads costs more in dispatch overhead than it saves. Above it,
// scoring candidates in parallel is worth doing even though the decision
// that follows - which candidate wins - stays strictly sequential.
const PARALLEL_CANDIDATE_THRESHOLD: usize = 16;

pub fn cluster_families(
    proteins: &[Vec<u8>],
    precluster_identity: f64,
    family_identity: f64,
    length_ratio_min: f64,
) -> Vec<u32> {
    cluster_families_scheduled(
        proteins,
        precluster_identity,
        family_identity,
        length_ratio_min,
        true,
    )
}

/// `chunked` scores several still-unordered proteins against one frozen
/// representative index. A protein that becomes a new representative ends
/// the chunk: everything after it has to see that representative, which is
/// the same rule as walking the proteins one at a time.
fn cluster_families_scheduled(
    proteins: &[Vec<u8>],
    precluster_identity: f64,
    family_identity: f64,
    length_ratio_min: f64,
    chunked: bool,
) -> Vec<u32> {
    if proteins.is_empty() {
        return Vec::new();
    }
    let (unique, back) = unique_proteins(proteins);
    let (pre, pre_reps) = greedy(
        &unique,
        proteins,
        precluster_identity,
        length_ratio_min,
        chunked,
    );
    let rep_indices: Vec<usize> = pre_reps.iter().map(|&slot| unique[slot]).collect();
    let (families, _) = if rep_indices.len() <= 1 {
        (vec![0u32; rep_indices.len()], Vec::new())
    } else {
        greedy(
            &rep_indices,
            proteins,
            family_identity,
            length_ratio_min,
            chunked,
        )
    };
    back.into_iter()
        .map(|unique_slot| {
            let precluster = pre[unique_slot] as usize;
            families[precluster]
        })
        .collect()
}

/// Distinct proteins by content, kept as indices into `proteins` rather than
/// owned copies. `unique[slot]` is the first original index with that exact
/// sequence; `back[i]` is the slot `proteins[i]` dedupes to. This used to
/// clone every unique protein's bytes twice (once into a `Vec<Vec<u8>>`, once
/// again as the hash map key); borrowing `&[u8]` slices of `proteins` that
/// already own the data removes both copies.
fn unique_proteins(proteins: &[Vec<u8>]) -> (Vec<usize>, Vec<usize>) {
    let mut unique = Vec::new();
    let mut index_of: HashMap<&[u8], usize> = HashMap::new();
    let mut back = Vec::with_capacity(proteins.len());
    for (index, protein) in proteins.iter().enumerate() {
        if let Some(&slot) = index_of.get(protein.as_slice()) {
            back.push(slot);
            continue;
        }
        let slot = unique.len();
        unique.push(index);
        index_of.insert(protein.as_slice(), slot);
        back.push(slot);
    }
    (unique, back)
}

/// Greedy longest-first clustering over `indices`, a set of positions into
/// `proteins`. Returns, per local position in `indices`, its cluster id, and
/// the local positions that became representatives.
fn greedy(
    indices: &[usize],
    proteins: &[Vec<u8>],
    identity: f64,
    length_ratio_min: f64,
    chunked: bool,
) -> (Vec<u32>, Vec<usize>) {
    let mut order: Vec<usize> = (0..indices.len()).collect();
    // The comparator is a total order, so a parallel sort cannot reorder a tie.
    order.par_sort_by(|&a, &b| {
        let pa = &proteins[indices[a]];
        let pb = &proteins[indices[b]];
        pb.len()
            .cmp(&pa.len())
            .then_with(|| pa.cmp(pb))
            .then(a.cmp(&b))
    });
    // Word size 5 is right for a 98% precluster. At 70% identity those exact
    // 5-mers disappear, so the family pass seeds with 3-mers instead.
    let kmer_size = if identity >= 0.90 { 5 } else { 3 };
    // K-mer extraction for one protein never depends on another protein or
    // on any clustering decision, so unlike the decision loop below, this is
    // embarrassingly parallel: compute every protein's k-mers up front
    // across all cores instead of one at a time, interleaved, inside the
    // sequential loop that actually has to stay ordered.
    let kmer_lists: Vec<Vec<u64>> = (0..indices.len())
        .into_par_iter()
        .map(|pos| kmers(&proteins[indices[pos]], kmer_size))
        .collect();
    let mut assignment = vec![0u32; indices.len()];
    let mut reps: Vec<usize> = Vec::new();
    let mut postings = Postings::new(kmer_size);
    let total = order.len();
    let step = (total / 20).max(1);
    if total >= 100 {
        progress!("rustroo: {identity:.2} identity pass, {total} proteins, word size {kmer_size}");
    }
    let mut cursor = 0usize;
    let mut width = if chunked { 32usize } else { 1 };
    while cursor < total {
        let end = if chunked {
            (cursor + width).min(total)
        } else {
            cursor + 1
        };
        // Every protein in this span is scored against the representatives
        // that already exist. The span is thrown away at the first protein
        // that founds a new representative, because a later protein in the
        // same span might have joined it.
        let chosen: Vec<Option<u32>> = {
            let reps = reps.as_slice();
            if end - cursor == 1 {
                let pos = order[cursor];
                vec![choose_representative(
                    &proteins[indices[pos]],
                    &kmer_lists[pos],
                    reps,
                    indices,
                    proteins,
                    &postings,
                    identity,
                    length_ratio_min,
                    true,
                )]
            } else {
                order[cursor..end]
                    .par_iter()
                    .map(|&pos| {
                        choose_representative(
                            &proteins[indices[pos]],
                            &kmer_lists[pos],
                            reps,
                            indices,
                            proteins,
                            &postings,
                            identity,
                            length_ratio_min,
                            false,
                        )
                    })
                    .collect()
            }
        };
        let chunk_start = cursor;
        for choice in chosen {
            let pos = order[cursor];
            let founded = if let Some(rep_slot) = choice {
                assignment[pos] = assignment[reps[rep_slot as usize]];
                false
            } else {
                let rep_slot = reps.len() as u32;
                assignment[pos] = rep_slot;
                reps.push(pos);
                for &kmer in &kmer_lists[pos] {
                    postings.push(kmer, rep_slot);
                }
                true
            };
            cursor += 1;
            if total >= 100 && (cursor == total || cursor % step == 0) {
                progress!(
                    "rustroo: {identity:.2} identity {cursor}/{total} proteins, {} clusters",
                    reps.len()
                );
            }
            if founded {
                break;
            }
        }
        if !chunked {
            continue;
        }
        if cursor == end {
            width = (width * 2).min(512);
        } else if cursor == chunk_start + 1 {
            width = (width / 2).max(4);
        }
    }
    (assignment, reps)
}

fn choose_representative(
    protein: &[u8],
    protein_kmers: &[u64],
    reps: &[usize],
    indices: &[usize],
    proteins: &[Vec<u8>],
    postings: &Postings,
    identity: f64,
    length_ratio_min: f64,
    parallel_scores: bool,
) -> Option<u32> {
    let candidates = candidates_for(
        protein,
        reps,
        indices,
        proteins,
        postings,
        length_ratio_min,
        protein_kmers,
    );
    let scored = score_candidates(
        protein,
        &candidates,
        reps,
        indices,
        proteins,
        parallel_scores,
    );
    let mut best: Option<(u32, f64)> = None;
    for (rep_slot, id) in scored.into_iter().flatten() {
        if id + f64::EPSILON < identity {
            continue;
        }
        let replace = match best {
            None => true,
            Some((_, best_id)) => id > best_id + 1e-12,
        };
        if replace {
            best = Some((rep_slot, id));
        }
    }
    best.map(|(rep_slot, _)| rep_slot)
}

fn candidates_for(
    protein: &[u8],
    reps: &[usize],
    indices: &[usize],
    proteins: &[Vec<u8>],
    postings: &Postings,
    length_ratio_min: f64,
    protein_kmers: &[u64],
) -> Vec<u32> {
    if reps.is_empty() {
        return Vec::new();
    }
    HITS.with(|cell| {
        let mut hits = cell.borrow_mut();
        hits.begin(reps.len());
        for kmer in protein_kmers {
            let Some(list) = postings.get(*kmer) else {
                continue;
            };
            if list.len() >= SATURATED_POSTING {
                continue;
            }
            for &rep in list {
                hits.bump(rep);
            }
        }
        let mut ranked = Vec::new();
        for index in 0..hits.touched.len() {
            let rep = hits.touched[index];
            let count = hits.count[rep as usize];
            if count < 2 {
                continue;
            }
            let rep_protein = &proteins[indices[reps[rep as usize]]];
            if length_ratio(protein, rep_protein) + f64::EPSILON >= length_ratio_min {
                ranked.push((count, rep));
            }
        }
        ranked.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
        ranked.truncate(MAX_CANDIDATES);
        ranked.into_iter().map(|(_, rep)| rep).collect()
    })
}

/// Dense stamp table instead of a `HashMap` per query. `generation` lets the
/// table be reused without clearing every representative slot.
struct HitScratch {
    stamp: Vec<u32>,
    count: Vec<u16>,
    generation: u32,
    touched: Vec<u32>,
}

impl HitScratch {
    fn begin(&mut self, n_reps: usize) {
        if self.stamp.len() < n_reps {
            self.stamp.resize(n_reps, 0);
            self.count.resize(n_reps, 0);
        }
        self.generation = self.generation.wrapping_add(1);
        if self.generation == 0 {
            self.stamp.fill(0);
            self.generation = 1;
        }
        self.touched.clear();
    }

    fn bump(&mut self, rep: u32) {
        let index = rep as usize;
        if self.stamp[index] != self.generation {
            self.stamp[index] = self.generation;
            self.count[index] = 1;
            self.touched.push(rep);
        } else {
            self.count[index] += 1;
        }
    }
}

thread_local! {
    static HITS: RefCell<HitScratch> = RefCell::new(HitScratch {
        stamp: Vec::new(),
        count: Vec::new(),
        generation: 0,
        touched: Vec::new(),
    });
}

/// Word index. A 3-mer fits in 15 bits, so the family pass uses a flat table
/// and skips hashing. A 5-mer does not, and those keys are already mixed by
/// `KmerHasher` instead of the default hasher.
struct Postings {
    dense: Vec<Vec<u32>>,
    sparse: HashMap<u64, Vec<u32>, BuildHasherDefault<KmerHasher>>,
}

impl Postings {
    fn new(kmer_size: usize) -> Self {
        let dense = if kmer_size <= 3 {
            vec![Vec::new(); 1 << (5 * kmer_size)]
        } else {
            Vec::new()
        };
        Self {
            dense,
            sparse: HashMap::with_hasher(BuildHasherDefault::default()),
        }
    }

    fn get(&self, kmer: u64) -> Option<&[u32]> {
        if self.dense.is_empty() {
            self.sparse.get(&kmer).map(Vec::as_slice)
        } else {
            self.dense.get(kmer as usize).map(Vec::as_slice)
        }
    }

    fn push(&mut self, kmer: u64, rep: u32) {
        let list = if self.dense.is_empty() {
            self.sparse.entry(kmer).or_default()
        } else {
            &mut self.dense[kmer as usize]
        };
        if list.len() < SATURATED_POSTING {
            list.push(rep);
        }
    }
}

#[derive(Clone, Default)]
struct KmerHasher(u64);

impl Hasher for KmerHasher {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        for &byte in bytes {
            self.write_u64(u64::from(byte));
        }
    }

    fn write_u64(&mut self, value: u64) {
        self.0 = (self.0.rotate_left(5) ^ value).wrapping_mul(0x517cc1b727220a95);
    }
}

/// Identity of `protein` against every candidate representative, in the same
/// order as `candidates`. Above `PARALLEL_CANDIDATE_THRESHOLD` this runs
/// across all cores: each alignment only reads `reps`/`proteins`, so doing
/// them concurrently cannot change any result, only how long it takes. The
/// caller still walks this output strictly in order to decide a winner,
/// which is what keeps a near-tie resolving the same way every run.
fn score_candidates(
    protein: &[u8],
    candidates: &[u32],
    reps: &[usize],
    indices: &[usize],
    proteins: &[Vec<u8>],
    parallel: bool,
) -> Vec<Option<(u32, f64)>> {
    let score_one = |&rep_slot: &u32| -> Option<(u32, f64)> {
        let rep_index = indices[reps[rep_slot as usize]];
        global_identity(protein, &proteins[rep_index]).map(|id| (rep_slot, id))
    };
    // A query that is already one task of a parallel chunk scores its own
    // candidates here. Nesting another parallel iterator on a short list
    // only pays the scheduler.
    if parallel && candidates.len() >= PARALLEL_CANDIDATE_THRESHOLD {
        candidates.par_iter().map(score_one).collect()
    } else {
        candidates.iter().map(score_one).collect()
    }
}

fn kmers(sequence: &[u8], k: usize) -> Vec<u64> {
    if sequence.len() < k {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(sequence.len() - k + 1);
    for i in 0..=sequence.len() - k {
        let window = &sequence[i..i + k];
        if let Some(packed) = pack(window) {
            out.push(packed);
        }
    }
    out
}

fn pack(window: &[u8]) -> Option<u64> {
    let mut distinct = [false; 20];
    let mut n_distinct = 0u8;
    let mut value = 0u64;
    for &residue in window {
        let code = match residue {
            b'A' => 0,
            b'C' => 1,
            b'D' => 2,
            b'E' => 3,
            b'F' => 4,
            b'G' => 5,
            b'H' => 6,
            b'I' => 7,
            b'K' => 8,
            b'L' => 9,
            b'M' => 10,
            b'N' => 11,
            b'P' => 12,
            b'Q' => 13,
            b'R' => 14,
            b'S' => 15,
            b'T' => 16,
            b'V' => 17,
            b'W' => 18,
            b'Y' => 19,
            _ => return None,
        };
        if !distinct[code] {
            distinct[code] = true;
            n_distinct += 1;
        }
        value = (value << 5) | code as u64;
    }
    // A homopolymer seed matches every low-complexity protein and hides the
    // real hits. One repeated residue is enough to drop the seed.
    if n_distinct <= 1 {
        None
    } else {
        Some(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn protein(length: usize, phase: usize) -> Vec<u8> {
        let alphabet = b"ACDEFGHIKLMNPQRSTVWY";
        let mut seq: Vec<u8> = (0..length)
            .map(|i| alphabet[(i + phase) % alphabet.len()])
            .collect();
        seq[0] = b'M';
        seq
    }

    #[test]
    fn identical_sequences_share_a_family() {
        let seq = protein(40, 1);
        let ids = cluster_families(&[seq.clone(), seq], 0.98, 0.7, 0.95);
        assert_eq!(ids[0], ids[1]);
    }

    /// `score_candidates` switches to a parallel path above
    /// `PARALLEL_CANDIDATE_THRESHOLD`. Each candidate's alignment is
    /// independent of the others, and rayon's indexed `collect` is
    /// documented to preserve source order, so this must match a plain
    /// sequential scoring of the same candidates exactly - this test proves
    /// it with enough candidates (more than the threshold) to force the
    /// parallel branch.
    #[test]
    fn parallel_candidate_scoring_matches_serial() {
        let query = protein(120, 0);
        let count = PARALLEL_CANDIDATE_THRESHOLD * 2 + 3;
        let proteins: Vec<Vec<u8>> = (0..count).map(|salt| protein(120, salt)).collect();
        let indices: Vec<usize> = (0..count).collect();
        let reps: Vec<usize> = (0..count).collect();
        let candidates: Vec<u32> = (0..count as u32).collect();

        let parallel = score_candidates(&query, &candidates, &reps, &indices, &proteins, true);
        let serial: Vec<Option<(u32, f64)>> = candidates
            .iter()
            .map(|&rep_slot| {
                let rep_index = indices[reps[rep_slot as usize]];
                global_identity(&query, &proteins[rep_index]).map(|id| (rep_slot, id))
            })
            .collect();
        assert_eq!(parallel, serial);
    }

    #[test]
    fn chunked_scoring_matches_one_protein_at_a_time() {
        let mut proteins = Vec::new();
        for phase in 0..40 {
            let full = protein(80, phase);
            proteins.push(full.clone());
            let mut near = full.clone();
            near[10] = if near[10] == b'A' { b'C' } else { b'A' };
            proteins.push(near);
            let mut mid = full.clone();
            for index in (6..mid.len()).step_by(5) {
                mid[index] = if mid[index] == b'A' { b'L' } else { b'A' };
            }
            proteins.push(mid);
        }
        proteins.push(block(b"ACDEFGHIKL", 90));
        proteins.push(block(b"NPQRSTVWY", 70));
        let chunked = cluster_families(&proteins, 0.98, 0.70, 0.95);
        let serial = cluster_families_scheduled(&proteins, 0.98, 0.70, 0.95, false);
        assert_eq!(chunked, serial);
        assert_eq!(chunked, cluster_families(&proteins, 0.98, 0.70, 0.95));
    }

    #[test]
    fn diverged_sequences_stay_apart() {
        // Two compositions with no shared residue. A phase shift of one
        // 20-letter cycle is still the same peptide under global alignment,
        // so that is not a diverged pair.
        let a = block(b"ACDEFGHIKL", 40);
        let b = block(b"NPQRSTVWY", 40);
        let id = crate::align::global_identity(&a, &b).unwrap();
        assert!(id < 0.7, "{id}");
        let ids = cluster_families(&[a, b], 0.98, 0.7, 0.95);
        assert_ne!(ids[0], ids[1]);
    }

    fn block(residues: &[u8], length: usize) -> Vec<u8> {
        let mut seq = Vec::with_capacity(length);
        seq.push(b'M');
        for index in 1..length {
            seq.push(residues[index % residues.len()]);
        }
        seq
    }

    #[test]
    fn a_short_fragment_is_not_pulled_in_by_homology() {
        let full = protein(80, 2);
        let mut fragment = full[..25].to_vec();
        fragment[0] = b'M';
        let ids = cluster_families(&[full, fragment], 0.98, 0.7, 0.95);
        assert_ne!(ids[0], ids[1]);
    }

    #[test]
    fn seventy_percent_identity_joins_at_the_family_step() {
        let full = protein(60, 3);
        let mut near = full.clone();
        for i in (4..near.len()).step_by(4) {
            near[i] = if near[i] == b'A' { b'L' } else { b'A' };
        }
        let id = crate::align::global_identity(&full, &near).unwrap();
        assert!(id >= 0.7 && id < 0.98, "{id}");
        let ids = cluster_families(&[full, near], 0.98, 0.7, 0.95);
        assert_eq!(ids[0], ids[1]);
    }
}
