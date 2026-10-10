//! Two-step homology.
//!
//! Identical proteins collapse by hash. What remains is clustered longest-first,
//! the same order CD-HIT uses, so a fragment cannot become the representative
//! of a full gene. The first pass is 98% identity. Families are a second pass
//! on those representatives at the family threshold. Length ratio stays strict
//! here; a half-length fragment is the synteny step's problem, not this one's.

use std::collections::HashMap;

use rayon::prelude::*;

use crate::align::{global_identity, length_ratio};

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
    if proteins.is_empty() {
        return Vec::new();
    }
    let (unique, back) = unique_proteins(proteins);
    let (pre, pre_reps) = greedy(&unique, proteins, precluster_identity, length_ratio_min);
    let rep_indices: Vec<usize> = pre_reps.iter().map(|&slot| unique[slot]).collect();
    let (families, _) = if rep_indices.len() <= 1 {
        (vec![0u32; rep_indices.len()], Vec::new())
    } else {
        greedy(&rep_indices, proteins, family_identity, length_ratio_min)
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
) -> (Vec<u32>, Vec<usize>) {
    let mut order: Vec<usize> = (0..indices.len()).collect();
    order.sort_by(|&a, &b| {
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
    let mut postings: HashMap<u64, Vec<u32>> = HashMap::new();
    let total = order.len();
    let step = (total / 20).max(1);
    if total >= 100 {
        eprintln!(
            "rustroo: {identity:.2} identity pass, {total} proteins, word size {kmer_size}, 1 thread for decisions"
        );
    }
    for (placed, pos) in order.into_iter().enumerate() {
        let protein = &proteins[indices[pos]];
        let protein_kmers = &kmer_lists[pos];
        let candidates = candidates_for(
            protein,
            &reps,
            indices,
            proteins,
            &postings,
            length_ratio_min,
            protein_kmers,
        );
        // Scoring each candidate against the query is independent of every
        // other candidate - none of them changes while we evaluate this one
        // query - so the alignments themselves can run in parallel. Which
        // candidate *wins* still has to be decided by walking the results in
        // the same fixed order the serial version used, so a near-tie keeps
        // resolving to the earlier-ranked candidate exactly as before.
        let scored = score_candidates(protein, &candidates, &reps, indices, proteins);
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
        if let Some((rep_slot, _)) = best {
            assignment[pos] = assignment[reps[rep_slot as usize]];
        } else {
            let cluster = reps.len() as u32;
            assignment[pos] = cluster;
            let rep_slot = reps.len() as u32;
            reps.push(pos);
            for &kmer in protein_kmers {
                let list = postings.entry(kmer).or_default();
                if list.len() < SATURATED_POSTING {
                    list.push(rep_slot);
                }
            }
        }
        let seen = placed + 1;
        if total >= 100 && (seen == total || seen % step == 0) {
            eprintln!(
                "rustroo: {identity:.2} identity {seen}/{total} proteins, {} clusters",
                reps.len()
            );
        }
    }
    (assignment, reps)
}

fn candidates_for(
    protein: &[u8],
    reps: &[usize],
    indices: &[usize],
    proteins: &[Vec<u8>],
    postings: &HashMap<u64, Vec<u32>>,
    length_ratio_min: f64,
    protein_kmers: &[u64],
) -> Vec<u32> {
    if reps.is_empty() {
        return Vec::new();
    }
    let mut hits: HashMap<u32, u16> = HashMap::new();
    for kmer in protein_kmers {
        let Some(list) = postings.get(kmer) else {
            continue;
        };
        if list.len() >= SATURATED_POSTING {
            continue;
        }
        for &rep in list {
            *hits.entry(rep).or_default() += 1;
        }
    }
    let mut ranked: Vec<(u16, u32)> = hits
        .into_iter()
        .filter(|(_, count)| *count >= 2)
        .map(|(rep, count)| (count, rep))
        .filter(|(_, rep)| {
            let rep_protein = &proteins[indices[reps[*rep as usize]]];
            length_ratio(protein, rep_protein) + f64::EPSILON >= length_ratio_min
        })
        .collect();
    ranked.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
    ranked.truncate(MAX_CANDIDATES);
    ranked.into_iter().map(|(_, rep)| rep).collect()
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
) -> Vec<Option<(u32, f64)>> {
    let score_one = |&rep_slot: &u32| -> Option<(u32, f64)> {
        let rep_index = indices[reps[rep_slot as usize]];
        global_identity(protein, &proteins[rep_index]).map(|id| (rep_slot, id))
    };
    if candidates.len() >= PARALLEL_CANDIDATE_THRESHOLD {
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

        let parallel = score_candidates(&query, &candidates, &reps, &indices, &proteins);
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
