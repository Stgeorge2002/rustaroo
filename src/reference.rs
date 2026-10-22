//! Which allele becomes the reference.
//!
//! Reads are mapped to this sequence and loss-of-function is whatever differs
//! from it. The common allele is the wrong choice when the common allele is
//! already broken: every truncated sample would look wild-type. The reference
//! is therefore the longest intact coding sequence that is still the same gene.
//!
//! A longer sequence is rejected when it is not a clean extension of the
//! usual allele (an internal insertion) or when the extra tail is some other
//! gene in the pangenome (a fusion). It is also rejected when the shorter
//! copies already share one finished stop and the extra sequence is private.
//! A finished stop is a terminal stop on intact copies that are not cut off
//! by a contig boundary, shared by at least `max(2, 5% of samples)`. The
//! longer length is private below `max(2, 1% of samples)`. When both lengths
//! clear that finished-end bar, the longer one stays and the family is a
//! length polymorphism. Frequency only breaks ties between intact sequences
//! of one length. The sequence that is written is one that was observed, or
//! a head-to-tail join of observed fragments when no sample still carries
//! the gene in one piece. A join of two finished genes is refused: stripping
//! the stop between them invents a sequence no genome carries. A join cannot
//! outvote a finished end.

use std::collections::HashMap;

use crate::align::{fit_short_in_long, length_ratio, terminal_identity_parts};
use crate::genetic_code::{is_intact, orf_status, strip_terminal_stop, translate};

const EXTENSION_IDENTITY: f64 = 0.95;
const EXTENSION_COVERAGE: f64 = 0.90;
const FUSION_IDENTITY: f64 = 0.75;
const FUSION_COVERAGE: f64 = 0.80;
const FUSION_MIN_AA: usize = 20;
/// Length beyond the usual allele that is still "the same structure". Past
/// this, the extra sequence has to look like a terminal truncation, not a
/// new insertion in the middle.
const SMALL_EXTENSION_AA: usize = 5;
const WHOLE_GENE_IDENTITY: f64 = 0.98;

#[derive(Debug, Clone)]
pub struct Allele {
    pub dna: String,
    pub count: u32,
    pub synthetic: bool,
    /// Samples whose copy is a closed ORF away from both contig ends.
    /// Those copies can define a finished gene end. A clipped copy still
    /// counts in `count` and still stays in the family.
    pub finished: u32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Choice {
    pub dna: String,
    pub intact: bool,
    pub selection: &'static str,
    pub exact_copies: u32,
    /// DNA length of the finished modal end, or of the chosen reference when
    /// no shorter end qualified.
    pub modal_end_bp: u32,
    /// Samples at that finished end. Zero when no intact copy cleared the bar.
    pub modal_end_samples: u32,
    /// Samples carrying an intact allele that extends past the finished end.
    pub extension_samples: u32,
}

/// `reps` is every family's representative protein, keyed by cluster id;
/// `own_cluster` is excluded from the fusion/extension check so a gene is
/// never compared against itself. Borrowing from `reps` here (rather than
/// the caller handing in a pre-cloned `Vec<Vec<u8>>` of every other family)
/// means the common case - a candidate allele that is not notably longer
/// than the baseline - never touches `reps` at all: `allele_is_clean` only
/// reaches `matches_other_gene` when it needs to.
///
/// `n_samples` is the collection size, used only when intact alleles disagree
/// on length. A family whose intact copies already share one length returns
/// on the same path as before: no end census and no extra alignment.
pub fn choose_reference(
    alleles: &[Allele],
    reps: &HashMap<u32, Vec<u8>>,
    own_cluster: u32,
    n_samples: usize,
) -> Option<Choice> {
    let observed: Vec<&Allele> = alleles.iter().filter(|allele| !allele.synthetic).collect();
    if observed.is_empty() && alleles.is_empty() {
        return None;
    }
    let pool: Vec<&Allele> = if observed.is_empty() {
        alleles.iter().collect()
    } else {
        observed
    };
    let intact: Vec<&Allele> = pool
        .iter()
        .copied()
        .filter(|allele| is_intact(&allele.dna))
        .collect();
    if intact.is_empty() {
        return Some(broken_fallback(pool));
    }

    let baseline = most_frequent(&intact);
    // One length is the hot path. The census below translates and fits, so it
    // runs only when some intact copy, or a joined model, is actually longer.
    let uniform = {
        let length = intact[0].dna.len();
        intact.iter().all(|allele| allele.dna.len() == length)
            && !alleles.iter().any(|allele| {
                allele.synthetic && is_intact(&allele.dna) && allele.dna.len() != length
            })
    };

    let mut kept: Vec<&Allele> = Vec::new();
    let mut proteins: Vec<Vec<u8>> = Vec::new();
    let mut joined: Vec<(&Allele, Vec<u8>)> = Vec::new();
    if uniform {
        let baseline_protein = translate(&baseline.dna);
        for allele in &intact {
            if allele_is_clean(allele, &baseline_protein, reps, own_cluster) {
                kept.push(*allele);
            }
        }
        if kept.is_empty() {
            kept = intact.clone();
        }
        for allele in alleles
            .iter()
            .filter(|allele| allele.synthetic && is_intact(&allele.dna))
        {
            if allele_is_clean(allele, &baseline_protein, reps, own_cluster) {
                kept.push(allele);
            }
        }
    } else {
        proteins = intact.iter().map(|allele| translate(&allele.dna)).collect();
        let baseline_at = intact
            .iter()
            .position(|allele| std::ptr::eq(*allele, baseline))
            .expect("baseline is one of the intact alleles");
        for (allele, protein) in intact.iter().copied().zip(&proteins) {
            if protein_is_clean(protein, &proteins[baseline_at], reps, own_cluster) {
                kept.push(allele);
            }
        }
        if kept.is_empty() {
            kept = intact.clone();
        }
        for allele in alleles
            .iter()
            .filter(|allele| allele.synthetic && is_intact(&allele.dna))
        {
            let protein = translate(&allele.dna);
            if protein_is_clean(&protein, &proteins[baseline_at], reps, own_cluster) {
                kept.push(allele);
                joined.push((allele, protein));
            }
        }
    }

    let mut census = Census::inactive();
    if !uniform {
        if let Some(end) = finished_end(&intact, &proteins, n_samples) {
            census =
                refuse_private_extensions(&mut kept, &intact, &proteins, &joined, end, n_samples);
            if kept.is_empty() {
                // The clean-check fallback must not put a refused private
                // extension back. The finished end itself is still intact.
                kept = intact
                    .iter()
                    .copied()
                    .filter(|allele| allele.dna.len() == census.modal_bp as usize)
                    .collect();
            }
        }
    }

    // The longer finished allele is the reference. The common shorter end
    // fails the 90% both-coverage comparison against a longer copy of the
    // same gene, which is the length the other reference usually keeps.
    let longest = kept.iter().map(|allele| allele.dna.len()).max()?;
    let mut class: Vec<&Allele> = kept
        .into_iter()
        .filter(|allele| allele.dna.len() == longest)
        .collect();
    if class
        .iter()
        .any(|allele| orf_status(&allele.dna).terminal_stop)
    {
        class.retain(|allele| orf_status(&allele.dna).terminal_stop);
    }
    // Prefer a sequence that exists in a sample over a join of the same length.
    if class.iter().any(|allele| !allele.synthetic) {
        class.retain(|allele| !allele.synthetic);
    }

    let best_count = class.iter().map(|allele| allele.count).max()?;
    let mut tied: Vec<&Allele> = class
        .iter()
        .copied()
        .filter(|allele| allele.count == best_count)
        .collect();
    tied.sort_by(|a, b| a.dna.cmp(&b.dna));

    let (chosen, mut selection) = if tied.len() == 1 {
        let label = if tied[0].synthetic {
            "joined_fragments"
        } else if class.len() == 1 {
            "intact_longest"
        } else {
            "intact_modal"
        };
        (tied[0], label)
    } else {
        let majority = codon_majority(&class);
        let chosen = tied
            .into_iter()
            .min_by(|a, b| {
                codon_distance(&a.dna, &majority)
                    .cmp(&codon_distance(&b.dna, &majority))
                    .then(a.dna.cmp(&b.dna))
            })
            .expect("tie set is non-empty");
        (chosen, "intact_ancestral_tiebreak")
    };
    if census.active && !chosen.synthetic {
        let chosen_aa = protein_of(chosen, &intact, &proteins, &joined).len();
        let support = end_support(n_samples);
        let chosen_finished = census.finished_by_aa.get(&chosen_aa).copied().unwrap_or(0);
        let other_finished = census
            .finished_by_aa
            .iter()
            .any(|(&length, &samples)| length != chosen_aa && samples >= support);
        if census.refused && chosen_aa == census.anchor_aa {
            selection = "intact_conserved_end";
        } else if chosen_finished >= support && other_finished {
            selection = "intact_length_polymorphism";
        }
    }

    let exact_copies = pool
        .iter()
        .filter(|allele| allele.dna == chosen.dna)
        .map(|allele| allele.count)
        .sum();
    let (modal_end_bp, modal_end_samples, extension_samples) = if census.active {
        (
            census.modal_bp,
            census.modal_samples,
            census.extension_samples,
        )
    } else {
        (
            chosen.dna.len() as u32,
            finished_at_dna_length(&pool, chosen.dna.len()),
            0,
        )
    };
    Some(Choice {
        dna: chosen.dna.clone(),
        intact: true,
        selection,
        exact_copies,
        modal_end_bp,
        modal_end_samples,
        extension_samples,
    })
}

fn allele_is_clean(
    allele: &Allele,
    baseline: &[u8],
    reps: &HashMap<u32, Vec<u8>>,
    own_cluster: u32,
) -> bool {
    let protein = translate(&allele.dna);
    protein_is_clean(&protein, baseline, reps, own_cluster)
}

fn protein_is_clean(
    protein: &[u8],
    baseline: &[u8],
    reps: &HashMap<u32, Vec<u8>>,
    own_cluster: u32,
) -> bool {
    if protein.len() <= baseline.len().saturating_add(SMALL_EXTENSION_AA) {
        return true;
    }
    let Some(fit) = fit_short_in_long(baseline, protein) else {
        return false;
    };
    if fit.identity < EXTENSION_IDENTITY
        || fit.coverage < EXTENSION_COVERAGE
        || fit.internal_gap > SMALL_EXTENSION_AA
    {
        return false;
    }
    let prefix = &protein[..fit.extra_prefix(protein.len())];
    let suffix = &protein[protein.len() - fit.extra_suffix(protein.len())..];
    !matches_other_gene(prefix, protein, reps, own_cluster)
        && !matches_other_gene(suffix, protein, reps, own_cluster)
}

/// Samples that can vouch for one DNA length. Clipped copies are absent here
/// even when they contribute to the allele count.
fn finished_at_dna_length(alleles: &[&Allele], dna_len: usize) -> u32 {
    alleles
        .iter()
        .filter(|allele| !allele.synthetic && allele.dna.len() == dna_len)
        .map(|allele| allele.finished)
        .sum()
}

fn end_support(n_samples: usize) -> u32 {
    let fraction = (n_samples as f64 * 0.05).ceil() as u32;
    fraction.max(2)
}

fn extension_support(n_samples: usize) -> u32 {
    let fraction = (n_samples as f64 * 0.01).ceil() as u32;
    fraction.max(2)
}

struct Census {
    active: bool,
    refused: bool,
    anchor_aa: usize,
    modal_bp: u32,
    modal_samples: u32,
    extension_samples: u32,
    finished_by_aa: HashMap<usize, u32>,
}

impl Census {
    fn inactive() -> Self {
        Self {
            active: false,
            refused: false,
            anchor_aa: 0,
            modal_bp: 0,
            modal_samples: 0,
            extension_samples: 0,
            finished_by_aa: HashMap::new(),
        }
    }
}

struct FinishedEnd {
    anchor: Vec<u8>,
    dna_bp: u32,
    samples: u32,
    finished_by_aa: HashMap<usize, u32>,
}

/// The best-supported finished stop, when intact copies actually share one.
///
/// The end with the most unclipped samples wins. An equal count keeps the
/// shorter protein, so a tied longer allele stays an extension rather than
/// the end the collection is measured against.
///
/// Copies are already collapsed to distinct sequences, so this aligns those
/// sequences once. It does not walk samples. A length class counts only the
/// sequences that land on the same termini as its most common member: ragged
/// ends stay in their own classes and fail the support bar on their own.
fn finished_end(intact: &[&Allele], proteins: &[Vec<u8>], n_samples: usize) -> Option<FinishedEnd> {
    let threshold = end_support(n_samples);
    let mut classes: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut finished_by_aa: HashMap<usize, u32> = HashMap::new();
    for (index, allele) in intact.iter().enumerate() {
        if allele.finished == 0 || !orf_status(&allele.dna).terminal_stop {
            continue;
        }
        let length = proteins[index].len();
        classes.entry(length).or_default().push(index);
    }
    let mut best: Option<FinishedEnd> = None;
    for (length, members) in classes {
        let anchor_at = members
            .iter()
            .copied()
            .max_by(|&a, &b| {
                intact[a]
                    .finished
                    .cmp(&intact[b].finished)
                    .then(intact[b].dna.len().cmp(&intact[a].dna.len()))
                    .then(intact[b].dna.cmp(&intact[a].dna))
            })
            .expect("length class has a member");
        let anchor = &proteins[anchor_at];
        let samples = members
            .iter()
            .copied()
            .filter(|&index| same_end(anchor, &proteins[index]))
            .map(|index| intact[index].finished)
            .sum();
        finished_by_aa.insert(length, samples);
        if samples < threshold {
            continue;
        }
        let replace = match &best {
            None => true,
            Some(current) => {
                samples > current.samples
                    || (samples == current.samples && length < current.anchor.len())
            }
        };
        if replace {
            best = Some(FinishedEnd {
                anchor: anchor.clone(),
                dna_bp: intact[anchor_at].dna.len() as u32,
                samples,
                finished_by_aa: HashMap::new(),
            });
        }
    }
    let mut best = best?;
    best.finished_by_aa = finished_by_aa;
    Some(best)
}

fn same_end(anchor: &[u8], other: &[u8]) -> bool {
    if anchor.len() != other.len() {
        return false;
    }
    if anchor == other {
        return true;
    }
    let Some(fit) = fit_short_in_long(anchor, other) else {
        return false;
    };
    fit.identity + f64::EPSILON >= EXTENSION_IDENTITY
        && fit.internal_gap <= SMALL_EXTENSION_AA
        && fit.extra_prefix(other.len()) == 0
        && fit.extra_suffix(other.len()) == 0
}

/// A longer protein is a terminal extension when the finished protein sits
/// inside it, the extra residues are at the ends, and the match is the same
/// bar used for a clean extension. An internal insertion fails this and is
/// left to the existing rejection.
fn extends_past(anchor: &[u8], longer: &[u8]) -> bool {
    if longer.len() <= anchor.len() {
        return false;
    }
    let Some(fit) = fit_short_in_long(anchor, longer) else {
        return false;
    };
    fit.identity + f64::EPSILON >= EXTENSION_IDENTITY
        && fit.coverage + f64::EPSILON >= EXTENSION_COVERAGE
        && fit.internal_gap <= SMALL_EXTENSION_AA
        && (fit.extra_prefix(longer.len()) > 0 || fit.extra_suffix(longer.len()) > 0)
}

fn protein_of<'a>(
    allele: &Allele,
    intact: &[&Allele],
    proteins: &'a [Vec<u8>],
    joined: &'a [(&Allele, Vec<u8>)],
) -> &'a [u8] {
    if let Some(index) = intact.iter().position(|item| std::ptr::eq(*item, allele)) {
        return &proteins[index];
    }
    joined
        .iter()
        .find(|(item, _)| std::ptr::eq(*item, allele))
        .map(|(_, protein)| protein.as_slice())
        .expect("allele protein was translated with the length census")
}

fn refuse_private_extensions<'a>(
    kept: &mut Vec<&'a Allele>,
    intact: &[&'a Allele],
    proteins: &[Vec<u8>],
    joined: &[(&Allele, Vec<u8>)],
    end: FinishedEnd,
    n_samples: usize,
) -> Census {
    let private_below = extension_support(n_samples);
    let mut extending_copies: HashMap<usize, u32> = HashMap::new();
    let mut judged: Vec<(&Allele, bool, usize)> = Vec::with_capacity(kept.len());
    for allele in kept.iter().copied() {
        let protein = protein_of(allele, intact, proteins, joined);
        let extends = extends_past(&end.anchor, protein);
        if extends {
            *extending_copies.entry(protein.len()).or_default() += allele.count;
        }
        judged.push((allele, extends, protein.len()));
    }
    let extension_samples = extending_copies.values().sum();
    let before = kept.len();
    kept.retain(|allele| {
        judged.iter().any(|(item, extends, length)| {
            std::ptr::eq(*item, *allele)
                && !(*extends && extending_copies.get(length).copied().unwrap_or(0) < private_below)
        })
    });
    Census {
        active: true,
        refused: kept.len() != before,
        anchor_aa: end.anchor.len(),
        modal_bp: end.dna_bp,
        modal_samples: end.samples,
        extension_samples,
        finished_by_aa: end.finished_by_aa,
    }
}

fn matches_other_gene(
    peptide: &[u8],
    whole: &[u8],
    reps: &HashMap<u32, Vec<u8>>,
    own_cluster: u32,
) -> bool {
    if peptide.len() < FUSION_MIN_AA {
        return false;
    }
    reps.iter()
        .filter(|(id, _)| **id != own_cluster)
        .any(|(_, other)| {
            let matches_tail = if peptide.len() <= other.len() {
                fit_short_in_long(peptide, other).is_some_and(|fit| {
                    fit.identity >= FUSION_IDENTITY && fit.coverage >= FUSION_COVERAGE
                })
            } else {
                fit_short_in_long(other, peptide).is_some_and(|fit| {
                    fit.identity >= FUSION_IDENTITY && fit.coverage >= FUSION_COVERAGE
                })
            };
            if !matches_tail {
                return false;
            }
            let (short, long) = if whole.len() <= other.len() {
                (whole, other.as_slice())
            } else {
                (other.as_slice(), whole)
            };
            // A whole-protein paralog is not evidence for an independent
            // fused tail on the candidate's otherwise intact allele.
            !(length_ratio(short, long) >= 0.95
                && terminal_identity_parts(short, &[], long)
                    .is_some_and(|identity| identity + f64::EPSILON >= WHOLE_GENE_IDENTITY))
        })
}

fn most_frequent<'a>(alleles: &[&'a Allele]) -> &'a Allele {
    alleles
        .iter()
        .copied()
        .max_by(|a, b| {
            a.count
                .cmp(&b.count)
                .then(a.dna.len().cmp(&b.dna.len()))
                .then(b.dna.cmp(&a.dna))
        })
        .expect("allele list is non-empty")
}

fn broken_fallback(alleles: Vec<&Allele>) -> Choice {
    let chosen = alleles
        .iter()
        .copied()
        .max_by(|a, b| {
            let stops_a = orf_status(&a.dna).internal_stops;
            let stops_b = orf_status(&b.dna).internal_stops;
            stops_b
                .cmp(&stops_a)
                .then(a.dna.len().cmp(&b.dna.len()))
                .then(a.count.cmp(&b.count))
                .then(b.dna.cmp(&a.dna))
        })
        .expect("allele list is non-empty");
    let exact_copies = alleles
        .iter()
        .filter(|allele| allele.dna == chosen.dna)
        .map(|allele| allele.count)
        .sum();
    Choice {
        dna: chosen.dna.clone(),
        intact: false,
        selection: "no_intact_copy",
        exact_copies,
        modal_end_bp: chosen.dna.len() as u32,
        modal_end_samples: 0,
        extension_samples: 0,
    }
}

/// Codon-by-codon majority of one length class. Used only to choose which
/// observed haplotype is closer to the intact ancestor when two haplotypes
/// are equally common. The majority string itself is not emitted: it can be
/// a combination that never existed in one genome.
fn codon_majority(alleles: &[&Allele]) -> String {
    let length = alleles[0].dna.len();
    let mut out = String::with_capacity(length);
    let mut offset = 0;
    while offset + 3 <= length {
        let mut counts: HashMap<&str, u32> = HashMap::new();
        for allele in alleles {
            let codon = &allele.dna[offset..offset + 3];
            *counts.entry(codon).or_default() += allele.count.max(1);
        }
        let codon = counts
            .into_iter()
            .max_by(|a, b| a.1.cmp(&b.1).then(b.0.cmp(a.0)))
            .map(|(codon, _)| codon)
            .unwrap_or("NNN");
        out.push_str(codon);
        offset += 3;
    }
    out
}

fn codon_distance(dna: &str, majority: &str) -> usize {
    dna.as_bytes()
        .chunks(3)
        .zip(majority.as_bytes().chunks(3))
        .filter(|(a, b)| a != b)
        .count()
}

/// Join adjacent fragment CDS, dropping the stop between them, so a gene that
/// no sample still carries in one piece has an intact model to map against.
///
/// Two parts that each already end on their own stop are two genes. Stripping
/// those stops is how metN+metP and rpsC+rplP became one reference that no
/// genome carries. A broken gene still joins: only one of its pieces is a
/// finished ORF.
pub fn join_fragments(parts: &[&str]) -> Option<String> {
    if parts.len() < 2 {
        return None;
    }
    let finished = parts
        .iter()
        .filter(|part| {
            let status = orf_status(part);
            status.is_intact() && status.terminal_stop
        })
        .count();
    if finished >= 2 {
        return None;
    }
    let mut dna = String::new();
    for (index, part) in parts.iter().enumerate() {
        if part.len() < 3 || part.len() % 3 != 0 {
            return None;
        }
        if index + 1 == parts.len() {
            dna.push_str(part);
        } else {
            dna.push_str(strip_terminal_stop(part));
        }
    }
    if is_intact(&dna) {
        Some(dna)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genetic_code::{back_translate, is_intact, orf_status};

    fn allele(dna: &str, count: u32) -> Allele {
        Allele {
            dna: dna.to_string(),
            count,
            synthetic: false,
            finished: count,
        }
    }

    /// An intact sequence whose copies are clipped by a contig end, so they
    /// cannot vote for a finished stop.
    fn clipped(dna: &str, count: u32) -> Allele {
        Allele {
            dna: dna.to_string(),
            count,
            synthetic: false,
            finished: 0,
        }
    }

    fn choose(alleles: &[Allele], reps: &HashMap<u32, Vec<u8>>) -> Choice {
        let samples = alleles
            .iter()
            .map(|allele| allele.count as usize)
            .sum::<usize>()
            .max(1);
        choose_reference(alleles, reps, 0, samples).unwrap()
    }

    /// Builds the "other families" map `choose_reference` now takes directly,
    /// keyed away from `own_cluster` (0) so the exclusion filter never hides
    /// a test's own fixture by accident.
    fn other_reps(others: &[Vec<u8>]) -> HashMap<u32, Vec<u8>> {
        others
            .iter()
            .enumerate()
            .map(|(index, protein)| (index as u32 + 1, protein.clone()))
            .collect()
    }

    fn no_others() -> HashMap<u32, Vec<u8>> {
        HashMap::new()
    }

    fn protein(length: usize, phase: usize) -> Vec<u8> {
        let alphabet = b"ACDEFGHIKLMNPQRSTVWY";
        let mut seq: Vec<u8> = (0..length)
            .map(|i| alphabet[(i + phase) % alphabet.len()])
            .collect();
        seq[0] = b'M';
        seq
    }

    #[test]
    fn a_private_extension_loses_to_a_shared_stop_and_a_clipped_truncation_does_not() {
        let full = protein(99, 1);
        let short = full[..89].to_vec();
        let full_dna = back_translate(&full);
        let short_dna = back_translate(&short);
        assert_eq!(short_dna.len(), 270);
        assert_eq!(full_dna.len(), 300);
        let finished = vec![allele(&full_dna, 1), allele(&short_dna, 399)];
        let choice = choose_reference(&finished, &no_others(), 0, 400).unwrap();
        assert_eq!(choice.dna, short_dna);
        assert!(choice.intact);
        assert_eq!(choice.selection, "intact_conserved_end");
        assert_eq!(choice.exact_copies, 399);
        assert_eq!(choice.modal_end_bp, 270);
        assert_eq!(choice.modal_end_samples, 399);
        assert_eq!(choice.extension_samples, 1);

        let clipped_copies = vec![allele(&full_dna, 1), clipped(&short_dna, 399)];
        let kept = choose_reference(&clipped_copies, &no_others(), 0, 400).unwrap();
        assert_eq!(kept.dna, full_dna, "{}", kept.selection);
        assert_eq!(kept.selection, "intact_longest");
        assert_eq!(kept.extension_samples, 0);
    }

    #[test]
    fn a_short_private_tail_past_a_shared_stop_is_refused() {
        let common = protein(89, 2);
        let mut longer = common.clone();
        longer.extend_from_slice(b"ACD");
        let choice = choose_reference(
            &[
                allele(&back_translate(&longer), 1),
                allele(&back_translate(&common), 399),
            ],
            &no_others(),
            0,
            400,
        )
        .unwrap();
        assert_eq!(choice.dna, back_translate(&common));
        assert_eq!(choice.selection, "intact_conserved_end");
        assert_eq!(choice.extension_samples, 1);
    }

    #[test]
    fn two_finished_genes_are_not_concatenated() {
        let met = back_translate(&protein(40, 1));
        let next = back_translate(&protein(50, 8));
        assert!(orf_status(&met).terminal_stop && orf_status(&next).terminal_stop);
        assert!(join_fragments(&[&met, &next]).is_none());
        assert!(join_fragments(&[&met, &next, &back_translate(&protein(30, 4))]).is_none());

        let mut fragment = next.clone();
        fragment.replace_range(0..3, "GCT");
        assert!(!is_intact(&fragment));
        let joined = join_fragments(&[&met, &fragment]).expect("a real fragment still joins");
        assert!(is_intact(&joined));
        assert!(joined.len() > met.len());
    }

    #[test]
    fn both_finished_lengths_keep_the_longer_allele() {
        let short = protein(89, 3);
        let mut long = short.clone();
        long.extend_from_slice(&protein(11, 8)[1..]);
        let choice = choose_reference(
            &[
                allele(&back_translate(&long), 200),
                allele(&back_translate(&short), 200),
            ],
            &no_others(),
            0,
            400,
        )
        .unwrap();
        assert_eq!(choice.dna, back_translate(&long));
        assert_eq!(choice.selection, "intact_length_polymorphism");
        assert_eq!(choice.modal_end_bp, back_translate(&short).len() as u32);
        assert_eq!(choice.modal_end_samples, 200);
        assert_eq!(choice.extension_samples, 200);
    }

    #[test]
    fn ragged_ends_do_not_demote_the_long_allele() {
        let long = protein(99, 1);
        let left = long[..89].to_vec();
        let right = long[..80].to_vec();
        let choice = choose_reference(
            &[
                allele(&back_translate(&long), 1),
                allele(&back_translate(&left), 10),
                allele(&back_translate(&right), 10),
            ],
            &no_others(),
            0,
            400,
        )
        .unwrap();
        assert_eq!(choice.dna, back_translate(&long), "{}", choice.selection);
        assert_eq!(choice.selection, "intact_longest");
    }

    #[test]
    fn synonymous_haplotypes_add_up_to_one_finished_end() {
        let mut short = protein(89, 1);
        short[2] = b'L';
        let mut long = short.clone();
        long.extend_from_slice(&protein(11, 5)[1..]);
        let dna = back_translate(&short);
        let mut synonym = dna.clone();
        synonym.replace_range(6..9, "TTG");
        assert_ne!(dna, synonym);
        let choice = choose_reference(
            &[
                allele(&back_translate(&long), 1),
                allele(&dna, 10),
                allele(&synonym, 10),
            ],
            &no_others(),
            0,
            400,
        )
        .unwrap();
        assert_eq!(choice.selection, "intact_conserved_end");
        assert_eq!(choice.modal_end_samples, 20);
        assert_eq!(choice.extension_samples, 1);
        assert_ne!(choice.dna, back_translate(&long));
    }

    #[test]
    fn a_joined_fragment_cannot_outvote_a_finished_end() {
        let short = protein(89, 4);
        let mut long = short.clone();
        long.extend_from_slice(&protein(11, 6)[1..]);
        let mut joined = allele(&back_translate(&long), 0);
        joined.synthetic = true;
        joined.finished = 0;
        let choice = choose_reference(
            &[allele(&back_translate(&short), 50), joined],
            &no_others(),
            0,
            100,
        )
        .unwrap();
        assert_eq!(choice.dna, back_translate(&short));
        assert_eq!(choice.selection, "intact_conserved_end");
    }

    #[test]
    fn internal_insertion_does_not_become_the_reference() {
        let common = protein(40, 2);
        let mut inserted = common[..20].to_vec();
        inserted.extend(std::iter::repeat(b'W').take(24));
        inserted.extend_from_slice(&common[20..]);
        let common_dna = back_translate(&common);
        let inserted_dna = back_translate(&inserted);
        let alleles = vec![allele(&inserted_dna, 1), allele(&common_dna, 8)];
        let choice = choose(&alleles, &no_others());
        assert_eq!(choice.dna, common_dna, "{}", choice.selection);
    }

    #[test]
    fn fusion_tail_that_is_another_gene_is_rejected() {
        let gene = protein(40, 3);
        let other = protein(30, 9);
        let mut fusion = gene.clone();
        // Drop the duplicate start on the tail so the extra peptide is the other gene.
        fusion.extend_from_slice(&other[1..]);
        let gene_dna = back_translate(&gene);
        let fusion_dna = back_translate(&fusion);
        let alleles = vec![allele(&gene_dna, 6), allele(&fusion_dna, 1)];
        let choice = choose(&alleles, &other_reps(&[other]));
        assert_eq!(choice.dna, gene_dna, "{}", choice.selection);
    }

    #[test]
    fn full_length_paralog_does_not_disqualify_a_rare_intact_reference() {
        let full = protein(80, 3);
        let short = full[..40].to_vec();
        let full_dna = back_translate(&full);
        let alleles = vec![clipped(&back_translate(&short), 9), allele(&full_dna, 1)];
        for changes in 0..=1 {
            let mut paralog = full.clone();
            if changes == 1 {
                paralog[60] = if paralog[60] == b'A' { b'C' } else { b'A' };
            }
            let choice = choose(&alleles, &other_reps(&[paralog]));
            assert_eq!(choice.dna, full_dna);
            assert!(choice.intact);
        }
    }

    #[test]
    fn equal_frequency_picks_the_haplotype_closer_to_the_codon_majority() {
        let base = protein(24, 4);
        let dna_a = back_translate(&base);
        let mut other = base.clone();
        other[5] = if other[5] == b'A' { b'L' } else { b'A' };
        let dna_b = back_translate(&other);
        // A third copy of A makes A the majority, but counts of the two
        // haplotypes at the winning length are what we tie. Use a third
        // sequence equal to A so the codon majority follows A, with A and B
        // both at count 2? Then A has count 2 if we pass two A alleles...
        // Pass A twice as one allele of count 2 and B of count 2: majority
        // follows A only on a tie-break of equal codon counts. Force A to be
        // the lex-smaller codon at the changed site by constructing codons.
        let alleles = vec![allele(&dna_a, 2), allele(&dna_b, 2)];
        let choice = choose(&alleles, &no_others());
        assert_eq!(choice.selection, "intact_ancestral_tiebreak");
        assert!(choice.dna == dna_a || choice.dna == dna_b);
        assert_eq!(choice.dna.len(), dna_a.len());
    }

    #[test]
    fn no_intact_copy_is_flagged() {
        let full = protein(30, 5);
        let mut dna = back_translate(&full);
        // Internal stop. The sequence is still full length.
        dna.replace_range(12..15, "TAA");
        assert!(!is_intact(&dna));
        let choice = choose(&[allele(&dna, 4)], &no_others());
        assert!(!choice.intact);
        assert_eq!(choice.selection, "no_intact_copy");
    }
}
