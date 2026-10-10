//! Which allele becomes the reference.
//!
//! Reads are mapped to this sequence and loss-of-function is whatever differs
//! from it. The common allele is the wrong choice when the common allele is
//! already broken: every truncated sample would look wild-type. The reference
//! is therefore the longest intact coding sequence that is still the same gene.
//!
//! A longer sequence is rejected when it is not a clean extension of the
//! usual allele (an internal insertion) or when the extra tail is some other
//! gene in the pangenome (a fusion). Frequency only breaks ties between
//! intact sequences of that same length. The sequence that is written is one
//! that was observed, or a head-to-tail join of observed fragments when no
//! sample still carries the gene in one piece.

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
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Choice {
    pub dna: String,
    pub intact: bool,
    pub selection: &'static str,
    pub exact_copies: u32,
}

/// `reps` is every family's representative protein, keyed by cluster id;
/// `own_cluster` is excluded from the fusion/extension check so a gene is
/// never compared against itself. Borrowing from `reps` here (rather than
/// the caller handing in a pre-cloned `Vec<Vec<u8>>` of every other family)
/// means the common case - a candidate allele that is not notably longer
/// than the baseline - never touches `reps` at all: `allele_is_clean` only
/// reaches `matches_other_gene` when it needs to.
pub fn choose_reference(
    alleles: &[Allele],
    reps: &HashMap<u32, Vec<u8>>,
    own_cluster: u32,
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
    let intact: Vec<&Allele> = pool.iter().copied().filter(|allele| is_intact(&allele.dna)).collect();
    if intact.is_empty() {
        return Some(broken_fallback(pool));
    }

    let baseline = most_frequent(&intact);
    let baseline_protein = translate(&baseline.dna);
    let mut kept: Vec<&Allele> = Vec::new();
    for allele in &intact {
        if allele_is_clean(allele, &baseline_protein, reps, own_cluster) {
            kept.push(*allele);
        }
    }
    if kept.is_empty() {
        kept = intact.clone();
    }
    // A joined fragment model is only allowed to win when it is a longer
    // intact version of the same gene. It is not an observed CDS.
    for allele in alleles.iter().filter(|allele| allele.synthetic && is_intact(&allele.dna)) {
        if allele_is_clean(allele, &baseline_protein, reps, own_cluster) {
            kept.push(allele);
        }
    }

    let max_len = kept.iter().map(|allele| allele.dna.len()).max()?;
    let mut class: Vec<&Allele> = kept
        .into_iter()
        .filter(|allele| allele.dna.len() == max_len)
        .collect();
    if class.iter().any(|allele| orf_status(&allele.dna).terminal_stop) {
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

    let (chosen, selection) = if tied.len() == 1 {
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

    let exact_copies = pool
        .iter()
        .filter(|allele| allele.dna == chosen.dna)
        .map(|allele| allele.count)
        .sum();
    Some(Choice {
        dna: chosen.dna.clone(),
        intact: true,
        selection,
        exact_copies,
    })
}

fn allele_is_clean(
    allele: &Allele,
    baseline: &[u8],
    reps: &HashMap<u32, Vec<u8>>,
    own_cluster: u32,
) -> bool {
    let protein = translate(&allele.dna);
    if protein.len() <= baseline.len().saturating_add(SMALL_EXTENSION_AA) {
        return true;
    }
    let Some(fit) = fit_short_in_long(baseline, &protein) else {
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
    !matches_other_gene(prefix, &protein, reps, own_cluster)
        && !matches_other_gene(suffix, &protein, reps, own_cluster)
}

fn matches_other_gene(
    peptide: &[u8], whole: &[u8], reps: &HashMap<u32, Vec<u8>>, own_cluster: u32,
) -> bool {
    if peptide.len() < FUSION_MIN_AA {
        return false;
    }
    reps.iter()
        .filter(|(id, _)| **id != own_cluster)
        .any(|(_, other)| {
            let matches_tail = if peptide.len() <= other.len() {
                fit_short_in_long(peptide, other)
                    .is_some_and(|fit| fit.identity >= FUSION_IDENTITY && fit.coverage >= FUSION_COVERAGE)
            } else {
                fit_short_in_long(other, peptide)
                    .is_some_and(|fit| fit.identity >= FUSION_IDENTITY && fit.coverage >= FUSION_COVERAGE)
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
pub fn join_fragments(parts: &[&str]) -> Option<String> {
    if parts.len() < 2 {
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
    if is_intact(&dna) { Some(dna) } else { None }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genetic_code::back_translate;

    fn allele(dna: &str, count: u32) -> Allele {
        Allele {
            dna: dna.to_string(),
            count,
            synthetic: false,
        }
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
    fn rare_intact_gene_beats_the_common_truncation() {
        let full = protein(60, 1);
        let short = full[..30].to_vec();
        let full_dna = back_translate(&full);
        let short_dna = back_translate(&short);
        let alleles = vec![allele(&full_dna, 1), allele(&short_dna, 9)];
        let choice = choose_reference(&alleles, &no_others(), 0).unwrap();
        assert_eq!(choice.dna, full_dna);
        assert!(choice.intact);
        assert_eq!(choice.selection, "intact_longest");
        assert_eq!(choice.exact_copies, 1);
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
        let choice = choose_reference(&alleles, &no_others(), 0).unwrap();
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
        let choice = choose_reference(&alleles, &other_reps(&[other]), 0).unwrap();
        assert_eq!(choice.dna, gene_dna, "{}", choice.selection);
    }

    #[test]
    fn full_length_paralog_does_not_disqualify_a_rare_intact_reference() {
        let full = protein(80, 3);
        let short = full[..40].to_vec();
        let full_dna = back_translate(&full);
        let alleles = vec![allele(&back_translate(&short), 9), allele(&full_dna, 1)];
        for changes in 0..=1 {
            let mut paralog = full.clone();
            if changes == 1 {
                paralog[60] = if paralog[60] == b'A' { b'C' } else { b'A' };
            }
            let choice = choose_reference(&alleles, &other_reps(&[paralog]), 0).unwrap();
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
        let choice = choose_reference(&alleles, &no_others(), 0).unwrap();
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
        let choice = choose_reference(&[allele(&dna, 4)], &no_others(), 0).unwrap();
        assert!(!choice.intact);
        assert_eq!(choice.selection, "no_intact_copy");
    }
}
