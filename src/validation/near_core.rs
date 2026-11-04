use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::Write;
use std::path::Path;

use rayon::prelude::*;

use super::{
    assess, bundle_evidence, fingerprint, list, load_in, loci, writer, Input, Locus, LogGuard,
};
use crate::log::{self, progress};
use crate::{Error, Result};

const MIN_IDENTITY: f64 = 0.98;
const MIN_COVERAGE: f64 = 0.95;
const MIN_AGREEMENT: f64 = 0.90;
const MIN_CONTEXT_SAMPLES: usize = 3;
const MIN_PROTEIN_AA: usize = 30;

#[derive(Clone, Debug)]
pub struct NearCoreOptions {
    pub threads: usize,
    pub core_fraction: f64,
    pub near_core_fraction: f64,
    pub max_alignment_cells: usize,
    pub test_merges: bool,
    pub max_merges: usize,
    pub comparison_core_count: Option<usize>,
}

impl Default for NearCoreOptions {
    fn default() -> Self {
        Self {
            threads: 0,
            core_fraction: 0.99,
            near_core_fraction: 0.95,
            max_alignment_cells: 5_000_000,
            test_merges: false,
            max_merges: 10,
            comparison_core_count: None,
        }
    }
}

#[derive(Debug)]
pub struct NearCoreSummary {
    pub candidate_pairs: usize,
    pub strong_pairs: usize,
    pub unassessed_pairs: usize,
    pub simulated_merges: usize,
    pub core_before: usize,
    pub core_after: Option<usize>,
}

#[derive(Default)]
struct Evidence {
    samples: BTreeSet<usize>,
    loci_by_sample: BTreeMap<usize, Vec<usize>>,
    contexts: BTreeMap<(usize, usize), BTreeSet<usize>>,
    complete_context_samples: BTreeSet<usize>,
    multilocus_samples: usize,
    uncertain_samples: BTreeSet<usize>,
    annotated_reference: bool,
}

fn family_evidence(input: &Input, loci: &[Locus]) -> Vec<Evidence> {
    let mut evidence: Vec<_> = (0..input.families.len()).map(|_| Evidence::default()).collect();
    for (index, locus) in loci.iter().enumerate() {
        let family = &input.families[locus.family];
        let row = &mut evidence[locus.family];
        row.samples.insert(locus.sample);
        row.loci_by_sample.entry(locus.sample).or_default().push(index);
        let uncertain = bundle_evidence(locus, input).ends_with("requires_review")
            || locus.members.iter().any(|&member| {
                let member = &input.members[member];
                u64::from(member.end) - u64::from(member.start) + 1
                    != input.alleles[member.allele].dna.len() as u64
            });
        if uncertain {
            row.uncertain_samples.insert(locus.sample);
        }
        let annotated = locus.members.iter().any(|&member| {
            let member = &input.members[member];
            let allele = &input.alleles[member.allele];
            member.origin == "annotated"
                && allele.orf.is_intact()
                && allele.protein.len() as f64 >= family.reference.protein.len() as f64 * MIN_COVERAGE
        });
        row.annotated_reference |= locus.members.iter().any(|&member| {
            let member = &input.members[member];
            member.origin == "annotated"
                && input.alleles[member.allele].dna == family.reference.dna
                && input.alleles[member.allele].orf.is_intact()
        });
        if annotated && !uncertain {
            if let (Some(upstream), Some(downstream)) = locus.flanks {
                if upstream != downstream && upstream != locus.family && downstream != locus.family {
                    row.contexts
                        .entry((upstream, downstream))
                        .or_default()
                        .insert(locus.sample);
                    row.complete_context_samples.insert(locus.sample);
                }
            }
        }
    }
    for row in &mut evidence {
        row.multilocus_samples = row.loci_by_sample.values().filter(|loci| loci.len() > 1).count();
    }
    evidence
}

#[derive(Clone, Copy, Default)]
struct Neighborhood {
    context: Option<(usize, usize)>,
    support_a: usize,
    support_b: usize,
    agreement_a: f64,
    agreement_b: f64,
}

fn neighborhood(a: &Evidence, b: &Evidence) -> Neighborhood {
    let mut best = Neighborhood::default();
    for (&context, samples_a) in &a.contexts {
        if let Some(samples_b) = b.contexts.get(&context) {
            let current = Neighborhood {
                context: Some(context),
                support_a: samples_a.len(),
                support_b: samples_b.len(),
                agreement_a: samples_a.len() as f64 / a.samples.len() as f64,
                agreement_b: samples_b.len() as f64 / b.samples.len() as f64,
            };
            if current.agreement_a.min(current.agreement_b)
                > best.agreement_a.min(best.agreement_b)
                || (current.agreement_a.min(current.agreement_b)
                    == best.agreement_a.min(best.agreement_b)
                    && current.support_a + current.support_b > best.support_a + best.support_b)
            {
                best = current;
            }
        }
    }
    best
}

#[derive(Clone, Copy)]
struct SequenceMetrics {
    identity: f64,
    coverage_a: f64,
    coverage_b: f64,
}

struct Candidate {
    a: usize,
    b: usize,
    union_samples: usize,
    cooccurring: Vec<usize>,
    distinct_locus: Vec<usize>,
    sequence_status: &'static str,
    sequence: Option<SequenceMetrics>,
    neighborhood: Neighborhood,
    rejection_reasons: Vec<&'static str>,
    simulated: bool,
}

impl Candidate {
    fn strong(&self) -> bool {
        self.rejection_reasons.is_empty()
    }

    fn sequence_score(&self) -> f64 {
        self.sequence
            .map(|row| row.identity * row.coverage_a.min(row.coverage_b))
            .unwrap_or(0.0)
    }
}

fn distinct_loci(a: &Evidence, b: &Evidence, sample: usize, loci: &[Locus]) -> bool {
    a.loci_by_sample[&sample].iter().any(|&x| {
        b.loci_by_sample[&sample].iter().any(|&y| {
            let (x, y) = (&loci[x], &loci[y]);
            x.contig != y.contig || x.end < y.start || y.end < x.start
        })
    })
}

fn candidate(
    a: usize,
    b: usize,
    input: &Input,
    loci: &[Locus],
    evidence: &[Evidence],
    core_samples: usize,
    cap: usize,
) -> Candidate {
    let (x, y) = (&evidence[a], &evidence[b]);
    let cooccurring: Vec<_> = x.samples.intersection(&y.samples).copied().collect();
    let distinct_locus = cooccurring
        .iter()
        .copied()
        .filter(|&sample| distinct_loci(x, y, sample, loci))
        .collect();
    let reference_a = &input.families[a].reference;
    let reference_b = &input.families[b].reference;
    let assessment = assess(reference_a, 0, reference_b, cap);
    let sequence = (assessment.profile.is_some() && assessment.paired > 0).then(|| SequenceMetrics {
        identity: assessment.matches as f64 / assessment.paired as f64,
        coverage_a: assessment.paired as f64 / reference_a.protein.len() as f64,
        coverage_b: assessment.paired as f64 / reference_b.protein.len() as f64,
    });
    let mut row = Candidate {
        a,
        b,
        union_samples: x.samples.len() + y.samples.len() - cooccurring.len(),
        cooccurring,
        distinct_locus,
        sequence_status: assessment.status,
        sequence,
        neighborhood: neighborhood(x, y),
        rejection_reasons: Vec::new(),
        simulated: false,
    };
    if row.union_samples < core_samples {
        row.rejection_reasons.push("union_below_core");
    }
    if !row.distinct_locus.is_empty() {
        row.rejection_reasons.push("distinct_locus_cooccurrence");
    }
    if !row.cooccurring.is_empty() {
        row.rejection_reasons.push("sample_sets_not_disjoint");
    }
    if x.multilocus_samples > 0 || y.multilocus_samples > 0 {
        row.rejection_reasons.push("within_family_multiple_loci");
    }
    if !x.uncertain_samples.is_empty() || !y.uncertain_samples.is_empty() {
        row.rejection_reasons.push("uncertain_physical_locus");
    }
    if !x.annotated_reference || !y.annotated_reference {
        row.rejection_reasons.push("reference_not_intact_annotated");
    }
    if reference_a.protein.len().min(reference_b.protein.len()) < MIN_PROTEIN_AA {
        row.rejection_reasons.push("short_reference");
    }
    match sequence {
        None => row.rejection_reasons.push("sequence_unassessed"),
        Some(metrics) => {
            if metrics.identity < MIN_IDENTITY {
                row.rejection_reasons.push("identity_below_threshold");
            }
            if metrics.coverage_a.min(metrics.coverage_b) < MIN_COVERAGE {
                row.rejection_reasons.push("coverage_below_threshold");
            }
        }
    }
    if row.neighborhood.support_a.min(row.neighborhood.support_b) < MIN_CONTEXT_SAMPLES {
        row.rejection_reasons.push("insufficient_replicated_neighborhood");
    }
    if row.neighborhood.agreement_a.min(row.neighborhood.agreement_b) < MIN_AGREEMENT {
        row.rejection_reasons.push("neighborhood_disagreement");
    }
    row
}

fn rank_candidates(rows: &mut [Candidate], families: usize) {
    let mut partners = vec![0usize; families];
    for row in rows.iter().filter(|row| row.strong()) {
        partners[row.a] += 1;
        partners[row.b] += 1;
    }
    for row in rows.iter_mut() {
        if row.strong() && (partners[row.a] > 1 || partners[row.b] > 1) {
            row.rejection_reasons.push("competing_candidate");
        }
    }
    rows.sort_by(|a, b| {
        b.strong()
            .cmp(&a.strong())
            .then_with(|| a.distinct_locus.len().cmp(&b.distinct_locus.len()))
            .then_with(|| a.cooccurring.len().cmp(&b.cooccurring.len()))
            .then_with(|| b.union_samples.cmp(&a.union_samples))
            .then_with(|| b.sequence_score().total_cmp(&a.sequence_score()))
            .then_with(|| {
                b.neighborhood.agreement_a.min(b.neighborhood.agreement_b).total_cmp(
                    &a.neighborhood.agreement_a.min(a.neighborhood.agreement_b),
                )
            })
            .then_with(|| (a.a, a.b).cmp(&(b.a, b.b)))
    });
}

fn candidate_pairs(evidence: &[Evidence], core_samples: usize, near_samples: usize) -> Vec<(usize, usize)> {
    let mut order: Vec<_> = (0..evidence.len())
        .filter(|&index| evidence[index].samples.len() < core_samples)
        .collect();
    order.sort_by_key(|&index| (std::cmp::Reverse(evidence[index].samples.len()), index));
    let mut pairs = Vec::new();
    for (position, &a) in order.iter().enumerate() {
        let x = &evidence[a].samples;
        if x.len() * 2 < near_samples {
            break;
        }
        for &b in &order[position + 1..] {
            let y = &evidence[b].samples;
            if x.len() + y.len() < near_samples {
                break;
            }
            let shared = x.intersection(y).count();
            if x.len() + y.len() - shared >= near_samples {
                pairs.push((a.min(b), a.max(b)));
            }
        }
    }
    pairs.sort_unstable();
    pairs
}

fn write_reports(
    input: &Input,
    loci: &[Locus],
    evidence: &[Evidence],
    rows: &[Candidate],
    output: &Path,
    test_merges: bool,
) -> Result<()> {
    let mut out = writer(&output.join("candidates.tsv"))?;
    writeln!(out, "rank\tgene_a\tgene_b\tsamples_a\tsamples_b\tunion_samples\tunion_fraction\tcooccurring_samples\tdistinct_locus_cooccurring_samples\tprotein_aa_a\tprotein_aa_b\tsequence_status\tidentity_paired_pct\tcoverage_a_pct\tcoverage_b_pct\tupstream\tdownstream\tcontext_support_a\tcontext_support_b\tcomplete_context_samples_a\tcomplete_context_samples_b\tneighborhood_agreement_a\tneighborhood_agreement_b\tmultilocus_samples_a\tmultilocus_samples_b\tuncertain_locus_samples_a\tuncertain_locus_samples_b\tstrong_candidate\trejection_reasons\tmerge_test")?;
    let mut union = writer(&output.join("candidate_presence.tsv"))?;
    write!(union, "candidate")?;
    for sample in &input.samples {
        write!(union, "\t{sample}")?;
    }
    writeln!(union)?;
    let mut overlap = writer(&output.join("cooccurrence.tsv"))?;
    writeln!(overlap, "rank\tgene_a\tgene_b\tsample\tdistinct_locus\tloci_a\tloci_b")?;
    for (index, row) in rows.iter().enumerate() {
        let rank = index + 1;
        let (a, b) = (&evidence[row.a], &evidence[row.b]);
        let name_a = &input.families[row.a].reference.id;
        let name_b = &input.families[row.b].reference.id;
        let metric = |value: Option<f64>| value.map(|x| format!("{:.6}", x * 100.0)).unwrap_or_default();
        let context_name = |side: bool| {
            row.neighborhood.context.map(|(up, down)| {
                input.families[if side { down } else { up }].reference.id.as_str()
            }).unwrap_or("")
        };
        let state = if !row.strong() {
            "rejected"
        } else if !test_merges {
            "not_requested"
        } else if row.simulated {
            "simulated"
        } else {
            "limit_reached"
        };
        let fields = vec![
            rank.to_string(), name_a.clone(), name_b.clone(),
            a.samples.len().to_string(), b.samples.len().to_string(),
            row.union_samples.to_string(),
            format!("{:.6}", row.union_samples as f64 / input.samples.len() as f64),
            row.cooccurring.len().to_string(), row.distinct_locus.len().to_string(),
            input.families[row.a].reference.protein.len().to_string(),
            input.families[row.b].reference.protein.len().to_string(),
            row.sequence_status.to_string(),
            metric(row.sequence.map(|x| x.identity)),
            metric(row.sequence.map(|x| x.coverage_a)),
            metric(row.sequence.map(|x| x.coverage_b)),
            context_name(false).to_string(), context_name(true).to_string(),
            row.neighborhood.support_a.to_string(), row.neighborhood.support_b.to_string(),
            a.complete_context_samples.len().to_string(), b.complete_context_samples.len().to_string(),
            format!("{:.6}", row.neighborhood.agreement_a),
            format!("{:.6}", row.neighborhood.agreement_b),
            a.multilocus_samples.to_string(), b.multilocus_samples.to_string(),
            a.uncertain_samples.len().to_string(), b.uncertain_samples.len().to_string(),
            row.strong().to_string(), row.rejection_reasons.join(";"), state.to_string(),
        ];
        writeln!(out, "{}", fields.join("\t"))?;
        write!(union, "pair_{rank:05}")?;
        for sample in 0..input.samples.len() {
            write!(union, "\t{}", u8::from(a.samples.contains(&sample) || b.samples.contains(&sample)))?;
        }
        writeln!(union)?;
        for &sample in &row.cooccurring {
            let descriptions = |family: &Evidence| {
                let values: Vec<_> = family.loci_by_sample[&sample].iter().map(|&index| {
                    let locus = &loci[index];
                    format!("{}:{}:{}-{}:{}", input.contigs[locus.contig], locus.id, locus.start, locus.end, locus.strand)
                }).collect();
                list(values.iter().map(String::as_str))
            };
            writeln!(overlap, "{rank}\t{name_a}\t{name_b}\t{}\t{}\t{}\t{}",
                input.samples[sample], row.distinct_locus.contains(&sample),
                descriptions(a), descriptions(b))?;
        }
    }
    out.flush()?;
    union.flush()?;
    overlap.flush()?;
    Ok(())
}

fn simulate(
    input: &Input,
    evidence: &[Evidence],
    rows: &mut [Candidate],
    output: &Path,
    core_samples: usize,
    max_merges: usize,
) -> Result<(usize, usize)> {
    let mut assigned: Vec<_> = (0..input.families.len()).collect();
    let mut merged_presence: Vec<_> = evidence.iter().map(|row| row.samples.clone()).collect();
    let mut selected_ranks = BTreeMap::new();
    let mut merges = 0;
    for (index, row) in rows.iter_mut().enumerate() {
        if !row.strong() || merges == max_merges {
            continue;
        }
        if assigned[row.a] != row.a || assigned[row.b] != row.b
            || selected_ranks.contains_key(&row.a) || selected_ranks.contains_key(&row.b)
            || !merged_presence[row.a].is_disjoint(&merged_presence[row.b])
        {
            return Err(Error::new("near-core simulation safety invariant failed: reused or co-occurring families"));
        }
        let presence_b = std::mem::take(&mut merged_presence[row.b]);
        merged_presence[row.a].extend(presence_b);
        assigned[row.b] = row.a;
        selected_ranks.insert(row.a, index + 1);
        selected_ranks.insert(row.b, index + 1);
        row.simulated = true;
        merges += 1;
    }
    let mut mapping = writer(&output.join("simulated_family_map.tsv"))?;
    writeln!(mapping, "original_gene\tsimulated_gene\tcandidate_rank")?;
    for (index, &target) in assigned.iter().enumerate() {
        let rank = selected_ranks.get(&index).map(usize::to_string).unwrap_or_default();
        writeln!(mapping, "{}\t{}\t{rank}",
            input.families[index].reference.id, input.families[target].reference.id)?;
    }
    mapping.flush()?;
    let mut presence = writer(&output.join("simulated_presence.tsv"))?;
    write!(presence, "gene")?;
    for sample in &input.samples {
        write!(presence, "\t{sample}")?;
    }
    writeln!(presence)?;
    let mut core_after = 0;
    for (index, samples) in merged_presence.iter().enumerate() {
        if assigned[index] != index {
            continue;
        }
        core_after += usize::from(samples.len() >= core_samples);
        write!(presence, "{}", input.families[index].reference.id)?;
        for sample in 0..input.samples.len() {
            write!(presence, "\t{}", u8::from(samples.contains(&sample)))?;
        }
        writeln!(presence)?;
    }
    presence.flush()?;
    Ok((merges, core_after))
}

/// Diagnose pairwise near-core splitting; optional merges are a presence-only counterfactual.
pub fn near_core_splits(
    directory: &Path,
    output: &Path,
    options: &NearCoreOptions,
) -> Result<NearCoreSummary> {
    if !options.core_fraction.is_finite() || !options.near_core_fraction.is_finite()
        || options.near_core_fraction <= 0.0 || options.core_fraction > 1.0
        || options.near_core_fraction > options.core_fraction
    {
        return Err(Error::new("require 0 < near-core-fraction <= core-fraction <= 1"));
    }
    if options.max_merges == 0 {
        return Err(Error::new("max-merges must be positive"));
    }
    if output.exists() {
        return Err(Error::new("near-core output already exists; use a fresh directory"));
    }
    let mut fingerprints = Vec::new();
    for name in ["genes.fa", "presence.tsv", "members.tsv", "member_alleles.fa"] {
        let path = directory.join(name);
        if !path.is_file() {
            return Err(Error::new(format!(
                "missing near-core input {}; build with --audit-members", path.display()
            )));
        }
        fingerprints.push((name, fingerprint(&path)?));
    }
    let pool = rayon::ThreadPoolBuilder::new().num_threads(options.threads).build()?;
    let input = load_in(directory, Some(&pool))?;
    let loci = loci(&input)?;
    let evidence = family_evidence(&input, &loci);
    let core_samples = (options.core_fraction * input.samples.len() as f64).ceil() as usize;
    let near_samples = (options.near_core_fraction * input.samples.len() as f64).ceil() as usize;
    if let Some(parent) = output.parent().filter(|path| !path.as_os_str().is_empty()) {
        fs::create_dir_all(parent)?;
    }
    fs::create_dir(output)?;
    log::open(&output.join("rustroo.log"))?;
    let _log_guard = LogGuard;
    progress!("rustroo near-core-splits: enumerating sub-core pairs by sample union; input assignments remain unchanged");
    let pairs = candidate_pairs(&evidence, core_samples, near_samples);
    progress!("rustroo near-core-splits: aligning {} candidate reference pairs on {} threads", pairs.len(), pool.current_num_threads());
    let mut rows: Vec<_> = pool.install(|| {
        pairs.par_iter().map(|&(a, b)| {
            candidate(a, b, &input, &loci, &evidence, core_samples, options.max_alignment_cells)
        }).collect()
    });
    rank_candidates(&mut rows, input.families.len());
    let core_before = evidence.iter().filter(|row| row.samples.len() >= core_samples).count();
    let (simulated_merges, core_after) = if options.test_merges {
        let (merges, after) = simulate(
            &input, &evidence, &mut rows, output, core_samples, options.max_merges,
        )?;
        (merges, Some(after))
    } else {
        (0, None)
    };
    write_reports(&input, &loci, &evidence, &rows, output, options.test_merges)?;
    let mut provenance = writer(&output.join("input_fingerprints.tsv"))?;
    writeln!(provenance, "file\tfnv1a64_noncryptographic")?;
    for (name, before) in fingerprints {
        if fingerprint(&directory.join(name))? != before {
            return Err(Error::new(format!(
                "near-core input {name} changed during analysis; reports are incomplete"
            )));
        }
        writeln!(provenance, "{name}\t{before:016x}")?;
    }
    provenance.flush()?;
    let summary = NearCoreSummary {
        candidate_pairs: rows.len(),
        strong_pairs: rows.iter().filter(|row| row.strong()).count(),
        unassessed_pairs: rows.iter().filter(|row| row.sequence.is_none()).count(),
        simulated_merges,
        core_before,
        core_after,
    };
    let mut out = writer(&output.join("summary.tsv"))?;
    writeln!(out, "metric\tvalue")?;
    for (name, value) in [
        ("samples", input.samples.len()),
        ("families_before", input.families.len()),
        ("core_min_samples", core_samples),
        ("near_core_min_samples", near_samples),
        ("candidate_pairs", summary.candidate_pairs),
        ("strong_pairs", summary.strong_pairs),
        ("unassessed_pairs", summary.unassessed_pairs),
        ("pairs_with_distinct_locus_cooccurrence", rows.iter().filter(|row| !row.distinct_locus.is_empty()).count()),
        ("core_before", summary.core_before),
        ("simulated_merges", summary.simulated_merges),
        ("max_merges", options.max_merges),
        ("threads", pool.current_num_threads()),
        ("max_alignment_cells", options.max_alignment_cells),
        ("min_context_samples_per_family", MIN_CONTEXT_SAMPLES),
        ("min_reference_protein_aa", MIN_PROTEIN_AA),
    ] {
        writeln!(out, "{name}\t{value}")?;
    }
    writeln!(out, "core_fraction\t{}\nnear_core_fraction\t{}\nmin_identity\t{MIN_IDENTITY}\nmin_coverage_each_reference\t{MIN_COVERAGE}\nmin_neighborhood_agreement_each_family\t{MIN_AGREEMENT}\nmerge_test_requested\t{}",
        options.core_fraction, options.near_core_fraction, options.test_merges)?;
    let optional_count = |value: Option<usize>| value.map(|x| x.to_string()).unwrap_or_else(|| "not_assessed".into());
    let unsafe_merges = rows.iter().filter(|row| row.simulated && !row.distinct_locus.is_empty()).count();
    writeln!(out, "core_after_simulation\t{}\nfamilies_after_simulation\t{}\ncore_gain\t{}\nsimulated_distinct_locus_cooccurring_merges\t{}",
        optional_count(summary.core_after),
        optional_count(summary.core_after.map(|_| input.families.len() - simulated_merges)),
        optional_count(summary.core_after.map(|after| after - core_before)),
        optional_count(summary.core_after.map(|_| unsafe_merges)))?;
    writeln!(out, "comparison_core_count\t{}", optional_count(options.comparison_core_count))?;
    if let Some(target) = options.comparison_core_count {
        let gap = target.saturating_sub(core_before);
        writeln!(out, "comparison_gap_before\t{gap}\ncomparison_gap_after_simulation\t{}\ncomparison_surplus_after_simulation\t{}",
            optional_count(core_after.map(|after| target.saturating_sub(after))),
            optional_count(core_after.map(|after| after.saturating_sub(target))))?;
        let closed = core_after.filter(|_| gap > 0)
            .map(|after| format!("{:.6}", (after - core_before).min(gap) as f64 / gap as f64))
            .unwrap_or_else(|| "not_assessed".into());
        writeln!(out, "comparison_gap_fraction_closed\t{closed}")?;
    }
    let mut rejections = BTreeMap::new();
    for reason in rows.iter().flat_map(|row| row.rejection_reasons.iter().copied()) {
        *rejections.entry(reason).or_insert(0usize) += 1;
    }
    for (reason, count) in rejections {
        writeln!(out, "rejected_pairs_{reason}\t{count}")?;
    }
    writeln!(out, "method_version\t1\nrustroo_version\t{}\nprotein_alignment_scoring\tmatch=2;mismatch=-1;gap_open=-5;gap_extend=-1\nneighborhood_gap_bp\t-15..15000\nassignments_modified\tfalse\nsimulation_scope\tpresence_only_no_references_or_memberships_rewritten\ncomparison_cohort_and_threshold\tuser_must_match\ncomplete\ttrue", env!("CARGO_PKG_VERSION"))?;
    out.flush()?;
    if summary.unassessed_pairs > 0 {
        progress!("rustroo near-core-splits: {} sequence comparisons unassessed; inspect sequence_status, no heuristic replacements were used", summary.unassessed_pairs);
    }
    progress!("rustroo near-core-splits: complete; {} candidates, {} strong pairs, {} simulated merges; original clustering unchanged",
        summary.candidate_pairs, summary.strong_pairs, summary.simulated_merges);
    Ok(summary)
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::{Family, Member, Sequence};
    use crate::genetic_code::back_translate;
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicUsize, Ordering};

    static NEXT: AtomicUsize = AtomicUsize::new(0);

    struct Directory(PathBuf);

    impl Directory {
        fn new() -> Self {
            let path = std::env::temp_dir().join(format!(
                "rustroo-near-core-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed),
            ));
            fs::create_dir(&path).unwrap();
            Self(path)
        }
    }

    impl Drop for Directory {
        fn drop(&mut self) {
            fs::remove_dir_all(&self.0).expect("remove only this test's directory");
        }
    }

    fn add_locus(
        input: &mut Input,
        family: usize,
        sample: usize,
        contig: usize,
        start: u32,
        strand: i8,
    ) {
        let index = input.members.len();
        input.members.push(Member {
            family,
            sample,
            contig,
            annotation: format!("annotation_{index}"),
            start,
            end: start + input.alleles[0].dna.len() as u32 - 1,
            strand,
            locus: index as u32,
            origin: "annotated",
            allele: 0,
        });
        input.families[family].members.push(index);
        *input.families[family].alleles.entry(0).or_default() += 1;
    }

    fn fixture() -> Input {
        fixture_with_samples(6)
    }

    fn fixture_with_samples(samples: usize) -> Input {
        let dna = back_translate(b"MACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY");
        let mut input = Input {
            families: ["a", "b", "left", "right"].iter().map(|name| Family {
                reference: Sequence::new((*name).into(), dna.clone()),
                members: Vec::new(),
                alleles: BTreeMap::new(),
            }).collect(),
            alleles: vec![Sequence::new("dna_0".into(), dna)],
            aliases: BTreeMap::from([("dna_0".into(), 0)]),
            samples: (0..samples).map(|index| format!("sample_{index}")).collect(),
            contigs: (0..samples).map(|index| format!("contig_{index}")).collect(),
            members: Vec::new(),
        };
        for sample in 0..samples {
            add_locus(&mut input, 2, sample, sample, 1, 1);
            add_locus(&mut input, usize::from(sample >= samples / 2), sample, sample, 500, 1);
            add_locus(&mut input, 3, sample, sample, 1000, 1);
        }
        input
    }

    fn analyzed(input: &Input, cap: usize) -> (Vec<Evidence>, Vec<Candidate>) {
        let loci = loci(input).unwrap();
        let evidence = family_evidence(input, &loci);
        let mut rows: Vec<_> = candidate_pairs(&evidence, 6, 6).into_iter()
            .map(|(a, b)| candidate(a, b, input, &loci, &evidence, 6, cap))
            .collect();
        rank_candidates(&mut rows, input.families.len());
        (evidence, rows)
    }

    fn replace_protein(input: &mut Input, family: usize, protein: &[u8]) {
        let dna = back_translate(protein);
        let allele = input.alleles.len();
        input.families[family].reference = Sequence::new(
            input.families[family].reference.id.clone(), dna.clone(),
        );
        input.alleles.push(Sequence::new(format!("dna_{allele}"), dna.clone()));
        input.families[family].alleles.clear();
        for &member in &input.families[family].members {
            input.members[member].allele = allele;
            input.members[member].end = input.members[member].start + dna.len() as u32 - 1;
        }
        let count = input.families[family].members.len();
        input.families[family].alleles.insert(allele, count);
    }

    fn write_input(input: &Input, path: &Path) {
        fs::create_dir(path).unwrap();
        let genes: String = input.families.iter().map(|family| {
            format!(">{}\n{}\n", family.reference.id, family.reference.dna)
        }).collect();
        let alleles: String = input.alleles.iter().map(|allele| {
            format!(">{}\n{}\n", allele.id, allele.dna)
        }).collect();
        fs::write(path.join("genes.fa"), genes).unwrap();
        fs::write(path.join("member_alleles.fa"), alleles).unwrap();
        let mut members = "gene\tsample\tcontig\tannotation_id\tstart\tend\tstrand\tlocus\torigin\tallele\n".to_string();
        for member in &input.members {
            members.push_str(&format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                input.families[member.family].reference.id,
                input.samples[member.sample], input.contigs[member.contig], member.annotation,
                member.start, member.end, member.strand, member.locus, member.origin,
                input.alleles[member.allele].id));
        }
        fs::write(path.join("members.tsv"), members).unwrap();
        let mut presence = format!("gene\t{}\n", input.samples.join("\t"));
        for family in &input.families {
            presence.push_str(&family.reference.id);
            for sample in 0..input.samples.len() {
                presence.push_str(if family.members.iter().any(|&member| {
                    input.members[member].sample == sample
                }) { "\t1" } else { "\t0" });
            }
            presence.push('\n');
        }
        fs::write(path.join("presence.tsv"), presence).unwrap();
    }

    fn options() -> NearCoreOptions {
        NearCoreOptions { threads: 1, ..NearCoreOptions::default() }
    }

    fn summary_table(path: &Path) -> BTreeMap<String, String> {
        fs::read_to_string(path).unwrap().lines().skip(1).map(|line| {
            let (name, value) = line.split_once('\t').unwrap();
            (name.to_string(), value.to_string())
        }).collect()
    }

    #[test]
    fn disjoint_complements_have_full_union_and_replicated_context() {
        let input = fixture();
        let (_, rows) = analyzed(&input, 0);
        assert_eq!(rows.len(), 1);
        let row = &rows[0];
        assert!(row.strong(), "{:?}", row.rejection_reasons);
        assert_eq!((row.a, row.b, row.union_samples), (0, 1, 6));
        assert!(row.cooccurring.is_empty());
        assert_eq!((row.neighborhood.support_a, row.neighborhood.support_b), (3, 3));
        let sequence = row.sequence.unwrap();
        assert_eq!((sequence.identity, sequence.coverage_a, sequence.coverage_b), (1.0, 1.0, 1.0));
    }

    #[test]
    fn reverse_strand_context_is_normalized() {
        let mut input = fixture();
        for member in input.members.iter_mut().filter(|member| member.sample >= 3) {
            member.strand = -1;
            if member.family == 2 || member.family == 3 {
                let length = member.end - member.start;
                member.start = if member.family == 2 { 1000 } else { 1 };
                member.end = member.start + length;
            }
        }
        let (_, rows) = analyzed(&input, 0);
        assert!(rows[0].strong(), "{:?}", rows[0].rejection_reasons);
        assert_eq!(rows[0].neighborhood.context, Some((2, 3)));
    }

    #[test]
    fn distinct_and_overlapping_cooccurrence_are_both_vetoed() {
        for same_locus_coordinates in [false, true] {
            let mut input = fixture();
            let contig = if same_locus_coordinates {
                0
            } else {
                input.contigs.push("another_contig".into());
                input.contigs.len() - 1
            };
            add_locus(&mut input, 1, 0, contig, 500, 1);
            let (_, rows) = analyzed(&input, 0);
            assert_eq!(rows[0].cooccurring, vec![0]);
            assert_eq!(rows[0].distinct_locus.is_empty(), same_locus_coordinates);
            assert!(rows[0].rejection_reasons.contains(&"sample_sets_not_disjoint"));
            assert!(!rows[0].strong());
        }
    }

    #[test]
    fn missing_context_does_not_inflate_agreement_and_refound_context_is_excluded() {
        let mut input = fixture();
        let mut loci = loci(&input).unwrap();
        for locus in loci.iter_mut().filter(|locus| locus.family == 0 && locus.sample != 0) {
            locus.flanks = (None, None);
        }
        let evidence = family_evidence(&input, &loci);
        let row = candidate(0, 1, &input, &loci, &evidence, 6, 0);
        assert_eq!(row.neighborhood.support_a, 1);
        assert_eq!(row.neighborhood.agreement_a, 1.0 / 3.0);
        assert!(row.rejection_reasons.contains(&"neighborhood_disagreement"));
        for member in input.members.iter_mut().filter(|member| member.family == 0) {
            member.origin = "refound";
            member.annotation.clear();
        }
        let (_, rows) = analyzed(&input, 0);
        assert_eq!(rows[0].neighborhood.support_a, 0);
        assert!(rows[0].rejection_reasons.contains(&"reference_not_intact_annotated"));
    }

    #[test]
    fn within_family_multiple_loci_are_not_merged() {
        let mut input = fixture();
        add_locus(&mut input, 0, 0, 0, 2000, 1);
        let (evidence, rows) = analyzed(&input, 0);
        assert_eq!(evidence[0].multilocus_samples, 1);
        assert!(rows[0].rejection_reasons.contains(&"within_family_multiple_loci"));
        let mut input = fixture();
        input.members[1].end += 1;
        let (_, rows) = analyzed(&input, 0);
        assert!(rows[0].rejection_reasons.contains(&"uncertain_physical_locus"));
    }

    #[test]
    fn neighborhood_requires_ninety_percent_and_three_samples_per_side() {
        let input = fixture_with_samples(20);
        for support in [9, 8] {
            let mut loci = loci(&input).unwrap();
            for locus in &mut loci {
                if locus.family < 2 && locus.sample % 10 >= support {
                    locus.flanks = (None, None);
                }
            }
            let evidence = family_evidence(&input, &loci);
            let row = candidate(0, 1, &input, &loci, &evidence, 20, 0);
            assert_eq!(row.neighborhood.agreement_a, support as f64 / 10.0);
            assert_eq!(row.neighborhood.agreement_b, support as f64 / 10.0);
            assert_eq!(row.strong(), support == 9, "{:?}", row.rejection_reasons);
        }
        let input = fixture_with_samples(4);
        let loci = loci(&input).unwrap();
        let evidence = family_evidence(&input, &loci);
        let row = candidate(0, 1, &input, &loci, &evidence, 4, 0);
        assert_eq!(row.neighborhood.agreement_a, 1.0);
        assert_eq!(row.rejection_reasons, vec!["insufficient_replicated_neighborhood"]);
    }

    #[test]
    fn identity_and_bidirectional_coverage_thresholds_are_inclusive() {
        let mut input = fixture();
        let mut protein = vec![b'A'; 100];
        protein[0] = b'M';
        replace_protein(&mut input, 0, &protein);
        for mismatches in [2, 3] {
            let mut other = protein.clone();
            for index in 0..mismatches {
                other[20 + index * 10] = b'W';
            }
            replace_protein(&mut input, 1, &other);
            let (_, rows) = analyzed(&input, 0);
            assert_eq!(rows[0].strong(), mismatches == 2, "{:?}", rows[0].rejection_reasons);
            assert_eq!(rows[0].sequence.unwrap().identity, (100 - mismatches) as f64 / 100.0);
        }
        for length in [95, 94] {
            replace_protein(&mut input, 1, &protein[..length]);
            let (_, rows) = analyzed(&input, 0);
            assert_eq!(rows[0].strong(), length == 95, "{:?}", rows[0].rejection_reasons);
            let sequence = rows[0].sequence.unwrap();
            assert_eq!(sequence.coverage_a, length as f64 / 100.0);
            assert_eq!(sequence.coverage_b, 1.0);
        }
    }

    #[test]
    fn work_limited_pairs_remain_explicitly_unassessed() {
        let mut input = fixture();
        replace_protein(&mut input, 1, b"MACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWA");
        let (_, rows) = analyzed(&input, 1);
        assert_eq!(rows[0].sequence_status, "alignment_work_limit");
        assert!(rows[0].sequence.is_none());
        assert!(rows[0].rejection_reasons.contains(&"sequence_unassessed"));
    }

    #[test]
    fn competing_partners_are_not_resolved_by_rank_or_transitive_chaining() {
        let mut input = fixture();
        let dna = input.alleles[0].dna.clone();
        input.families.push(Family {
            reference: Sequence::new("third".into(), dna),
            members: Vec::new(),
            alleles: BTreeMap::new(),
        });
        for sample in 3..6 {
            let contig = input.contigs.len();
            input.contigs.push(format!("extra_{sample}"));
            add_locus(&mut input, 2, sample, contig, 1, 1);
            add_locus(&mut input, 4, sample, contig, 500, 1);
            add_locus(&mut input, 3, sample, contig, 1000, 1);
        }
        let (_, rows) = analyzed(&input, 0);
        assert_eq!(rows.len(), 2);
        assert!(rows.iter().all(|row| row.rejection_reasons == vec!["competing_candidate"]));
    }

    #[test]
    fn report_is_read_only_and_simulation_measures_core_gap() {
        let directory = Directory::new();
        let input = directory.0.join("input");
        write_input(&fixture(), &input);
        let names = ["genes.fa", "presence.tsv", "members.tsv", "member_alleles.fa"];
        let before: Vec<_> = names.iter().map(|name| fingerprint(&input.join(name)).unwrap()).collect();
        let output = directory.0.join("report");
        let summary = near_core_splits(&input, &output, &options()).unwrap();
        assert_eq!(summary.candidate_pairs, 1);
        assert_eq!(summary.strong_pairs, 1);
        assert_eq!(summary.core_before, 2);
        assert_eq!(summary.core_after, None);
        assert!(!output.join("simulated_presence.tsv").exists());
        let report = fs::read_to_string(output.join("candidates.tsv")).unwrap();
        let mut lines = report.lines();
        let header: Vec<_> = lines.next().unwrap().split('\t').collect();
        let values: Vec<_> = lines.next().unwrap().split('\t').collect();
        assert_eq!(header.len(), values.len());
        let fields: BTreeMap<_, _> = header.into_iter().zip(values).collect();
        assert_eq!(fields["merge_test"], "not_requested");
        assert_eq!(fields["union_fraction"], "1.000000");
        assert_eq!(fields["identity_paired_pct"], "100.000000");
        assert!(fs::read_to_string(output.join("candidate_presence.tsv")).unwrap()
            .contains("pair_00001\t1\t1\t1\t1\t1\t1\n"));
        let simulation = directory.0.join("simulation");
        let summary = near_core_splits(&input, &simulation, &NearCoreOptions {
            test_merges: true,
            comparison_core_count: Some(4),
            ..options()
        }).unwrap();
        assert_eq!(summary.simulated_merges, 1);
        assert_eq!(summary.core_after, Some(3));
        let metrics = summary_table(&simulation.join("summary.tsv"));
        assert_eq!(metrics["core_gain"], "1");
        assert_eq!(metrics["comparison_gap_before"], "2");
        assert_eq!(metrics["comparison_gap_after_simulation"], "1");
        assert_eq!(metrics["comparison_gap_fraction_closed"], "0.500000");
        assert_eq!(metrics["simulated_distinct_locus_cooccurring_merges"], "0");
        assert_eq!(metrics["complete"], "true");
        let mapping = fs::read_to_string(simulation.join("simulated_family_map.tsv")).unwrap();
        assert!(mapping.contains("a\ta\t1\n"));
        assert!(mapping.contains("b\ta\t1\n"));
        let presence = fs::read_to_string(simulation.join("simulated_presence.tsv")).unwrap();
        assert!(presence.contains("a\t1\t1\t1\t1\t1\t1\n"));
        assert!(!presence.lines().any(|line| line.starts_with("b\t")));
        assert!(near_core_splits(&input, &output, &options()).is_err());
        let after: Vec<_> = names.iter().map(|name| fingerprint(&input.join(name)).unwrap()).collect();
        assert_eq!(before, after);
    }

    #[test]
    fn simulation_retains_cooccurring_families_and_records_their_loci() {
        let directory = Directory::new();
        let mut data = fixture();
        add_locus(&mut data, 1, 0, 0, 1500, 1);
        let input = directory.0.join("input");
        write_input(&data, &input);
        let output = directory.0.join("simulation");
        let summary = near_core_splits(&input, &output, &NearCoreOptions {
            test_merges: true, ..options()
        }).unwrap();
        assert_eq!(summary.simulated_merges, 0);
        assert_eq!(summary.core_after, Some(summary.core_before));
        let overlap = fs::read_to_string(output.join("cooccurrence.tsv")).unwrap();
        assert!(overlap.contains("1\ta\tb\tsample_0\ttrue\t"));
        let mapping = fs::read_to_string(output.join("simulated_family_map.tsv")).unwrap();
        assert!(mapping.contains("b\tb\t\n"));
    }

    #[test]
    fn simulation_respects_the_cap_without_combining_independent_candidates() {
        let directory = Directory::new();
        let mut data = fixture();
        let dna = data.alleles[0].dna.clone();
        for name in ["c", "d", "second_left", "second_right"] {
            data.families.push(Family {
                reference: Sequence::new(name.into(), dna.clone()),
                members: Vec::new(),
                alleles: BTreeMap::new(),
            });
        }
        for sample in 0..6 {
            let contig = data.contigs.len();
            data.contigs.push(format!("second_{sample}"));
            add_locus(&mut data, 6, sample, contig, 1, 1);
            add_locus(&mut data, if sample < 3 { 4 } else { 5 }, sample, contig, 500, 1);
            add_locus(&mut data, 7, sample, contig, 1000, 1);
        }
        let input = directory.0.join("input");
        write_input(&data, &input);
        let output = directory.0.join("simulation");
        let summary = near_core_splits(&input, &output, &NearCoreOptions {
            test_merges: true, max_merges: 1, ..options()
        }).unwrap();
        assert_eq!(summary.strong_pairs, 2);
        assert_eq!(summary.simulated_merges, 1);
        assert_eq!(summary.core_after, Some(summary.core_before + 1));
        let report = fs::read_to_string(output.join("candidates.tsv")).unwrap();
        assert_eq!(report.lines().filter(|line| line.ends_with("\tsimulated")).count(), 1);
        assert_eq!(report.lines().filter(|line| line.ends_with("\tlimit_reached")).count(), 1);
    }

    #[test]
    fn empty_candidate_universe_keeps_the_simulated_matrix_unchanged() {
        let directory = Directory::new();
        let input = directory.0.join("input");
        write_input(&fixture(), &input);
        let output = directory.0.join("simulation");
        let summary = near_core_splits(&input, &output, &NearCoreOptions {
            core_fraction: 0.5,
            near_core_fraction: 0.5,
            test_merges: true,
            ..options()
        }).unwrap();
        assert_eq!(summary.candidate_pairs, 0);
        assert_eq!(summary.core_before, 4);
        assert_eq!(summary.core_after, Some(4));
        assert_eq!(
            fs::read_to_string(input.join("presence.tsv")).unwrap(),
            fs::read_to_string(output.join("simulated_presence.tsv")).unwrap(),
        );
    }

    #[test]
    fn pair_enumeration_is_complete_and_does_not_require_shared_context() {
        let mut evidence: Vec<_> = (0..5).map(|_| Evidence::default()).collect();
        for (row, samples) in evidence.iter_mut().zip([
            vec![0, 1, 2, 3, 4], vec![4, 5], vec![0, 1, 2], vec![3, 4, 5], vec![0, 1, 2, 3, 4, 5],
        ]) {
            row.samples.extend(samples);
        }
        let actual: BTreeSet<_> = candidate_pairs(&evidence, 6, 5).into_iter().collect();
        let mut expected = BTreeSet::new();
        for a in 0..evidence.len() {
            for b in a + 1..evidence.len() {
                if evidence[a].samples.len() < 6 && evidence[b].samples.len() < 6
                    && evidence[a].samples.union(&evidence[b].samples).count() >= 5
                {
                    expected.insert((a, b));
                }
            }
        }
        assert_eq!(actual, expected);
        assert!(actual.contains(&(1, 2)));
    }

    #[test]
    fn invalid_options_fail_before_creating_output() {
        let directory = Directory::new();
        for (index, (near, core)) in [
            (0.0, 0.99), (0.99, 0.95), (0.95, 1.1),
            (f64::NAN, 0.99), (0.95, f64::INFINITY),
        ].into_iter().enumerate() {
            let output = directory.0.join(format!("invalid_{index}"));
            assert!(near_core_splits(&directory.0, &output, &NearCoreOptions {
                near_core_fraction: near, core_fraction: core, ..options()
            }).is_err());
            assert!(!output.exists());
        }
        assert!(near_core_splits(&directory.0, &directory.0.join("zero"), &NearCoreOptions {
            max_merges: 0, ..options()
        }).is_err());
    }
}
