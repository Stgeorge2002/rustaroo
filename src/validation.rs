use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

use rayon::prelude::*;

use crate::align::length_ratio;
use crate::genetic_code::{orf_status, translate, OrfStatus};
use crate::log::{self, progress};
use crate::phylogeny::{bootstrap, common_sites, distances, neighbor_joining, variable_sites};
use crate::validation_align::{align, Profile, Step};
use crate::{Error, Result};

#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum TreeMode {
    Uncertain,
    All,
    None,
}

#[derive(Clone, Debug)]
pub struct ValidationOptions {
    pub threads: usize,
    pub tree_mode: TreeMode,
    pub bootstrap_replicates: usize,
    pub max_tree_alleles: usize,
    pub max_alignment_cells: usize,
    pub min_tree_sites: usize,
}

impl Default for ValidationOptions {
    fn default() -> Self {
        Self {
            threads: 0,
            tree_mode: TreeMode::Uncertain,
            bootstrap_replicates: 100,
            max_tree_alleles: 128,
            max_alignment_cells: 5_000_000,
            min_tree_sites: 90,
        }
    }
}

#[derive(Debug)]
pub struct ValidationSummary {
    pub families: usize,
    pub unique_member_alleles: usize,
    pub physical_loci: usize,
    pub families_with_multiple_loci_per_sample: usize,
    pub unassessed_alignments: usize,
    pub unassessed_related_pairs: usize,
    pub trees_written: usize,
}

struct Sequence {
    id: String,
    dna: String,
    protein: Vec<u8>,
    orf: OrfStatus,
}

impl Sequence {
    fn new(id: String, dna: String) -> Self {
        let protein = translate(&dna);
        let orf = orf_status(&dna);
        Self {
            id,
            dna,
            protein,
            orf,
        }
    }
}

struct Family {
    reference: Sequence,
    members: Vec<usize>,
    alleles: BTreeMap<usize, usize>,
}

struct Member {
    family: usize,
    sample: usize,
    contig: usize,
    annotation: String,
    start: u32,
    end: u32,
    strand: i8,
    locus: u32,
    origin: &'static str,
    allele: usize,
}

struct Input {
    families: Vec<Family>,
    alleles: Vec<Sequence>,
    aliases: BTreeMap<String, usize>,
    samples: Vec<String>,
    contigs: Vec<String>,
    members: Vec<Member>,
}

fn fasta(path: &Path) -> Result<Vec<Sequence>> {
    let mut records = BTreeMap::new();
    let mut id: Option<String> = None;
    let mut dna = String::new();
    let finish = |id: String, dna: String, records: &mut BTreeMap<String, String>| -> Result<()> {
        if dna.is_empty() || records.insert(id.clone(), dna).is_some() {
            return Err(Error::new(format!(
                "{}: empty sequence or duplicate FASTA ID {id}",
                path.display()
            )));
        }
        Ok(())
    };
    for (line_number, line) in BufReader::new(File::open(path)?).lines().enumerate() {
        let line = line?;
        if let Some(header) = line.strip_prefix('>') {
            if let Some(previous) = id.take() {
                finish(previous, std::mem::take(&mut dna), &mut records)?;
            }
            if header.is_empty()
                || header
                    .chars()
                    .any(|character| character.is_whitespace() || character.is_control())
            {
                return Err(Error::new(format!(
                    "{}:{}: expected a nonempty FASTA ID without whitespace",
                    path.display(),
                    line_number + 1
                )));
            }
            id = Some(header.to_string());
        } else if !line.is_empty() {
            if id.is_none() {
                return Err(Error::new(format!(
                    "{}: sequence before FASTA header",
                    path.display()
                )));
            }
            for byte in line.bytes() {
                let base = byte.to_ascii_uppercase();
                if !b"ACGTURYSWKMBDHVN".contains(&base) {
                    return Err(Error::new(format!(
                        "{}:{}: invalid DNA character",
                        path.display(),
                        line_number + 1
                    )));
                }
                dna.push(if base == b'U' { 'T' } else { char::from(base) });
            }
        }
    }
    if let Some(last) = id {
        finish(last, dna, &mut records)?;
    }
    if records.is_empty() {
        return Err(Error::new(format!("{}: no FASTA records", path.display())));
    }
    let records: Vec<_> = records.into_iter().collect();
    Ok(records
        .into_par_iter()
        .map(|(id, dna)| Sequence::new(id, dna))
        .collect())
}

fn number<T: std::str::FromStr>(value: &str, field: &str, line: usize) -> Result<T> {
    value
        .parse()
        .map_err(|_| Error::new(format!("members.tsv:{line}: invalid {field}: {value}")))
}

fn load(directory: &Path) -> Result<Input> {
    let families: Vec<_> = fasta(&directory.join("genes.fa"))?
        .into_iter()
        .map(|reference| Family {
            reference,
            members: Vec::new(),
            alleles: BTreeMap::new(),
        })
        .collect();
    let family_of: HashMap<_, _> = families
        .iter()
        .enumerate()
        .map(|(index, family)| (family.reference.id.clone(), index))
        .collect();
    let mut alleles: Vec<Sequence> = Vec::new();
    let mut aliases = BTreeMap::new();
    let mut dna_of: HashMap<u64, Vec<usize>> = HashMap::new();
    for sequence in fasta(&directory.join("member_alleles.fa"))? {
        let hash = crate::model::content_hash(sequence.dna.as_bytes());
        if let Some(bucket) = dna_of.get(&hash) {
            if let Some(&index) = bucket
                .iter()
                .find(|&&index| alleles[index].dna == sequence.dna)
            {
                aliases.insert(sequence.id, index);
                continue;
            }
        }
        let index = alleles.len();
        dna_of.entry(hash).or_default().push(index);
        aliases.insert(sequence.id.clone(), index);
        alleles.push(sequence);
    }
    let mut presence_lines = BufReader::new(File::open(directory.join("presence.tsv"))?).lines();
    let header = presence_lines
        .next()
        .ok_or_else(|| Error::new("presence.tsv is empty"))??;
    let mut columns = header.split('\t');
    if columns.next() != Some("gene") {
        return Err(Error::new("invalid presence.tsv header"));
    }
    let samples: Vec<_> = columns.map(str::to_string).collect();
    if samples.is_empty()
        || samples
            .iter()
            .any(|name| name.is_empty() || name.chars().any(char::is_control))
        || samples.iter().collect::<HashSet<_>>().len() != samples.len()
    {
        return Err(Error::new(
            "presence.tsv requires unique nonempty sample names",
        ));
    }
    let sample_of: HashMap<_, _> = samples
        .iter()
        .enumerate()
        .map(|(index, name)| (name.as_str(), index))
        .collect();
    let mut expected = HashSet::new();
    let mut seen = HashSet::new();
    for line in presence_lines {
        let line = line?;
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() != samples.len() + 1 {
            return Err(Error::new("wrong presence.tsv column count"));
        }
        let family = *family_of
            .get(fields[0])
            .ok_or_else(|| Error::new("presence.tsv has an unknown gene"))?;
        if !seen.insert(family) {
            return Err(Error::new("duplicate presence.tsv gene"));
        }
        for (sample, value) in fields[1..].iter().enumerate() {
            match *value {
                "1" => {
                    expected.insert((family, sample));
                }
                "0" => {}
                _ => return Err(Error::new("presence.tsv values must be 0 or 1")),
            }
        }
    }
    if seen.len() != families.len() {
        return Err(Error::new("presence.tsv and genes.fa gene sets differ"));
    }
    let mut input = Input {
        families,
        alleles,
        aliases,
        samples: samples.clone(),
        contigs: Vec::new(),
        members: Vec::new(),
    };
    let mut contig_of = HashMap::new();
    let mut observed = HashSet::new();
    let mut unique_rows = HashSet::new();
    let mut lines = BufReader::new(File::open(directory.join("members.tsv"))?).lines();
    let header = lines
        .next()
        .ok_or_else(|| Error::new("members.tsv is empty"))??;
    if header != "gene\tsample\tcontig\tannotation_id\tstart\tend\tstrand\tlocus\torigin\tallele" {
        return Err(Error::new(
            "unsupported members.tsv header; regenerate with --audit-members",
        ));
    }
    for (offset, line) in lines.enumerate() {
        let line = line?;
        let line_number = offset + 2;
        if !unique_rows.insert(line.clone()) {
            return Err(Error::new(format!(
                "duplicate members.tsv row {line_number}"
            )));
        }
        let f: Vec<_> = line.split('\t').collect();
        if f.len() != 10 {
            return Err(Error::new(format!(
                "members.tsv:{line_number}: expected 10 columns"
            )));
        }
        if f.iter().any(|field| field.chars().any(char::is_control)) {
            return Err(Error::new(format!(
                "members.tsv:{line_number}: control character in a field"
            )));
        }
        let family = *family_of
            .get(f[0])
            .ok_or_else(|| Error::new(format!("members.tsv:{line_number}: unknown gene")))?;
        let sample = *sample_of
            .get(f[1])
            .ok_or_else(|| Error::new(format!("members.tsv:{line_number}: unknown sample")))?;
        if f[2].is_empty() {
            return Err(Error::new("empty member contig"));
        }
        let next_contig = input.contigs.len();
        let contig = *contig_of
            .entry((sample, f[2].to_string()))
            .or_insert_with(|| {
                input.contigs.push(f[2].to_string());
                next_contig
            });
        let start = number(f[4], "start", line_number)?;
        let end = number(f[5], "end", line_number)?;
        let strand = number(f[6], "strand", line_number)?;
        if start == 0 || end < start || !matches!(strand, -1 | 1) {
            return Err(Error::new(format!(
                "members.tsv:{line_number}: invalid coordinates or strand"
            )));
        }
        let origin = match f[8] {
            "annotated" => "annotated",
            "fusion_piece" => "fusion_piece",
            "refound" => "refound",
            _ => {
                return Err(Error::new(format!(
                    "members.tsv:{line_number}: unsupported origin"
                )))
            }
        };
        let allele = *input.aliases.get(f[9]).ok_or_else(|| {
            Error::new(format!(
                "members.tsv:{line_number}: missing allele {}",
                f[9]
            ))
        })?;
        let index = input.members.len();
        input.members.push(Member {
            family,
            sample,
            contig,
            annotation: f[3].to_string(),
            start,
            end,
            strand,
            locus: number(f[7], "locus", line_number)?,
            origin,
            allele,
        });
        input.families[family].members.push(index);
        *input.families[family].alleles.entry(allele).or_default() += 1;
        observed.insert((family, sample));
    }
    if observed != expected
        || input
            .families
            .iter()
            .any(|family| family.members.is_empty())
    {
        return Err(Error::new(
            "presence.tsv and members.tsv disagree; validation requires complete memberships",
        ));
    }
    Ok(input)
}

#[derive(Debug)]
struct Locus {
    family: usize,
    sample: usize,
    contig: usize,
    id: u32,
    start: u32,
    end: u32,
    strand: i8,
    members: Vec<usize>,
    flanks: (Option<usize>, Option<usize>),
}

fn loci(input: &Input) -> Result<Vec<Locus>> {
    let mut grouped: BTreeMap<_, Locus> = BTreeMap::new();
    for (index, member) in input.members.iter().enumerate() {
        let locus = grouped
            .entry((member.sample, member.contig, member.locus))
            .or_insert_with(|| Locus {
                family: member.family,
                sample: member.sample,
                contig: member.contig,
                id: member.locus,
                start: member.start,
                end: member.end,
                strand: member.strand,
                members: Vec::new(),
                flanks: (None, None),
            });
        if locus.family != member.family || locus.strand != member.strand {
            return Err(Error::new(
                "one exported physical locus has conflicting families or strands",
            ));
        }
        locus.start = locus.start.min(member.start);
        locus.end = locus.end.max(member.end);
        locus.members.push(index);
    }
    let mut result: Vec<_> = grouped.into_values().collect();
    result.sort_by_key(|locus| (locus.sample, locus.contig, locus.start, locus.end, locus.id));
    for position in 0..result.len() {
        let previous = position
            .checked_sub(1)
            .filter(|&neighbor| {
                result[neighbor].sample == result[position].sample
                    && result[neighbor].contig == result[position].contig
                    && (-15..=15_000).contains(
                        &(i64::from(result[position].start) - i64::from(result[neighbor].end) - 1),
                    )
            })
            .map(|neighbor| result[neighbor].family);
        let next = (position + 1 < result.len())
            .then_some(position + 1)
            .filter(|&neighbor| {
                result[neighbor].sample == result[position].sample
                    && result[neighbor].contig == result[position].contig
                    && (-15..=15_000).contains(
                        &(i64::from(result[neighbor].start) - i64::from(result[position].end) - 1),
                    )
            })
            .map(|neighbor| result[neighbor].family);
        result[position].flanks = if result[position].strand == 1 {
            (previous, next)
        } else {
            (next, previous)
        };
    }
    Ok(result)
}

fn writer(path: &Path) -> Result<BufWriter<File>> {
    Ok(BufWriter::with_capacity(
        1 << 20,
        OpenOptions::new().write(true).create_new(true).open(path)?,
    ))
}

fn escaped(value: &str) -> String {
    value
        .replace('%', "%25")
        .replace(';', "%3B")
        .replace(':', "%3A")
}

fn list<'a>(values: impl IntoIterator<Item = &'a str>) -> String {
    values
        .into_iter()
        .map(escaped)
        .collect::<Vec<_>>()
        .join(";")
}

fn bundle_evidence(locus: &Locus, input: &Input) -> &'static str {
    if locus.members.len() == 1 {
        return "single_model";
    }
    let mut members: Vec<_> = locus
        .members
        .iter()
        .map(|&index| &input.members[index])
        .collect();
    members.sort_by_key(|member| (member.start, member.end));
    for pair in members.windows(2) {
        let gap = i64::from(pair[1].start) - i64::from(pair[0].end) - 1;
        let shorter = (pair[0].end - pair[0].start + 1).min(pair[1].end - pair[1].start + 1);
        let overlap = pair[0]
            .end
            .min(pair[1].end)
            .saturating_sub(pair[0].start.max(pair[1].start))
            + 1;
        if gap > 3 || (gap < -3 && u64::from(overlap) * 100 < u64::from(shorter) * 95) {
            return "bundle_requires_review";
        }
        let reference_len = input.families[locus.family].reference.protein.len();
        if (-3..=3).contains(&gap)
            && pair.iter().all(|member| {
                input.alleles[member.allele].protein.len() * 100 >= reference_len * 95
            })
        {
            return "touching_full_length_models_requires_review";
        }
    }
    "overlapping_or_adjacent_models"
}

fn write_loci(input: &Input, loci: &[Locus], directory: &Path) -> Result<BTreeSet<usize>> {
    let mut out = writer(&directory.join("loci.tsv"))?;
    writeln!(out, "gene\tsample\tcontig\tlocus\tstart\tend\tstrand\tmodels\tupstream\tdownstream\tbundle_evidence\talleles\tannotation_ids\torigins\tspan_length_disagreements")?;
    let mut by_sample: BTreeMap<_, Vec<&Locus>> = BTreeMap::new();
    let mut uncertain = BTreeSet::new();
    for locus in loci {
        let family = &input.families[locus.family];
        let name = |index: Option<usize>| {
            index
                .map(|index| input.families[index].reference.id.as_str())
                .unwrap_or("")
        };
        let alleles: BTreeSet<_> = locus
            .members
            .iter()
            .map(|&index| input.alleles[input.members[index].allele].id.as_str())
            .collect();
        let annotations: Vec<_> = locus
            .members
            .iter()
            .map(|&index| input.members[index].annotation.as_str())
            .collect();
        let origins: BTreeSet<_> = locus
            .members
            .iter()
            .map(|&index| input.members[index].origin)
            .collect();
        let span_disagreements = locus
            .members
            .iter()
            .filter(|&&index| {
                let member = &input.members[index];
                u64::from(member.end) - u64::from(member.start) + 1
                    != input.alleles[member.allele].dna.len() as u64
            })
            .count();
        let evidence = bundle_evidence(locus, input);
        if evidence.ends_with("requires_review") || span_disagreements > 0 {
            uncertain.insert(locus.family);
        }
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            family.reference.id,
            input.samples[locus.sample],
            input.contigs[locus.contig],
            locus.id,
            locus.start,
            locus.end,
            locus.strand,
            locus.members.len(),
            name(locus.flanks.0),
            name(locus.flanks.1),
            evidence,
            list(alleles),
            list(annotations),
            list(origins),
            span_disagreements
        )?;
        by_sample
            .entry((locus.family, locus.sample))
            .or_default()
            .push(locus);
    }
    out.flush()?;
    let mut copies = writer(&directory.join("copies.tsv"))?;
    writeln!(copies, "gene\tsample\tphysical_loci\tintact_near_full_length_loci\tspatially_separated_pairs\tevidence\tlocus_ids")?;
    for ((family, sample), loci) in by_sample {
        if loci.len() < 2 {
            continue;
        }
        uncertain.insert(family);
        let reference_len = input.families[family].reference.protein.len();
        let full = loci
            .iter()
            .filter(|locus| {
                locus.members.iter().any(|&index| {
                    let allele = &input.alleles[input.members[index].allele];
                    allele.orf.is_intact() && allele.protein.len() * 100 >= reference_len * 95
                })
            })
            .count();
        let mut separated = 0;
        for (position, a) in loci.iter().enumerate() {
            for b in &loci[position + 1..] {
                separated +=
                    usize::from(a.contig != b.contig || a.end < b.start || b.end < a.start);
            }
        }
        writeln!(
            copies,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            input.families[family].reference.id,
            input.samples[sample],
            loci.len(),
            full,
            separated,
            if full >= 2 {
                "cooccurring_intact_near_full_length_models"
            } else {
                "multiple_loci_with_partial_or_damaged_models"
            },
            loci.iter()
                .map(|locus| format!("{}:{}", escaped(&input.contigs[locus.contig]), locus.id))
                .collect::<Vec<_>>()
                .join(";")
        )?;
    }
    copies.flush()?;
    Ok(uncertain)
}

#[derive(Clone)]
struct Assessment {
    allele: usize,
    profile: Option<Profile>,
    status: &'static str,
    paired: usize,
    matches: usize,
    internal_loss: usize,
    internal_insert: usize,
    terminal_loss: usize,
    terminal_insert: usize,
    blocks: Vec<String>,
}

fn assess(reference: &Sequence, allele_index: usize, allele: &Sequence, cap: usize) -> Assessment {
    let mut result = Assessment {
        allele: allele_index,
        profile: None,
        status: "aligned",
        paired: 0,
        matches: 0,
        internal_loss: 0,
        internal_insert: 0,
        terminal_loss: 0,
        terminal_insert: 0,
        blocks: Vec::new(),
    };
    let alignment = match align(&reference.protein, &allele.protein, cap) {
        Ok(alignment) => alignment,
        Err(status) => {
            result.status = status;
            return result;
        }
    };
    let mut r = 0;
    let mut a = 0;
    let mut cursor = 0;
    while cursor < alignment.steps.len() {
        let step = alignment.steps[cursor];
        if step == Step::Both {
            let x = reference.protein[r];
            let y = allele.protein[a];
            if !matches!(x, b'X' | b'*') && !matches!(y, b'X' | b'*') {
                result.paired += 1;
                result.matches += usize::from(x == y);
            }
            r += 1;
            a += 1;
            cursor += 1;
            continue;
        }
        let start = cursor;
        while cursor < alignment.steps.len() && alignment.steps[cursor] == step {
            cursor += 1;
        }
        let length = cursor - start;
        if step == Step::ReferenceOnly {
            let kind = if r == 0 {
                "reference_n_terminal_gap"
            } else if r + length == reference.protein.len() {
                "reference_c_terminal_gap"
            } else {
                "reference_internal_gap"
            };
            if kind == "reference_internal_gap" {
                result.internal_loss += length;
            } else {
                result.terminal_loss += length;
            }
            result
                .blocks
                .push(format!("{kind}\t{}\t{}\t\t\t\t{length}", r + 1, r + length));
            r += length;
        } else {
            let kind = if r == 0 {
                "allele_n_terminal_extension"
            } else if r == reference.protein.len() {
                "allele_c_terminal_extension"
            } else {
                "allele_internal_insertion"
            };
            if kind == "allele_internal_insertion" {
                result.internal_insert += length;
            } else {
                result.terminal_insert += length;
            }
            result.blocks.push(format!(
                "{kind}\t\t\t{}\t{}\t{r}\t{length}",
                a + 1,
                a + length
            ));
            a += length;
        }
    }
    result.profile = Some(Profile::new(reference.protein.len(), &alignment));
    if result.paired == 0 {
        result.status = "aligned_no_canonical_paired_residues";
    }
    result
}

struct FamilyAssessment {
    alleles: Vec<Assessment>,
    uncertain: bool,
}

fn protein_msa(reference: &[u8], rows: &[(&[u8], &Profile)]) -> Vec<Vec<u8>> {
    let widths: Vec<_> = (0..=reference.len())
        .map(|position| {
            rows.iter()
                .map(|(_, profile)| profile.insertions[position].len())
                .max()
                .unwrap_or(0)
        })
        .collect();
    let mut result = Vec::new();
    let mut reference_row = Vec::new();
    for position in 0..=reference.len() {
        reference_row.extend(std::iter::repeat(b'-').take(widths[position]));
        if position < reference.len() {
            reference_row.push(reference[position]);
        }
    }
    result.push(reference_row);
    for &(protein, profile) in rows {
        let mut row = Vec::new();
        for position in 0..=reference.len() {
            row.extend(
                profile.insertions[position]
                    .iter()
                    .map(|&index| protein[index]),
            );
            row.extend(
                std::iter::repeat(b'-').take(widths[position] - profile.insertions[position].len()),
            );
            if position < reference.len() {
                row.push(
                    profile.positions[position]
                        .map(|index| protein[index])
                        .unwrap_or(b'-'),
                );
            }
        }
        result.push(row);
    }
    result
}

fn write_alignment(path: &Path, labels: &[String], rows: &[Vec<u8>]) -> Result<()> {
    if labels.len() != rows.len()
        || rows
            .first()
            .is_some_and(|first| rows.iter().any(|row| row.len() != first.len()))
        || labels.iter().collect::<HashSet<_>>().len() != labels.len()
    {
        return Err(Error::new("invalid or nonunique validation alignment rows"));
    }
    let mut out = writer(path)?;
    for (label, row) in labels.iter().zip(rows) {
        writeln!(out, ">{label}")?;
        for chunk in row.chunks(80) {
            out.write_all(chunk)?;
            writeln!(out)?;
        }
    }
    out.flush()?;
    Ok(())
}

fn family_assessment(family: &Family, input: &Input, cap: usize) -> FamilyAssessment {
    let indexes: Vec<usize> = family.alleles.keys().copied().collect();
    let mut first_of: HashMap<&[u8], usize> = HashMap::new();
    let mut unique = Vec::new();
    let mut slot_of = Vec::with_capacity(indexes.len());
    for &index in &indexes {
        let protein = input.alleles[index].protein.as_slice();
        if let Some(&slot) = first_of.get(protein) {
            slot_of.push(slot);
            continue;
        }
        let slot = unique.len();
        first_of.insert(protein, slot);
        unique.push(index);
        slot_of.push(slot);
    }
    let scored: Vec<Assessment> = if unique.len() >= 2 {
        unique
            .par_iter()
            .map(|&index| assess(&family.reference, index, &input.alleles[index], cap))
            .collect()
    } else {
        unique
            .iter()
            .map(|&index| assess(&family.reference, index, &input.alleles[index], cap))
            .collect()
    };
    let alleles: Vec<_> = indexes
        .into_iter()
        .zip(slot_of)
        .map(|(index, slot)| {
            let mut report = scored[slot].clone();
            report.allele = index;
            report
        })
        .collect();
    let uncertain = !family.reference.orf.is_intact()
        || alleles.iter().any(|row| {
            row.status != "aligned"
                || row.internal_loss > 0
                || row.internal_insert > 0
                || row.terminal_loss > 0
                || row.terminal_insert > 0
                || row
                    .profile
                    .as_ref()
                    .is_some_and(|profile| profile.ambiguous)
                || !input.alleles[row.allele].orf.is_intact()
                || input.alleles[row.allele].protein.contains(&b'X')
                || row.matches * 100 < row.paired * 70
                || row.paired * 100
                    < family
                        .reference
                        .protein
                        .len()
                        .min(input.alleles[row.allele].protein.len())
                        * 80
        });
    FamilyAssessment { alleles, uncertain }
}

fn write_assessments(
    input: &Input,
    reports: &[FamilyAssessment],
    directory: &Path,
) -> Result<usize> {
    let mut out = writer(&directory.join("alleles.tsv"))?;
    let mut blocks = writer(&directory.join("structural_blocks.tsv"))?;
    let mut families = writer(&directory.join("families.tsv"))?;
    let mut aliases = writer(&directory.join("allele_aliases.tsv"))?;
    writeln!(aliases, "input_allele\tcanonical_allele")?;
    for (alias, &index) in &input.aliases {
        writeln!(aliases, "{alias}\t{}", input.alleles[index].id)?;
    }
    aliases.flush()?;
    writeln!(out, "gene\tallele\toccurrences\tstatus\tprotein_aa\treference_aa\tpaired_canonical_aa\tidentity_paired_pct\treference_coverage_pct\tallele_coverage_pct\tinternal_loss_aa\tinternal_insert_aa\tterminal_loss_aa\tterminal_insert_aa\talignment_ambiguous\tin_frame\thas_start\tinternal_stops\tterminal_stop\tunknown_aa\tnamed_domains")?;
    writeln!(blocks, "gene\tallele\tfeature\treference_start_aa\treference_end_aa\tallele_start_aa\tallele_end_aa\treference_boundary_aa\tlength_aa")?;
    writeln!(families, "family_file_id\tgene\toccurrences\tunique_dna_alleles\taligned_alleles\tsequence_review_required\tnamed_domains")?;
    let mut unassessed = 0;
    for (index, (family, report)) in input.families.iter().zip(reports).enumerate() {
        let file_id = format!("f{:06}", index + 1);
        let mut rows = Vec::new();
        let mut reference_label = "reference".to_string();
        while input.aliases.contains_key(&reference_label) {
            reference_label.push('_');
        }
        let mut labels = vec![reference_label];
        for row in &report.alleles {
            let allele = &input.alleles[row.allele];
            let ratio = |numerator: usize, denominator: usize| {
                if denominator == 0 {
                    String::new()
                } else {
                    format!("{:.6}", numerator as f64 * 100.0 / denominator as f64)
                }
            };
            let identity = if row.profile.is_some() {
                ratio(row.matches, row.paired)
            } else {
                String::new()
            };
            let ref_coverage = if row.profile.is_some() {
                ratio(row.paired, family.reference.protein.len())
            } else {
                String::new()
            };
            let allele_coverage = if row.profile.is_some() {
                ratio(row.paired, allele.protein.len())
            } else {
                String::new()
            };
            let measured = |value: usize| {
                if row.profile.is_some() {
                    value.to_string()
                } else {
                    String::new()
                }
            };
            writeln!(out, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tunassessed",
                family.reference.id, allele.id, family.alleles[&row.allele], row.status,
                allele.protein.len(), family.reference.protein.len(), measured(row.paired), identity, ref_coverage,
                allele_coverage, measured(row.internal_loss), measured(row.internal_insert), measured(row.terminal_loss), measured(row.terminal_insert),
                row.profile.as_ref().map(|profile| profile.ambiguous.to_string()).unwrap_or_default(),
                allele.orf.in_frame, allele.orf.has_start, allele.orf.internal_stops, allele.orf.terminal_stop,
                allele.protein.iter().filter(|&&aa| aa == b'X').count())?;
            for block in &row.blocks {
                writeln!(blocks, "{}\t{}\t{block}", family.reference.id, allele.id)?;
            }
            if let Some(profile) = &row.profile {
                rows.push((allele.protein.as_slice(), profile));
                labels.push(allele.id.clone());
            } else {
                unassessed += 1;
            }
        }
        if !family.reference.protein.is_empty() {
            let msa = protein_msa(&family.reference.protein, &rows);
            write_alignment(
                &directory.join("alignments").join(format!("{file_id}.faa")),
                &labels,
                &msa,
            )?;
            let mut codons = vec![family.reference.dna.as_bytes()
                [..family.reference.protein.len() * 3]
                .to_vec()];
            for row in &report.alleles {
                if let Some(profile) = &row.profile {
                    codons.push(profile.codon_core(input.alleles[row.allele].dna.as_bytes()));
                }
            }
            write_alignment(
                &directory
                    .join("alignments")
                    .join(format!("{file_id}.core.fna")),
                &labels,
                &codons,
            )?;
        }
        writeln!(
            families,
            "{file_id}\t{}\t{}\t{}\t{}\t{}\tunassessed",
            family.reference.id,
            family.members.len(),
            family.alleles.len(),
            rows.len(),
            report.uncertain
        )?;
    }
    out.flush()?;
    blocks.flush()?;
    families.flush()?;
    Ok(unassessed)
}

fn related_sets(
    input: &Input,
    loci: &[Locus],
    directory: &Path,
    cap: usize,
) -> Result<(Vec<Vec<usize>>, usize)> {
    let mut contexts: BTreeMap<_, BTreeSet<usize>> = BTreeMap::new();
    for locus in loci {
        for (side, neighbor) in [(0u8, locus.flanks.0), (1u8, locus.flanks.1)] {
            if let Some(neighbor) = neighbor {
                contexts
                    .entry((side, neighbor))
                    .or_default()
                    .insert(locus.family);
            }
        }
    }
    let mut candidates = BTreeSet::new();
    for groups in contexts.values() {
        let groups: Vec<_> = groups.iter().copied().collect();
        for (position, &a) in groups.iter().enumerate() {
            for &b in &groups[position + 1..] {
                candidates.insert((a, b));
            }
        }
    }
    let kmers: Vec<HashSet<[u8; 4]>> = input
        .families
        .par_iter()
        .map(|family| {
            family
                .reference
                .protein
                .windows(4)
                .filter(|word| !word.contains(&b'X') && !word.contains(&b'*'))
                .map(|word| [word[0], word[1], word[2], word[3]])
                .collect()
        })
        .collect();
    let mut parents: Vec<_> = (0..input.families.len()).collect();
    fn root(parents: &mut [usize], mut index: usize) -> usize {
        while parents[index] != index {
            parents[index] = parents[parents[index]];
            index = parents[index];
        }
        index
    }
    let candidates: Vec<_> = candidates.into_iter().collect();
    let evidence: Vec<_> = candidates
        .par_iter()
        .map(|&(a, b)| {
            let x = &input.families[a].reference;
            let y = &input.families[b].reference;
            if length_ratio(&x.protein, &y.protein) < 0.5 {
                return None;
            }
            let (small, large) = if kmers[a].len() <= kmers[b].len() {
                (&kmers[a], &kmers[b])
            } else {
                (&kmers[b], &kmers[a])
            };
            if small.is_empty()
                || small.iter().filter(|word| large.contains(*word)).count() * 10 < small.len()
            {
                return None;
            }
            let assessment = assess(x, 0, y, cap);
            if assessment.profile.is_none() {
                return Some((a, b, assessment.status, false, None, None));
            }
            if assessment.paired == 0 {
                return Some((a, b, "no_canonical_paired_residues", false, None, None));
            }
            let shorter = x.protein.len().min(y.protein.len());
            let identity = assessment.matches as f64 / assessment.paired as f64;
            let coverage = assessment.paired as f64 / shorter as f64;
            let related = identity >= 0.70 && coverage >= 0.80;
            Some((
                a,
                b,
                if related {
                    "candidate_homologues_not_an_orthology_call"
                } else {
                    "sequence_filter_rejected"
                },
                related,
                Some(identity * 100.0),
                Some(coverage * 100.0),
            ))
        })
        .collect();
    let mut out = writer(&directory.join("related_families.tsv"))?;
    writeln!(
        out,
        "gene_a\tgene_b\tstatus\tpaired_identity_pct\tshorter_coverage_pct"
    )?;
    let mut unassessed = 0;
    for (a, b, status, related, identity, coverage) in evidence.into_iter().flatten() {
        let format_metric =
            |value: Option<f64>| value.map(|value| format!("{value:.6}")).unwrap_or_default();
        if identity.is_none() {
            unassessed += 1;
        }
        writeln!(
            out,
            "{}\t{}\t{status}\t{}\t{}",
            input.families[a].reference.id,
            input.families[b].reference.id,
            format_metric(identity),
            format_metric(coverage)
        )?;
        if related {
            let r = root(&mut parents, a);
            let s = root(&mut parents, b);
            parents[s.max(r)] = r.min(s);
        }
    }
    out.flush()?;
    let mut groups: BTreeMap<_, Vec<usize>> = BTreeMap::new();
    for index in 0..parents.len() {
        groups
            .entry(root(&mut parents, index))
            .or_default()
            .push(index);
    }
    Ok((groups.into_values().collect(), unassessed))
}

struct TreeReport {
    status: String,
    anchor: usize,
    labels: Vec<String>,
    rows: Vec<Vec<u8>>,
    sites: usize,
    variable: usize,
    usable_bootstraps: usize,
    rejected_bootstraps: usize,
    negative_branches: usize,
    newick: Option<String>,
    branches: Vec<(Vec<usize>, f64, usize)>,
}

fn tree_report(
    group: &[usize],
    input: &Input,
    assessments: &[FamilyAssessment],
    options: &ValidationOptions,
) -> Result<TreeReport> {
    let anchor = *group
        .iter()
        .max_by_key(|&&index| {
            (
                input.families[index].reference.protein.len(),
                std::cmp::Reverse(index),
            )
        })
        .unwrap();
    let mut report = TreeReport {
        status: "unassessed".to_string(),
        anchor,
        labels: Vec::new(),
        rows: Vec::new(),
        sites: 0,
        variable: 0,
        usable_bootstraps: 0,
        rejected_bootstraps: 0,
        negative_branches: 0,
        newick: None,
        branches: Vec::new(),
    };
    let allele_ids: BTreeSet<_> = group
        .iter()
        .flat_map(|&index| input.families[index].alleles.keys().copied())
        .collect();
    report.labels = allele_ids
        .iter()
        .map(|&index| input.alleles[index].id.clone())
        .collect();
    if allele_ids.len() < 3 {
        report.status = "fewer_than_three_distinct_dna_alleles".into();
        return Ok(report);
    }
    if options.max_tree_alleles != 0 && allele_ids.len() > options.max_tree_alleles {
        report.status = "tree_allele_limit".into();
        return Ok(report);
    }
    if !input.families[anchor].reference.orf.in_frame {
        report.status = "reference_out_of_frame".into();
        return Ok(report);
    }
    let ids: Vec<_> = allele_ids.iter().copied().collect();
    // An out-of-frame member stops the tree at that allele. Aligning a later
    // member first would be able to report a different failure.
    let prefix_end = ids
        .iter()
        .position(|&index| !input.alleles[index].orf.in_frame)
        .unwrap_or(ids.len());
    let mut cache: HashMap<&[u8], Assessment> = assessments[anchor]
        .alleles
        .iter()
        .map(|row| (input.alleles[row.allele].protein.as_slice(), row.clone()))
        .collect();
    let mut missing = Vec::new();
    let mut seen = HashSet::new();
    for &index in &ids[..prefix_end] {
        let protein = input.alleles[index].protein.as_slice();
        if cache.contains_key(protein) || !seen.insert(protein) {
            continue;
        }
        missing.push(index);
    }
    let computed: Vec<Assessment> = if missing.len() >= 2 {
        missing
            .par_iter()
            .map(|&index| {
                assess(
                    &input.families[anchor].reference,
                    index,
                    &input.alleles[index],
                    options.max_alignment_cells,
                )
            })
            .collect()
    } else {
        missing
            .iter()
            .map(|&index| {
                assess(
                    &input.families[anchor].reference,
                    index,
                    &input.alleles[index],
                    options.max_alignment_cells,
                )
            })
            .collect()
    };
    for (index, row) in missing.into_iter().zip(computed) {
        cache.insert(input.alleles[index].protein.as_slice(), row);
    }
    for &index in &ids[..prefix_end] {
        let allele = &input.alleles[index];
        let row = &cache[allele.protein.as_slice()];
        let Some(profile) = &row.profile else {
            report.status = format!("tree_{}", row.status);
            return Ok(report);
        };
        if profile.ambiguous {
            report.status = "ambiguous_alignment".into();
            return Ok(report);
        }
        report.rows.push(profile.codon_core(allele.dna.as_bytes()));
    }
    if prefix_end < ids.len() {
        report.status = "member_out_of_frame".into();
        return Ok(report);
    }
    let sites = common_sites(&report.rows)?;
    report.sites = sites.len();
    report.variable = variable_sites(&report.rows, &sites);
    if sites.len() < options.min_tree_sites {
        report.status = "insufficient_common_sites".into();
        return Ok(report);
    }
    if report.variable == 0 {
        report.status = "no_variable_common_sites".into();
        return Ok(report);
    }
    let matrix = match distances(&report.rows, &sites) {
        Ok(matrix) => matrix,
        Err(error) => {
            report.status = error.to_string();
            return Ok(report);
        }
    };
    let tree = neighbor_joining(&matrix)?;
    let (support, usable, rejected) =
        bootstrap(&report.rows, &sites, options.bootstrap_replicates)?;
    report.usable_bootstraps = usable;
    report.rejected_bootstraps = rejected;
    report.negative_branches = tree.negative_branches();
    report.branches = tree
        .splits()
        .into_iter()
        .map(|(split, length)| {
            let count = support.get(&split).copied().unwrap_or(0);
            (split, length, count)
        })
        .collect();
    report.newick = Some(tree.newick(
        &report.labels,
        &support,
        if usable == 0 {
            0
        } else {
            options.bootstrap_replicates
        },
    ));
    report.status = if report.negative_branches > 0 {
        "diagnostic_tree_negative_branches"
    } else if options.bootstrap_replicates > 0 && usable == 0 {
        "diagnostic_tree_no_usable_bootstraps"
    } else {
        "diagnostic_unrooted_tree_not_an_orthology_call"
    }
    .into();
    Ok(report)
}

fn write_trees(
    groups: &[Vec<usize>],
    input: &Input,
    reports: &[FamilyAssessment],
    uncertain: &BTreeSet<usize>,
    directory: &Path,
    options: &ValidationOptions,
) -> Result<usize> {
    let selected: Vec<_> = groups
        .iter()
        .filter(|group| match options.tree_mode {
            TreeMode::None => false,
            TreeMode::All => true,
            TreeMode::Uncertain => {
                group.len() > 1
                    || group
                        .iter()
                        .any(|index| uncertain.contains(index) || reports[*index].uncertain)
            }
        })
        .collect();
    let results: Vec<_> = selected
        .par_iter()
        .map(|group| tree_report(group, input, reports, options))
        .collect::<Result<_>>()?;
    let requested = selected.len();
    let mut out = writer(&directory.join("trees.tsv"))?;
    let mut branches = writer(&directory.join("tree_branches.tsv"))?;
    let mut members = writer(&directory.join("tree_members.tsv"))?;
    writeln!(out, "tree_id\tgenes\tanchor_reference\tstatus\tunique_alleles\tcommon_sites\tvariable_sites\trequested_bootstraps\tusable_bootstraps\trejected_bootstraps\tnegative_branches")?;
    writeln!(
        branches,
        "tree_id\tbipartition_alleles\tlength\tbootstrap_support_pct\tusable_bootstraps"
    )?;
    writeln!(
        members,
        "tree_id\tallele\tgene\tsample\tcontig\tlocus\tannotation_id\torigin"
    )?;
    let mut written = 0;
    for (index, (group, report)) in selected.into_iter().zip(results).enumerate() {
        let id = format!("t{:06}", index + 1);
        let names: Vec<_> = group
            .iter()
            .map(|&index| input.families[index].reference.id.as_str())
            .collect();
        writeln!(
            out,
            "{id}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            list(names),
            input.families[report.anchor].reference.id,
            report.status,
            report.labels.len(),
            report.sites,
            report.variable,
            options.bootstrap_replicates,
            report.usable_bootstraps,
            report.rejected_bootstraps,
            report.negative_branches
        )?;
        for &family in group {
            for &member in &input.families[family].members {
                let member = &input.members[member];
                writeln!(
                    members,
                    "{id}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    input.alleles[member.allele].id,
                    input.families[family].reference.id,
                    input.samples[member.sample],
                    input.contigs[member.contig],
                    member.locus,
                    member.annotation,
                    member.origin
                )?;
            }
        }
        if let Some(newick) = report.newick {
            let mut file = writer(&directory.join("trees").join(format!("{id}.nwk")))?;
            writeln!(file, "{newick}")?;
            file.flush()?;
            write_alignment(
                &directory.join("trees").join(format!("{id}.core.fna")),
                &report.labels,
                &report.rows,
            )?;
            for (split, length, support) in report.branches {
                let names: Vec<_> = split
                    .iter()
                    .map(|&index| report.labels[index].as_str())
                    .collect();
                let pct = if report.usable_bootstraps == 0 {
                    String::new()
                } else {
                    format!(
                        "{:.6}",
                        support as f64 * 100.0 / options.bootstrap_replicates as f64
                    )
                };
                writeln!(
                    branches,
                    "{id}\t{}\t{length:.8}\t{pct}\t{}",
                    list(names),
                    report.usable_bootstraps
                )?;
            }
            written += 1;
        }
    }
    out.flush()?;
    branches.flush()?;
    members.flush()?;
    progress!("rustroo validate: {written}/{requested} candidate sets have diagnostic trees; {} have explicit no-tree reasons in trees.tsv",
        requested - written);
    Ok(written)
}

fn fingerprint(path: &Path) -> Result<u64> {
    let mut input = BufReader::new(File::open(path)?);
    let mut buffer = [0u8; 65536];
    let mut hash = 0xcbf29ce484222325u64;
    loop {
        let count = input.read(&mut buffer)?;
        if count == 0 {
            break;
        }
        for &byte in &buffer[..count] {
            hash = (hash ^ u64::from(byte)).wrapping_mul(0x100000001b3);
        }
    }
    Ok(hash)
}

/// Read-only validation of a completed membership-enabled pangenome.
pub fn validate(
    directory: &Path,
    output: &Path,
    options: &ValidationOptions,
) -> Result<ValidationSummary> {
    if options.min_tree_sites == 0
        || (options.max_tree_alleles != 0 && options.max_tree_alleles < 3)
    {
        return Err(Error::new(
            "min-tree-sites must be positive and max-tree-alleles must be 0 or at least 3",
        ));
    }
    if output.exists() {
        return Err(Error::new(
            "validation output already exists; use a fresh directory",
        ));
    }
    let files = [
        "genes.fa",
        "presence.tsv",
        "members.tsv",
        "member_alleles.fa",
    ];
    let mut fingerprints = Vec::new();
    for name in files {
        let path = directory.join(name);
        if !path.is_file() {
            return Err(Error::new(format!(
                "missing validation input {}; build with --audit-members",
                path.display()
            )));
        }
        fingerprints.push((name, fingerprint(&path)?));
    }
    let input = load(directory)?;
    let loci = loci(&input)?;
    if let Some(parent) = output.parent().filter(|path| !path.as_os_str().is_empty()) {
        fs::create_dir_all(parent)?;
    }
    fs::create_dir(output)?;
    fs::create_dir(output.join("alignments"))?;
    fs::create_dir(output.join("trees"))?;
    log::open(&output.join("rustroo.log"))?;
    let _log_guard = LogGuard;
    progress!("rustroo validate: reading pooled alleles and sample-scoped memberships");
    let uncertain = write_loci(&input, &loci, output)?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(options.threads)
        .build()?;
    progress!("rustroo validate: aligning every distinct member DNA allele to its family reference on {} threads", pool.current_num_threads());
    let done_count = AtomicUsize::new(0);
    let reports: Vec<_> = pool.install(|| {
        input
            .families
            .par_iter()
            .map(|family| {
                let report = family_assessment(family, &input, options.max_alignment_cells);
                let done = done_count.fetch_add(1, Ordering::Relaxed) + 1;
                if done % 250 == 0 || done == input.families.len() {
                    progress!(
                        "rustroo validate: aligned {done}/{} families",
                        input.families.len()
                    );
                }
                report
            })
            .collect()
    });
    let unassessed_alignments = write_assessments(&input, &reports, output)?;
    if unassessed_alignments > 0 {
        progress!("rustroo validate: {unassessed_alignments} allele alignments unassessed; inspect alleles.tsv, no heuristic replacements were used");
    }
    progress!("rustroo validate: checking context-indexed candidate homologues");
    let (groups, unassessed_related_pairs) =
        pool.install(|| related_sets(&input, &loci, output, options.max_alignment_cells))?;
    if unassessed_related_pairs > 0 {
        progress!("rustroo validate: {unassessed_related_pairs} candidate-family comparisons unassessed; inspect related_families.tsv");
    }
    progress!("rustroo validate: constructing diagnostic gene trees; named domains and species-tree reconciliation remain unassessed");
    let trees_written =
        pool.install(|| write_trees(&groups, &input, &reports, &uncertain, output, options))?;
    let mut provenance = writer(&output.join("input_fingerprints.tsv"))?;
    writeln!(provenance, "file\tfnv1a64_noncryptographic")?;
    for (name, before) in fingerprints {
        if fingerprint(&directory.join(name))? != before {
            return Err(Error::new(format!(
                "validation input {name} changed during analysis; reports are incomplete"
            )));
        }
        writeln!(provenance, "{name}\t{before:016x}")?;
    }
    provenance.flush()?;
    let mut locus_counts = HashMap::new();
    for locus in &loci {
        *locus_counts
            .entry((locus.family, locus.sample))
            .or_insert(0usize) += 1;
    }
    let copies: BTreeSet<_> = locus_counts
        .into_iter()
        .filter(|(_, count)| *count > 1)
        .map(|((family, _), _)| family)
        .collect();
    let summary = ValidationSummary {
        families: input.families.len(),
        unique_member_alleles: input
            .families
            .iter()
            .flat_map(|family| family.alleles.keys().copied())
            .collect::<BTreeSet<_>>()
            .len(),
        physical_loci: loci.len(),
        families_with_multiple_loci_per_sample: copies.len(),
        unassessed_alignments,
        unassessed_related_pairs,
        trees_written,
    };
    let mut out = writer(&output.join("summary.tsv"))?;
    writeln!(out, "metric\tvalue")?;
    for (name, value) in [
        ("families", summary.families),
        ("unique_member_alleles", summary.unique_member_alleles),
        ("physical_loci", summary.physical_loci),
        (
            "families_with_multiple_loci_per_sample",
            summary.families_with_multiple_loci_per_sample,
        ),
        ("unassessed_alignments", summary.unassessed_alignments),
        ("trees_written", summary.trees_written),
        ("unassessed_related_pairs", summary.unassessed_related_pairs),
        ("threads", pool.current_num_threads()),
        ("bootstrap_replicates", options.bootstrap_replicates),
        ("max_tree_alleles", options.max_tree_alleles),
        ("max_alignment_cells", options.max_alignment_cells),
        ("min_tree_sites", options.min_tree_sites),
    ] {
        writeln!(out, "{name}\t{value}")?;
    }
    writeln!(out, "validator_version\t{}\ntree_mode\t{:?}\nnamed_domains\tunassessed\nspecies_tree_reconciliation\tunassessed\nassignments_modified\tfalse",
        env!("CARGO_PKG_VERSION"), options.tree_mode)?;
    writeln!(out, "protein_alignment_scoring\tmatch=2;mismatch=-1;gap_open=-5;gap_extend=-1\ntree_distance_model\tJC69\ntree_rooting\tunrooted\ntree_sites\tcomplete_case_reference_codon_columns\nbootstrap_denominator\trequested_replicates\nbootstrap_seed\t72a914bc5e378d01")?;
    writeln!(out, "complete\ttrue")?;
    out.flush()?;
    progress!("rustroo validate: complete; {trees_written} diagnostic trees; original assignments and references unchanged");
    Ok(summary)
}

struct LogGuard;

impl Drop for LogGuard {
    fn drop(&mut self) {
        log::close();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genetic_code::back_translate;
    use crate::{write_pangenome, Gene, GeneMember, MembershipAudit, Pangenome};
    use std::path::PathBuf;

    static NEXT: AtomicUsize = AtomicUsize::new(0);

    struct Directory(PathBuf);

    impl Directory {
        fn new() -> Self {
            let path = std::env::temp_dir().join(format!(
                "rustroo-validation-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
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

    fn fixture(path: &Path, truncated: bool) {
        let protein = b"MACDEFGHIKLMNPQRSTVWY";
        let full = back_translate(protein);
        let mut synonymous = full.clone();
        synonymous.replace_range(3..6, "GCC");
        let third = if truncated {
            back_translate(&protein[5..])
        } else {
            let mut third = full.clone();
            third.replace_range(3..6, "GCA");
            third
        };
        let sequences = vec![full.clone(), synonymous, third];
        let members = (0..3)
            .map(|index| {
                let sample = u32::from(index == 2);
                let contig = u32::from(index == 1);
                GeneMember {
                    gene: 0,
                    sample,
                    contig,
                    annotation_id: "reused;annotation:id".into(),
                    start: 1,
                    end: sequences[index].len() as u32,
                    strand: if index == 1 { -1 } else { 1 },
                    locus: 0,
                    dna: index as u32,
                    origin: "annotated",
                }
            })
            .collect();
        let pangenome = Pangenome {
            sample_names: vec!["alpha".into(), "beta".into()],
            genes: vec![Gene {
                id: "gene".into(),
                product: String::new(),
                reference_dna: full,
                intact: true,
                selection: "intact",
                presence: vec![true, true],
                lengths: vec![
                    (sequences[0].len() + sequences[1].len()) as u32,
                    sequences[2].len() as u32,
                ],
                samples_present: 2,
                length_evidence: "conserved_short",
                samples_with_exact_reference: 1,
                distinct_alleles: 3,
                intact_alleles: 3,
                modal_end_bp: sequences[0].len() as u32,
                modal_end_samples: 2,
                extension_samples: 0,
                support: "core",
            }],
            membership_audit: Some(MembershipAudit {
                contig_names: vec![
                    vec!["contig:1".into(), "contig;2".into()],
                    vec!["contig:1".into()],
                ],
                members,
                sequences,
            }),
        };
        write_pangenome(&pangenome, path).unwrap();
    }

    fn options() -> ValidationOptions {
        ValidationOptions {
            threads: 1,
            bootstrap_replicates: 10,
            min_tree_sites: 3,
            ..ValidationOptions::default()
        }
    }

    #[test]
    fn validation_covers_members_and_copies_without_modifying_inputs() {
        let directory = Directory::new();
        let input = directory.0.join("input");
        let output = directory.0.join("validation");
        fixture(&input, false);
        let names = [
            "genes.fa",
            "presence.tsv",
            "members.tsv",
            "member_alleles.fa",
        ];
        let before: Vec<_> = names
            .iter()
            .map(|name| fs::read(input.join(name)).unwrap())
            .collect();
        let summary = validate(&input, &output, &options()).unwrap();
        assert_eq!(summary.families, 1);
        assert_eq!(summary.unique_member_alleles, 3);
        assert_eq!(summary.physical_loci, 3);
        assert_eq!(summary.families_with_multiple_loci_per_sample, 1);
        assert_eq!(summary.unassessed_alignments, 0);
        assert_eq!(summary.trees_written, 1);
        assert!(fs::read_to_string(output.join("copies.tsv"))
            .unwrap()
            .contains("cooccurring_intact_near_full_length_models"));
        assert!(fs::read_to_string(output.join("loci.tsv"))
            .unwrap()
            .contains("reused%3Bannotation%3Aid"));
        assert_eq!(
            fs::read_to_string(output.join("alleles.tsv"))
                .unwrap()
                .lines()
                .count(),
            4
        );
        assert!(fs::read_to_string(output.join("alleles.tsv"))
            .unwrap()
            .lines()
            .skip(1)
            .all(|row| row.ends_with("\tunassessed")));
        assert_eq!(
            fs::read_to_string(output.join("tree_members.tsv"))
                .unwrap()
                .lines()
                .count(),
            4
        );
        assert!(fs::read_to_string(output.join("summary.tsv"))
            .unwrap()
            .contains("assignments_modified\tfalse"));
        for (name, bytes) in names.iter().zip(before) {
            assert_eq!(fs::read(input.join(name)).unwrap(), bytes);
        }
        assert!(validate(&input, &output, &options())
            .unwrap_err()
            .to_string()
            .contains("already exists"));
    }

    #[test]
    fn unavailable_alignment_metrics_are_blank_and_trees_are_not_fabricated() {
        let directory = Directory::new();
        let input = directory.0.join("input");
        let output = directory.0.join("validation");
        fixture(&input, true);
        let mut options = options();
        options.max_alignment_cells = 1;
        let summary = validate(&input, &output, &options).unwrap();
        assert_eq!(summary.unassessed_alignments, 1);
        assert_eq!(summary.trees_written, 0);
        let text = fs::read_to_string(output.join("alleles.tsv")).unwrap();
        let fields: Vec<_> = text
            .lines()
            .find(|line| line.contains("alignment_work_limit"))
            .unwrap()
            .split('\t')
            .collect();
        assert_eq!(fields[6], "");
        assert_eq!(fields[7], "");
        assert!(fs::read_to_string(output.join("trees.tsv"))
            .unwrap()
            .contains("tree_alignment_work_limit"));
        assert!(!output.join("trees").join("t000001.nwk").exists());
    }

    #[test]
    fn protein_alignment_retains_all_inserted_residues_but_tree_core_excludes_them() {
        let reference = b"MACDEFGHIKLMNPQRSTVWY";
        let mut allele = reference[..8].to_vec();
        allele.extend_from_slice(b"AAAAA");
        allele.extend_from_slice(&reference[8..]);
        let profile = Profile::new(reference.len(), &align(reference, &allele, 0).unwrap());
        let msa = protein_msa(reference, &[(&allele, &profile)]);
        assert_eq!(msa[0].len(), msa[1].len());
        assert_eq!(
            msa[1]
                .iter()
                .filter(|&&base| base != b'-')
                .copied()
                .collect::<Vec<_>>(),
            allele
        );
        assert_eq!(
            profile.codon_core(back_translate(&allele).as_bytes()).len(),
            reference.len() * 3
        );
    }

    #[test]
    fn exact_dna_aliases_are_pooled_and_preserved_in_a_mapping() {
        let directory = Directory::new();
        let input = directory.0.join("input");
        fixture(&input, false);
        let full = back_translate(b"MACDEFGHIKLMNPQRSTVWY");
        let mut file = OpenOptions::new()
            .append(true)
            .open(input.join("member_alleles.fa"))
            .unwrap();
        writeln!(file, ">alias\n{full}").unwrap();
        drop(file);
        let loaded = load(&input).unwrap();
        assert_eq!(loaded.alleles.len(), 3);
        assert_eq!(loaded.aliases["alias"], loaded.aliases["dna_0"]);
    }

    #[test]
    fn malformed_provenance_and_presence_fail_before_creating_reports() {
        let directory = Directory::new();
        let input = directory.0.join("input");
        let output = directory.0.join("validation");
        fixture(&input, false);
        fs::write(
            input.join("presence.tsv"),
            "gene\talpha\tbeta\ngene\t0\t1\n",
        )
        .unwrap();
        assert!(validate(&input, &output, &options())
            .unwrap_err()
            .to_string()
            .contains("disagree"));
        assert!(!output.exists());
        fs::write(input.join("bad.fa"), ">id\nATG?\n").unwrap();
        assert!(fasta(&input.join("bad.fa")).is_err());
        fs::write(input.join("bad.fa"), ">id\nATG\n>id\nATG\n").unwrap();
        assert!(fasta(&input.join("bad.fa")).is_err());
    }

    #[test]
    fn named_domains_are_not_invented_and_internal_gap_coordinates_are_explicit() {
        let protein = b"MACDEFGHIKLMNPQRSTVWY";
        let reference = Sequence::new("reference".into(), back_translate(protein));
        let mut deleted = protein[..7].to_vec();
        deleted.extend_from_slice(&protein[12..]);
        let allele = Sequence::new("allele".into(), back_translate(&deleted));
        let report = assess(&reference, 0, &allele, 0);
        assert_eq!(report.internal_loss, 5);
        assert!(report
            .blocks
            .iter()
            .any(|block| block.starts_with("reference_internal_gap\t8\t12")));
    }

    #[test]
    fn candidate_homology_sets_never_change_family_assignments() {
        let directory = Directory::new();
        let input_path = directory.0.join("input");
        fixture(&input_path, false);
        let mut input = load(&input_path).unwrap();
        let dna = input.families[0].reference.dna.clone();
        input.families.push(Family {
            reference: Sequence::new("other_copy".into(), dna),
            members: Vec::new(),
            alleles: BTreeMap::new(),
        });
        input.families.push(Family {
            reference: Sequence::new("unrelated".into(), back_translate(b"MYWVTSRQPNMLKIHGFEDCA")),
            members: Vec::new(),
            alleles: BTreeMap::new(),
        });
        let original: Vec<_> = input.members.iter().map(|member| member.family).collect();
        let loci: Vec<_> = (0..3)
            .map(|family| Locus {
                family,
                sample: 0,
                contig: 0,
                id: family as u32,
                start: 1 + family as u32 * 200,
                end: 66 + family as u32 * 200,
                strand: 1,
                members: Vec::new(),
                flanks: (Some(2), None),
            })
            .collect();
        let (groups, unassessed) = related_sets(&input, &loci, &directory.0, 1).unwrap();
        assert_eq!(groups, vec![vec![0, 1], vec![2]]);
        assert_eq!(unassessed, 0);
        assert_eq!(
            input
                .members
                .iter()
                .map(|member| member.family)
                .collect::<Vec<_>>(),
            original
        );
        assert!(fs::read_to_string(directory.0.join("related_families.tsv"))
            .unwrap()
            .contains("candidate_homologues_not_an_orthology_call"));
    }

    #[test]
    fn reports_are_identical_across_validation_thread_counts() {
        let directory = Directory::new();
        let input = directory.0.join("input");
        fixture(&input, false);
        let serial = directory.0.join("serial");
        let parallel = directory.0.join("parallel");
        let mut options = options();
        validate(&input, &serial, &options).unwrap();
        options.threads = 2;
        validate(&input, &parallel, &options).unwrap();
        for path in [
            "alleles.tsv",
            "loci.tsv",
            "copies.tsv",
            "families.tsv",
            "structural_blocks.tsv",
            "related_families.tsv",
            "trees.tsv",
            "tree_branches.tsv",
            "tree_members.tsv",
            "allele_aliases.tsv",
            "input_fingerprints.tsv",
        ] {
            assert_eq!(
                fs::read(serial.join(path)).unwrap(),
                fs::read(parallel.join(path)).unwrap()
            );
        }
        assert_eq!(
            fs::read(serial.join("trees").join("t000001.nwk")).unwrap(),
            fs::read(parallel.join("trees").join("t000001.nwk")).unwrap()
        );
    }

    #[test]
    fn inconsistent_exported_locus_assignments_are_rejected() {
        let directory = Directory::new();
        let input = directory.0.join("input");
        fixture(&input, false);
        let mut loaded = load(&input).unwrap();
        loaded.members[1].contig = loaded.members[0].contig;
        assert!(loci(&loaded)
            .unwrap_err()
            .to_string()
            .contains("conflicting families or strands"));
    }

    #[test]
    fn tree_resource_and_site_limits_have_exact_no_tree_outcomes() {
        let directory = Directory::new();
        let input = directory.0.join("input");
        fixture(&input, false);
        let loaded = load(&input).unwrap();
        let mut options = options();
        let reports = vec![family_assessment(
            &loaded.families[0],
            &loaded,
            options.max_alignment_cells,
        )];
        options.max_tree_alleles = 2;
        let report = tree_report(&[0], &loaded, &reports, &options).unwrap();
        assert_eq!(report.status, "tree_allele_limit");
        assert!(report.newick.is_none());
        options.max_tree_alleles = 3;
        let expected_sites = loaded.families[0].reference.protein.len() * 3;
        options.min_tree_sites = expected_sites + 1;
        let report = tree_report(&[0], &loaded, &reports, &options).unwrap();
        assert_eq!(report.sites, expected_sites);
        assert_eq!(report.status, "insufficient_common_sites");
        options.min_tree_sites = expected_sites;
        assert!(tree_report(&[0], &loaded, &reports, &options)
            .unwrap()
            .newick
            .is_some());
    }
}
