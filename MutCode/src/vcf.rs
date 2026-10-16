use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::codec::{apply, replacement, validate, Edit};
use crate::fasta::Reference;
use crate::io::{file_digest, identifier, is_dna, reader};
use crate::{fail, Result};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SampleInfo {
    pub sample: String,
    pub vcf_sha256: String,
}

pub struct Job {
    pub path: PathBuf,
    pub info: SampleInfo,
    column: usize,
    columns: usize,
}

pub struct Calls {
    pub genes: Vec<Vec<Edit>>,
    pub records: usize,
    pub applied: usize,
}

pub fn discover(
    paths: &[PathBuf],
    reference: &Reference,
    selected_sample: Option<&str>,
) -> Result<Vec<Job>> {
    let mut jobs = Vec::new();
    let mut names = HashSet::new();
    for path in paths {
        let hash = file_digest(path)?;
        let mut input = reader(path)?;
        let mut line = String::new();
        let mut found = false;
        let mut contigs = HashSet::new();
        let mut format_seen = false;
        loop {
            line.clear();
            if input.read_line(&mut line)? == 0 {
                break;
            }
            let text = line.trim_end_matches(['\r', '\n']);
            if text.starts_with("##fileformat=VCFv4.") {
                format_seen = true;
            } else if let Some(metadata) = text.strip_prefix("##contig=") {
                let (id, length) = contig_metadata(metadata)?;
                let index = reference.index.get(&id).ok_or_else(|| {
                    format!(
                        "VCF contig {id} is missing from reference {}",
                        path.display()
                    )
                })?;
                if !contigs.insert(id.clone()) {
                    return fail(format!("duplicate VCF contig header: {id}"));
                }
                if let Some(length) = length {
                    if length != reference.sequences[*index].dna.len() {
                        return fail(format!(
                            "VCF/reference length mismatch for {id} in {}",
                            path.display()
                        ));
                    }
                }
            } else if text.starts_with("#CHROM\t") {
                let columns: Vec<_> = text.split('\t').collect();
                if !format_seen
                    || columns.len() < 10
                    || columns[..9]
                        != [
                            "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                        ]
                {
                    return fail(format!(
                        "expected VCFv4 header with sample genotypes: {}",
                        path.display()
                    ));
                }
                let mut file_names = HashSet::new();
                for (column, name) in columns.iter().enumerate().skip(9) {
                    identifier(name)?;
                    if !file_names.insert(*name) {
                        return fail(format!("duplicate sample in VCF: {name}"));
                    }
                    if selected_sample.is_some_and(|selected| selected != *name) {
                        continue;
                    }
                    if !names.insert(name.to_string()) {
                        return fail(format!("sample {name} occurs in multiple VCFs; supply one final VCF per sample"));
                    }
                    jobs.push(Job {
                        path: path.clone(),
                        info: SampleInfo {
                            sample: name.to_string(),
                            vcf_sha256: hash.clone(),
                        },
                        column,
                        columns: columns.len(),
                    });
                }
                found = true;
                break;
            } else if !text.starts_with("##") {
                return fail(format!(
                    "unexpected content before VCF header: {}",
                    path.display()
                ));
            }
        }
        if !found {
            return fail(format!("missing VCF sample header: {}", path.display()));
        }
    }
    if jobs.is_empty() {
        return fail("no selected VCF samples");
    }
    jobs.sort_by(|a, b| a.info.sample.cmp(&b.info.sample));
    Ok(jobs)
}

fn contig_metadata(text: &str) -> Result<(String, Option<usize>)> {
    let body = text
        .strip_prefix('<')
        .and_then(|v| v.strip_suffix('>'))
        .ok_or("invalid VCF contig metadata")?;
    let mut fields = HashMap::new();
    let mut quoted = false;
    let mut escaped = false;
    let mut begin = 0;
    for (index, character) in body
        .char_indices()
        .chain(std::iter::once((body.len(), ',')))
    {
        if character == ',' && !quoted {
            let (key, value) = body[begin..index]
                .split_once('=')
                .ok_or("invalid VCF contig field")?;
            if fields.insert(key, value.trim_matches('"')).is_some() {
                return fail("duplicate VCF contig metadata field");
            }
            begin = index + 1;
        } else if character == '"' && !escaped {
            quoted = !quoted;
        }
        escaped = character == '\\' && !escaped;
    }
    if quoted {
        return fail("unclosed quote in VCF contig metadata");
    }
    let id = fields.get("ID").ok_or("VCF contig has no ID")?.to_string();
    identifier(&id)?;
    let length = fields
        .get("length")
        .map(|value| value.parse())
        .transpose()?;
    Ok((id, length))
}

pub fn read_calls(job: &Job, reference: &Reference, haplotype: Option<usize>) -> Result<Calls> {
    let mut genes = vec![Vec::new(); reference.sequences.len()];
    let mut input = reader(&job.path)?;
    let mut line = String::new();
    let mut records = 0;
    let mut applied = 0;
    let mut line_number = 0;
    let mut seen_header = false;
    loop {
        line.clear();
        if input.read_line(&mut line)? == 0 {
            break;
        }
        line_number += 1;
        let text = line.trim_end_matches(['\r', '\n']);
        if text.starts_with("#CHROM\t") {
            if seen_header {
                return fail("duplicate VCF sample header");
            }
            let header: Vec<_> = text.split('\t').collect();
            if header.len() != job.columns || header[job.column] != job.info.sample {
                return fail("VCF header changed since discovery");
            }
            seen_header = true;
            continue;
        }
        if text.starts_with("##") && !seen_header {
            continue;
        }
        if !seen_header || text.starts_with('#') || text.is_empty() {
            return fail(format!(
                "unexpected VCF content at {}:{line_number}",
                job.path.display()
            ));
        }
        let result = record(text, job, reference, haplotype);
        let (gene, edit) = result.map_err(|error| {
            format!(
                "{}:{line_number}, sample {}: {error}",
                job.path.display(),
                job.info.sample
            )
        })?;
        records += 1;
        if let Some(edit) = edit {
            genes[gene].push(edit);
            applied += 1;
        }
    }
    if !seen_header || file_digest(&job.path)? != job.info.vcf_sha256 {
        return fail(format!(
            "VCF missing header or changed during encoding: {}",
            job.path.display()
        ));
    }
    for (sequence, edits) in reference.sequences.iter().zip(&mut genes) {
        edits.sort_by_key(|edit| edit.0);
        validate(sequence.dna.len(), edits).map_err(|error| {
            format!("sample {}, gene {}: {error}", job.info.sample, sequence.id)
        })?;
    }
    Ok(Calls {
        genes,
        records,
        applied,
    })
}

fn record(
    text: &str,
    job: &Job,
    reference: &Reference,
    haplotype: Option<usize>,
) -> Result<(usize, Option<Edit>)> {
    let fields: Vec<_> = text.split('\t').collect();
    if fields.len() != job.columns {
        return fail("VCF record column count differs from header");
    }
    let gene = *reference
        .index
        .get(fields[0])
        .ok_or_else(|| format!("unknown reference ID: {}", fields[0]))?;
    let sequence = &reference.sequences[gene].dna;
    let start = fields[1]
        .parse::<usize>()?
        .checked_sub(1)
        .ok_or("VCF POS must be >= 1")?;
    let ref_allele = fields[3].as_bytes();
    if ref_allele.is_empty() || !is_dna(ref_allele) {
        return fail("VCF REF must be literal IUPAC DNA");
    }
    let end = start
        .checked_add(ref_allele.len())
        .ok_or("VCF REF coordinate overflow")?;
    let actual = sequence
        .get(start..end)
        .ok_or("VCF REF interval exceeds reference")?;
    if !actual.eq_ignore_ascii_case(ref_allele) {
        return fail(format!(
            "VCF REF does not match FASTA at {}:{}",
            fields[0], fields[1]
        ));
    }
    if !matches!(fields[6], "." | "PASS") {
        return fail("nonaccepted FILTER; supply the final Snippy VCF, not raw calls");
    }
    let formats: Vec<_> = fields[8].split(':').collect();
    if formats.iter().collect::<HashSet<_>>().len() != formats.len() {
        return fail("duplicate VCF FORMAT key");
    }
    let gt = formats
        .iter()
        .position(|field| *field == "GT")
        .ok_or("VCF FORMAT has no GT")?;
    let sample_values: Vec<_> = fields[job.column].split(':').collect();
    if sample_values.len() > formats.len() {
        return fail("sample FORMAT has excess values");
    }
    let genotype = sample_values.get(gt).ok_or("sample GT missing")?;
    let alleles: Vec<_> = fields[4].split(',').collect();
    let selected = select_allele(genotype, alleles.len(), haplotype)?;
    if selected == 0 {
        return Ok((gene, None));
    }
    let alternate = alleles[selected - 1];
    if alternate.is_empty() || !is_dna(alternate.as_bytes()) {
        return fail("selected ALT is not literal DNA (symbolic, spanning-deletion or breakend); no exact sequence can be inferred");
    }
    let edit = replacement(start, actual, alternate);
    if let Some(ref edit) = edit {
        let local = Edit(edit.0 - start, edit.1, edit.2.clone());
        if apply(actual, &[local])? != alternate.as_bytes() {
            return fail("internal error: trimmed VCF replacement does not reconstruct ALT");
        }
    }
    Ok((gene, edit))
}

fn select_allele(genotype: &str, alt_count: usize, haplotype: Option<usize>) -> Result<usize> {
    if genotype.contains('/') && genotype.contains('|') {
        return fail("mixed genotype separators");
    }
    let values: Vec<_> = genotype
        .split(['/', '|'])
        .map(|value| {
            if value == "." || value.is_empty() {
                return fail("missing genotype; cannot define a lossless called consensus");
            }
            let allele = value.parse::<usize>()?;
            if allele > alt_count {
                return fail("genotype allele index exceeds ALT count");
            }
            Ok(allele)
        })
        .collect::<Result<_>>()?;
    let first = *values.first().ok_or("empty genotype")?;
    if values.iter().all(|value| *value == first) {
        return Ok(first);
    }
    if !genotype.contains('|') {
        return fail(
            "unphased mixed genotype has no unique sequence; use a resolved haplotype VCF",
        );
    }
    let selected = haplotype.ok_or("phased mixed genotype requires --haplotype 1 or 2")?;
    values
        .get(selected - 1)
        .copied()
        .ok_or_else(|| "selected haplotype exceeds genotype ploidy".into())
}

pub fn final_vcfs(directory: &Path) -> Result<Vec<PathBuf>> {
    let mut files = Vec::new();
    for entry in std::fs::read_dir(directory)? {
        let path = entry?.path();
        if path.is_dir() {
            files.extend(final_vcfs(&path)?);
        } else if path.is_file() {
            let name = path
                .file_name()
                .and_then(|name| name.to_str())
                .ok_or("non-UTF8 VCF filename")?;
            let base = name.strip_suffix(".gz").unwrap_or(name);
            if base.ends_with(".vcf")
                && ![".raw.vcf", ".filt.vcf", ".subs.vcf"]
                    .iter()
                    .any(|suffix| base.ends_with(suffix))
                && !(name.ends_with(".gz") && path.with_file_name(base).is_file())
            {
                files.push(path);
            }
        }
    }
    files.sort();
    if files.is_empty() {
        return fail(format!(
            "no final .vcf/.vcf.gz files in {}",
            directory.display()
        ));
    }
    Ok(files)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn genotype_selection_is_explicit() {
        assert_eq!(select_allele("1/1", 2, None).unwrap(), 1);
        assert_eq!(select_allele("2", 2, None).unwrap(), 2);
        assert_eq!(select_allele("0/0", 2, None).unwrap(), 0);
        assert_eq!(select_allele("1|2", 2, Some(2)).unwrap(), 2);
        for gt in ["0/1", "./.", "3/3", "1|2", "1/2|1", ""] {
            assert!(select_allele(gt, 2, None).is_err(), "{gt}");
        }
    }

    #[test]
    fn metadata_quotes_and_lengths() {
        assert_eq!(
            contig_metadata("<ID=gene,length=100,Description=\"a,b\">").unwrap(),
            ("gene".to_string(), Some(100))
        );
        assert!(contig_metadata("<ID=gene,ID=other>").is_err());
    }
}
