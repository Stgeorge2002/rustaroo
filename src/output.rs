use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::error::Result;
use crate::types::Pangenome;

const OUTPUT_BUFFER: usize = 1 << 20;

pub fn write_pangenome(pangenome: &Pangenome, directory: &Path) -> Result<()> {
    fs::create_dir_all(directory)?;
    if pangenome.membership_audit.is_none()
        && (directory.join("members.tsv").exists() || directory.join("member_alleles.fa").exists())
    {
        return Err(crate::error::Error::new(
            "membership audit files already exist; use a fresh output directory or --audit-members to avoid a stale audit",
        ));
    }
    let fasta_path = directory.join("genes.fa");
    let presence_path = directory.join("presence.tsv");
    let lengths_path = directory.join("lengths.tsv");
    let report_path = directory.join("reference_report.tsv");
    let ((fasta, presence), (lengths, report)) = rayon::join(
        || {
            rayon::join(
                || write_fasta(&fasta_path, pangenome),
                || write_presence(&presence_path, pangenome),
            )
        },
        || {
            rayon::join(
                || write_lengths(&lengths_path, pangenome),
                || write_report(&report_path, pangenome),
            )
        },
    );
    fasta?;
    presence?;
    lengths?;
    report?;
    if pangenome.membership_audit.is_some() {
        write_members(&directory.join("members.tsv"), pangenome)?;
    }
    Ok(())
}

fn buffered(path: &Path) -> Result<BufWriter<File>> {
    Ok(BufWriter::with_capacity(OUTPUT_BUFFER, File::create(path)?))
}

fn write_fasta(path: &Path, pangenome: &Pangenome) -> Result<()> {
    let mut out = buffered(path)?;
    for gene in &pangenome.genes {
        // The header is only the id. SnpEff anchors on the whole header line,
        // so a product after a space would become part of the contig name.
        writeln!(out, ">{}", gene.id)?;
        for chunk in gene.reference_dna.as_bytes().chunks(60) {
            writeln!(out, "{}", std::str::from_utf8(chunk).unwrap_or(""))?;
        }
    }
    Ok(())
}

fn write_presence(path: &Path, pangenome: &Pangenome) -> Result<()> {
    let mut out = buffered(path)?;
    write!(out, "gene")?;
    for sample in &pangenome.sample_names {
        write!(out, "\t{sample}")?;
    }
    writeln!(out)?;
    for gene in &pangenome.genes {
        write!(out, "{}", gene.id)?;
        for present in &gene.presence {
            write!(out, "\t{}", if *present { "1" } else { "0" })?;
        }
        writeln!(out)?;
    }
    Ok(())
}

fn write_lengths(path: &Path, pangenome: &Pangenome) -> Result<()> {
    let mut out = buffered(path)?;
    writeln!(out, "gene\tsample\tbases")?;
    for gene in &pangenome.genes {
        for (sample, &length) in pangenome.sample_names.iter().zip(&gene.lengths) {
            if length > 0 {
                writeln!(out, "{}\t{sample}\t{length}", gene.id)?;
            }
        }
    }
    Ok(())
}

fn write_report(path: &Path, pangenome: &Pangenome) -> Result<()> {
    let mut out = buffered(path)?;
    writeln!(
        out,
        "gene\tproduct\treference_bp\tintact\tselection\tsamples_present\tlength_evidence\tsupport\tmodal_end_bp\tmodal_end_samples\textension_samples\tsamples_with_exact_reference\tdistinct_alleles\tintact_alleles"
    )?;
    for gene in &pangenome.genes {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            gene.id,
            sanitize_field(&gene.product),
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
        )?;
    }
    Ok(())
}

fn sanitize_field(value: &str) -> String {
    value
        .chars()
        .map(|character| match character {
            '\t' | '\n' | '\r' => ' ',
            _ => character,
        })
        .collect()
}

fn write_members(path: &Path, pangenome: &Pangenome) -> Result<()> {
    let audit = pangenome
        .membership_audit
        .as_ref()
        .expect("membership audit requested");
    let mut out = buffered(path)?;
    writeln!(
        out,
        "gene\tsample\tcontig\tannotation_id\tstart\tend\tstrand\tlocus\torigin\tallele"
    )?;
    for member in &audit.members {
        let gene = &pangenome.genes[member.gene].id;
        let sample = &pangenome.sample_names[member.sample as usize];
        let contig = &audit.contig_names[member.sample as usize][member.contig as usize];
        for value in [gene, sample, contig, &member.annotation_id] {
            if value
                .chars()
                .any(|character| matches!(character, '\t' | '\n' | '\r'))
            {
                return Err(crate::error::Error::new(
                    "membership identifiers contain TSV delimiters; cannot write an unambiguous audit",
                ));
            }
        }
        writeln!(
            out,
            "{gene}\t{sample}\t{contig}\t{}\t{}\t{}\t{}\t{}\t{}\tdna_{}",
            member.annotation_id,
            member.start,
            member.end,
            member.strand,
            member.locus,
            member.origin,
            member.dna
        )?;
    }
    out.flush()?;
    let used: std::collections::BTreeSet<_> =
        audit.members.iter().map(|member| member.dna).collect();
    let mut alleles = buffered(&path.with_file_name("member_alleles.fa"))?;
    for id in used {
        writeln!(alleles, ">dna_{id}")?;
        for chunk in audit.sequences[id as usize].as_bytes().chunks(60) {
            alleles.write_all(chunk)?;
            writeln!(alleles)?;
        }
    }
    alleles.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{Gene, GeneMember, MembershipAudit};
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicUsize, Ordering};

    static NEXT: AtomicUsize = AtomicUsize::new(0);

    struct TemporaryOutput(PathBuf);

    impl TemporaryOutput {
        fn new() -> Self {
            let path = std::env::temp_dir().join(format!(
                "rustroo-members-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed),
            ));
            fs::create_dir(&path).unwrap();
            Self(path)
        }
    }

    impl Drop for TemporaryOutput {
        fn drop(&mut self) {
            fs::remove_dir_all(&self.0).expect("remove this test's temporary output");
        }
    }

    fn example() -> Pangenome {
        Pangenome {
            sample_names: vec!["sample".to_string()],
            genes: vec![Gene {
                id: "gene".to_string(),
                product: String::new(),
                reference_dna: "ATGAAATAA".to_string(),
                intact: true,
                selection: "intact",
                presence: vec![true],
                lengths: vec![9],
                samples_present: 1,
                length_evidence: "unconfirmed_short",
                samples_with_exact_reference: 1,
                distinct_alleles: 1,
                intact_alleles: 1,
                modal_end_bp: 9,
                modal_end_samples: 1,
                extension_samples: 0,
                support: "core",
            }],
            membership_audit: Some(MembershipAudit {
                contig_names: vec![vec!["contig".to_string()]],
                members: vec![GeneMember {
                    gene: 0,
                    sample: 0,
                    contig: 0,
                    annotation_id: "cds1".to_string(),
                    start: 100,
                    end: 108,
                    strand: -1,
                    locus: 7,
                    dna: 0,
                    origin: "annotated",
                }],
                sequences: vec!["ATGAAATAA".to_string(), "ATGCCCTAA".to_string()],
            }),
        }
    }

    #[test]
    fn writes_members_and_only_used_alleles() {
        let directory = TemporaryOutput::new();
        write_pangenome(&example(), &directory.0).unwrap();
        assert_eq!(fs::read_to_string(directory.0.join("reference_report.tsv")).unwrap(),
            "gene\tproduct\treference_bp\tintact\tselection\tsamples_present\tlength_evidence\tsupport\tmodal_end_bp\tmodal_end_samples\textension_samples\tsamples_with_exact_reference\tdistinct_alleles\tintact_alleles\n\
gene\t\t9\tyes\tintact\t1\tunconfirmed_short\tcore\t9\t1\t0\t1\t1\t1\n");
        assert_eq!(
            fs::read_to_string(directory.0.join("members.tsv")).unwrap(),
            "gene\tsample\tcontig\tannotation_id\tstart\tend\tstrand\tlocus\torigin\tallele\n\
gene\tsample\tcontig\tcds1\t100\t108\t-1\t7\tannotated\tdna_0\n"
        );
        assert_eq!(
            fs::read_to_string(directory.0.join("member_alleles.fa")).unwrap(),
            ">dna_0\nATGAAATAA\n"
        );
    }

    #[test]
    fn normal_output_omits_audit_and_refuses_stale_memberships() {
        let directory = TemporaryOutput::new();
        let mut pangenome = example();
        pangenome.membership_audit = None;
        write_pangenome(&pangenome, &directory.0).unwrap();
        assert!(!directory.0.join("members.tsv").exists());
        assert!(!directory.0.join("member_alleles.fa").exists());
        fs::write(directory.0.join("members.tsv"), "stale").unwrap();
        let error = write_pangenome(&pangenome, &directory.0).unwrap_err();
        assert!(error.to_string().contains("stale audit"));
    }

    #[test]
    fn audit_identifiers_are_not_silently_sanitized() {
        let directory = TemporaryOutput::new();
        let mut pangenome = example();
        pangenome.membership_audit.as_mut().unwrap().members[0].annotation_id =
            "ambiguous\tid".to_string();
        let error = write_pangenome(&pangenome, &directory.0).unwrap_err();
        assert!(error.to_string().contains("TSV delimiters"));
    }
}
