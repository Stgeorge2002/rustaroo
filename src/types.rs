use crate::error::{Error, Result};

/// Knobs for homology: near-identical preclusters, then 70% families, with
/// fragments left for synteny rather than forced together by a loose length
/// cutoff. These defaults are not an equivalent of a Panaroo cleaning mode.
#[derive(Debug, Clone)]
pub struct Params {
    pub precluster_identity: f64,
    pub family_identity: f64,
    /// Minimum shorter/longer length ratio for two proteins to be homologs.
    pub length_ratio: f64,
    pub min_cds_bp: usize,
    pub refind: bool,
    /// `0` uses the default Rayon pool. Homology decisions stay in
    /// longest-first order; scoring of proteins that share one representative
    /// index runs concurrently.
    pub threads: usize,
    /// Emit original annotation memberships for an external concordance audit.
    pub audit_members: bool,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            precluster_identity: 0.98,
            family_identity: 0.70,
            length_ratio: 0.95,
            min_cds_bp: 90,
            refind: true,
            threads: 0,
            audit_members: false,
        }
    }
}

impl Params {
    pub fn validate(&self) -> Result<()> {
        if !(0.0 < self.family_identity && self.family_identity <= 1.0) {
            return Err(Error::new(
                "family identity must be between 0 and 1",
            ));
        }
        if !(self.family_identity <= self.precluster_identity
            && self.precluster_identity <= 1.0)
        {
            return Err(Error::new(
                "precluster identity must be between the family identity and 1",
            ));
        }
        if !(0.0 < self.length_ratio && self.length_ratio <= 1.0) {
            return Err(Error::new("length ratio must be between 0 and 1"));
        }
        if self.min_cds_bp < 6 {
            return Err(Error::new(
                "minimum CDS length must be at least 6 bp (start plus stop)",
            ));
        }
        Ok(())
    }
}

#[derive(Debug, Clone)]
pub struct Contig {
    pub name: String,
    pub sequence: String,
}

/// One annotated CDS, already oriented 5' to 3'. `start` and `end` are
/// 1-based inclusive positions on the contig as assembled, not on the
/// transcribed sequence.
#[derive(Debug, Clone)]
pub struct ParsedGene {
    /// Decoded GFF CDS ID; empty if the input supplied no ID.
    pub annotation_id: String,
    pub contig: usize,
    pub start: u32,
    pub end: u32,
    pub strand: i8,
    pub dna: String,
    pub protein: Vec<u8>,
    pub symbol: String,
    pub product: String,
}

#[derive(Debug, Clone)]
pub struct ParsedSample {
    pub name: String,
    pub contigs: Vec<Contig>,
    pub genes: Vec<ParsedGene>,
}

#[derive(Debug, Clone)]
pub struct Gene {
    pub id: String,
    pub product: String,
    pub reference_dna: String,
    pub intact: bool,
    /// Why this sequence was chosen. See `reference`.
    pub selection: &'static str,
    pub presence: Vec<bool>,
    /// Nucleotides covered by this gene in the sample. Zero means absent.
    pub lengths: Vec<u32>,
    pub samples_present: u32,
    /// `length_sufficient`, `conserved_short`, or `unconfirmed_short`.
    /// Below 250 bp one genome cannot vouch for a gene; the other two labels
    /// record whether the collection did.
    pub length_evidence: &'static str,
    pub samples_with_exact_reference: u32,
    pub distinct_alleles: u32,
    pub intact_alleles: u32,
    /// DNA length of the finished modal end, or the reference length when no
    /// shorter end qualified.
    pub modal_end_bp: u32,
    /// Samples whose intact, unclipped copies land on that end.
    pub modal_end_samples: u32,
    /// Samples carrying an intact allele longer than the finished end.
    pub extension_samples: u32,
}

#[derive(Debug, Clone)]
pub struct Pangenome {
    pub sample_names: Vec<String>,
    pub genes: Vec<Gene>,
    pub membership_audit: Option<MembershipAudit>,
}

#[derive(Debug, Clone)]
pub struct MembershipAudit {
    pub contig_names: Vec<Vec<String>>,
    pub members: Vec<GeneMember>,
    pub sequences: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct GeneMember {
    pub gene: usize,
    pub sample: u32,
    pub contig: u32,
    pub annotation_id: String,
    pub start: u32,
    pub end: u32,
    pub strand: i8,
    pub locus: u32,
    pub dna: u32,
    /// `annotated`, `fusion_piece`, or `refound`.
    pub origin: &'static str,
}
