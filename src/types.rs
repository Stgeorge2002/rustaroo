use crate::error::{Error, Result};

/// What a finished family must show before it is written.
///
/// `standard` is the repaired pangenome. `high` drops the two weak labels:
/// a family with no intact copy, and a short family that is rare and not
/// between neighboring genes. `strict` also drops a second reference that is
/// already the same gene: a nested or deep fragment down to 80% identity, a
/// shorter piece down to 80%, a copy at 90% identity, or the same allele in
/// other genomes, which is folded into the longer gene. Repeated copies of
/// one gene fold even when they match each other. A fragment that matches two
/// different longer genes stays, and so does a paralogue below 90% identity
/// at nearly the same length. The reference stays the longer finished allele.
/// Every gene is also labelled core, accessory, unsupported, or refound. The
/// label is a column. `strict` does not drop a gene because of it.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum Quality {
    Standard,
    High,
    Strict,
}

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
    /// `strict` drops weak and duplicate references. `standard` keeps every
    /// repaired family. The command-line default is `strict`.
    pub quality: Quality,
    /// Fail on any input the lenient reader would repair or skip
    /// (see `InputFallbacks`), instead of only counting it.
    pub strict_input: bool,
    /// Split a family member allele off when its protein identity to the
    /// family representative is below this and the synteny evidence says it
    /// is a different locus. `0` disables the split.
    pub paralog_identity: f64,
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
            quality: Quality::Strict,
            strict_input: false,
            paralog_identity: 0.70,
        }
    }
}

impl Params {
    pub fn validate(&self) -> Result<()> {
        if !(0.0 < self.family_identity && self.family_identity <= 1.0) {
            return Err(Error::new("family identity must be between 0 and 1"));
        }
        if !(self.family_identity <= self.precluster_identity && self.precluster_identity <= 1.0) {
            return Err(Error::new(
                "precluster identity must be between the family identity and 1",
            ));
        }
        if !(0.0 < self.length_ratio && self.length_ratio <= 1.0) {
            return Err(Error::new("length ratio must be between 0 and 1"));
        }
        if !(0.0..=1.0).contains(&self.paralog_identity) {
            return Err(Error::new("paralog identity must be between 0 and 1"));
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
    /// Lenient-parsing decisions taken while reading this sample.
    pub fallbacks: InputFallbacks,
}

/// Every place the GFF3/FASTA reader silently repairs or skips input.
///
/// The default (lenient) behaviour is unchanged; these counts only make it
/// visible in `rustroo.log` and on stderr. `--strict-input` turns any non-zero
/// count into an error instead.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct InputFallbacks {
    /// Non-comment feature lines with fewer than 9 tab-separated columns (skipped).
    pub short_lines: u64,
    /// CDS whose seqid has no FASTA record (skipped).
    pub unknown_contig: u64,
    /// CDS whose start/end is not an integer, is 0, or has end < start (skipped).
    pub bad_coordinates: u64,
    /// CDS whose strand is not `+` or `-` (`.`, `?`, ...; read as `+`).
    pub unknown_strand: u64,
    /// CDS whose phase is not 0, 1 or 2 (unparsable read as 0, larger clamped to 2).
    pub bad_phase: u64,
    /// CDS parts ending past the end of their contig (clipped to the contig end).
    pub end_past_contig: u64,
    /// CDS dropped because a part starts past the end of its contig.
    pub part_outside_contig: u64,
    /// CDS whose length is not a multiple of 3 (trailing bases trimmed).
    pub trimmed_length: u64,
    /// CDS dropped because fewer than 6 bp remained.
    pub too_short: u64,
    /// CDS dropped because the translation was empty.
    pub empty_protein: u64,
    /// FASTA bases other than A, C, G, T, N (any case; `u` reads as T) read as N.
    pub non_acgt_bases: u64,
    /// Repeated FASTA record names (features resolve to the last record).
    pub duplicate_contigs: u64,
    /// The first offending input, for error messages.
    pub first: Option<String>,
}

impl InputFallbacks {
    /// `(log key, count)` for every kind, in a fixed order.
    pub fn counts(&self) -> [(&'static str, u64); 12] {
        [
            ("short_lines_skipped", self.short_lines),
            ("cds_unknown_contig_skipped", self.unknown_contig),
            ("cds_bad_coordinates_skipped", self.bad_coordinates),
            ("cds_unknown_strand_read_as_plus", self.unknown_strand),
            ("cds_bad_phase_defaulted", self.bad_phase),
            ("cds_end_clipped_to_contig", self.end_past_contig),
            ("cds_outside_contig_dropped", self.part_outside_contig),
            ("cds_length_trimmed_to_codons", self.trimmed_length),
            ("cds_too_short_dropped", self.too_short),
            ("cds_empty_protein_dropped", self.empty_protein),
            ("fasta_non_acgt_read_as_n", self.non_acgt_bases),
            ("fasta_duplicate_contig_names", self.duplicate_contigs),
        ]
    }

    pub fn total(&self) -> u64 {
        self.counts().iter().map(|(_, count)| count).sum()
    }

    pub fn add(&mut self, other: &InputFallbacks) {
        self.short_lines += other.short_lines;
        self.unknown_contig += other.unknown_contig;
        self.bad_coordinates += other.bad_coordinates;
        self.unknown_strand += other.unknown_strand;
        self.bad_phase += other.bad_phase;
        self.end_past_contig += other.end_past_contig;
        self.part_outside_contig += other.part_outside_contig;
        self.trimmed_length += other.trimmed_length;
        self.too_short += other.too_short;
        self.empty_protein += other.empty_protein;
        self.non_acgt_bases += other.non_acgt_bases;
        self.duplicate_contigs += other.duplicate_contigs;
        if self.first.is_none() {
            self.first = other.first.clone();
        }
    }

    pub(crate) fn note(&mut self, what: impl FnOnce() -> String) {
        if self.first.is_none() {
            self.first = Some(what());
        }
    }
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
    /// `core`, `accessory`, `unsupported`, or `refound`.
    ///
    /// `core` is present in at least 95% of genomes. `accessory` is rarer and
    /// the copies that exist share one neighborhood. `unsupported` is rarer
    /// and the neighborhoods do not agree. `refound` was not annotated; it
    /// was found in a gap. The label is a column on every gene. No quality
    /// mode drops a gene because of it.
    pub support: &'static str,
    /// Lowest protein identity of any member allele to the reference
    /// (over the shorter sequence; fragments are fitted, not globally aligned).
    pub min_allele_identity: f64,
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
