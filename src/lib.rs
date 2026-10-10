//! Gene-level bacterial pangenome.
//!
//! Bakta or Prokka GFF3 in. One intact nucleotide reference per ortholog
//! group, and a matrix of which samples that gene was detected in. The
//! reference is the longest clean coding sequence, not the most common
//! allele, so a loss maps as a variant instead of disappearing into the
//! reference.

mod align;
mod assemble;
mod cluster;
mod error;
mod genetic_code;
mod gff;
mod model;
mod output;
mod phylogeny;
mod reference;
mod synteny;
mod types;
mod validation;
mod validation_align;

pub use assemble::build;
pub use error::{Error, Result};
pub use gff::{is_gff, parse_gff_file, sample_name};
pub use output::write_pangenome;
pub use types::{Gene, GeneMember, MembershipAudit, Pangenome, Params, ParsedGene, ParsedSample};
pub use validation::{validate, TreeMode, ValidationOptions, ValidationSummary};
