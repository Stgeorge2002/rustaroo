use std::path::PathBuf;

use clap::{Args, Parser, Subcommand};
use mutcode::{decode, encode, DecodeOptions, EncodeOptions, Result};

#[derive(Parser)]
#[command(
    version,
    about = "Lossless gene consensus codes from final Snippy VCFs"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Encode every reference gene for each selected VCF sample
    Encode(EncodeArgs),
    /// Verify and optionally reconstruct FASTAs from codes
    Decode(DecodeArgs),
}

#[derive(Args)]
struct EncodeArgs {
    #[arg(short, long)]
    reference: PathBuf,
    /// Recursively discover final .vcf/.vcf.gz, excluding raw/filt/subs duplicates
    #[arg(long, required_unless_present = "vcf", conflicts_with = "vcf")]
    snippy_dir: Option<PathBuf>,
    #[arg(long, num_args = 1.., required_unless_present = "snippy_dir")]
    vcf: Vec<PathBuf>,
    /// Fresh TSV output; .gz enables streaming gzip compression
    #[arg(short, long)]
    output: PathBuf,
    /// Select a sample by its exact VCF header name; default is all samples
    #[arg(long)]
    sample: Option<String>,
    /// Choose a haplotype for phased mixed genotypes; unphased mixtures fail
    #[arg(long, value_parser = clap::value_parser!(u8).range(1..=2))]
    haplotype: Option<u8>,
    #[arg(long, default_value_t = default_threads())]
    threads: usize,
    /// Independently supplied full-variant consensus; requires one selected sample
    #[arg(long)]
    check_consensus: Option<PathBuf>,
}

#[derive(Args)]
struct DecodeArgs {
    #[arg(short, long)]
    reference: PathBuf,
    #[arg(short, long)]
    codes: PathBuf,
    /// Fresh directory for reconstructed FASTAs and sample filename mapping
    #[arg(
        short,
        long,
        required_unless_present = "verify_only",
        conflicts_with = "verify_only"
    )]
    output: Option<PathBuf>,
    /// Verify every row and hash without writing FASTAs
    #[arg(long)]
    verify_only: bool,
    /// Write only this sample; all rows are still verified
    #[arg(long)]
    sample: Option<String>,
    #[arg(long, default_value_t = default_threads())]
    threads: usize,
}

fn default_threads() -> usize {
    std::thread::available_parallelism()
        .map_or(1, |count| count.get())
        .min(8)
}

fn main() {
    if let Err(error) = run() {
        eprintln!("mutcode: {error}");
        std::process::exit(1);
    }
}

fn run() -> Result<()> {
    let summary = match Cli::parse().command {
        Command::Encode(args) => encode(&EncodeOptions {
            reference: args.reference,
            vcfs: args.vcf,
            snippy_dir: args.snippy_dir,
            output: args.output,
            sample: args.sample,
            haplotype: args.haplotype.map(usize::from),
            threads: args.threads,
            check_consensus: args.check_consensus,
        })?,
        Command::Decode(args) => decode(&DecodeOptions {
            reference: args.reference,
            codes: args.codes,
            output: args.output,
            sample: args.sample,
            threads: args.threads,
        })?,
    };
    println!(
        "samples\t{}\ngenes\t{}\nverified_rows\t{}\nvcf_records\t{}\napplied_records\t{}",
        summary.samples,
        summary.genes,
        summary.rows,
        summary.variant_records,
        summary.applied_records
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encode_requires_one_input_source_and_accepts_haplotype() {
        assert!(
            Cli::try_parse_from(["mutcode", "encode", "-r", "ref.fa", "-o", "codes.tsv"]).is_err()
        );
        assert!(Cli::try_parse_from([
            "mutcode",
            "encode",
            "-r",
            "ref.fa",
            "-o",
            "codes.tsv",
            "--vcf",
            "sample.vcf",
            "--snippy-dir",
            "runs",
        ])
        .is_err());
        assert!(Cli::try_parse_from([
            "mutcode",
            "encode",
            "-r",
            "ref.fa",
            "-o",
            "codes.tsv",
            "--vcf",
            "sample.vcf",
            "--haplotype",
            "3",
        ])
        .is_err());
        let cli = Cli::try_parse_from([
            "mutcode",
            "encode",
            "-r",
            "ref.fa",
            "-o",
            "codes.tsv",
            "--snippy-dir",
            "runs",
            "--haplotype",
            "2",
        ])
        .unwrap();
        let Command::Encode(args) = cli.command else {
            panic!("expected encode");
        };
        assert_eq!(args.haplotype, Some(2));
    }

    #[test]
    fn decode_requires_output_or_verification() {
        assert!(
            Cli::try_parse_from(["mutcode", "decode", "-r", "ref.fa", "-c", "codes.tsv"]).is_err()
        );
        assert!(Cli::try_parse_from([
            "mutcode",
            "decode",
            "-r",
            "ref.fa",
            "-c",
            "codes.tsv",
            "--verify-only",
        ])
        .is_ok());
        assert!(Cli::try_parse_from([
            "mutcode",
            "decode",
            "-r",
            "ref.fa",
            "-c",
            "codes.tsv",
            "--verify-only",
            "-o",
            "decoded",
        ])
        .is_err());
    }
}
