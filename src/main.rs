use std::fs;
use std::path::PathBuf;

use clap::{Args, Parser, Subcommand};
use rayon::prelude::*;

use rustroo::{build, is_gff, parse_gff_file, validate, write_pangenome, Error, Params, Result, TreeMode, ValidationOptions};

#[derive(Parser)]
#[command(
    name = "rustroo",
    version,
    about = "Build a bacterial gene pangenome: one intact reference per gene, and which samples carry it",
    subcommand_negates_reqs = true,
    args_conflicts_with_subcommands = true
)]
struct Cli {
    #[command(subcommand)]
    command: Option<Command>,

    /// GFF3 files, directories of GFF3 files, or .txt lists of paths
    #[arg(short, long, num_args = 1.., required = true)]
    input: Vec<PathBuf>,

    /// Output directory (also members.tsv and member_alleles.fa with --audit-members)
    #[arg(short, long, required = true)]
    output_dir: Option<PathBuf>,

    /// Identity for the first, near-exact protein clustering step
    #[arg(long, default_value_t = 0.98)]
    precluster_identity: f64,

    /// Identity for joining preclusters into gene families
    #[arg(long, default_value_t = 0.70)]
    family_identity: f64,

    /// Minimum shorter/longer length ratio for a homology link
    #[arg(long, default_value_t = 0.95)]
    length_ratio: f64,

    /// Ignore CDS shorter than this many nucleotides
    #[arg(long, default_value_t = 90)]
    min_cds_bp: usize,

    /// Do not search intergenic gaps for a gene that the annotation missed
    #[arg(long)]
    no_refind: bool,

    /// Parser threads. 0 uses the default pool. Clustering stays serial.
    #[arg(long, default_value_t = 0)]
    threads: usize,

    /// Also write members.tsv with original CDS IDs, coordinates and physical loci
    #[arg(long)]
    audit_members: bool,
}

#[derive(Subcommand)]
enum Command {
    /// Independently inspect all member alleles, physical loci and diagnostic gene trees
    Validate(ValidateCli),
}

#[derive(Args)]
struct ValidateCli {
    /// Completed Rustroo output containing --audit-members exports
    #[arg(short, long)]
    input: PathBuf,
    /// Fresh validation directory; existing directories are never overwritten
    #[arg(short, long)]
    output_dir: PathBuf,
    /// Validation worker threads; 0 uses Rayon's default
    #[arg(long, default_value_t = 0)]
    threads: usize,
    /// Diagnostic tree scope; all also includes ordinary single-copy families
    #[arg(long, value_enum, default_value = "uncertain")]
    tree_mode: TreeMode,
    /// Deterministic nucleotide-site bootstrap replicates; 0 disables support estimation
    #[arg(long, default_value_t = 100)]
    bootstrap_replicates: usize,
    /// Maximum distinct DNA alleles per tree; 0 removes the limit
    #[arg(long, default_value_t = 128)]
    max_tree_alleles: usize,
    /// Maximum traceback matrix cells per nonidentical protein pair; 0 removes the limit
    #[arg(long, default_value_t = 5_000_000)]
    max_alignment_cells: usize,
    /// Minimum unambiguous nucleotide sites present in every tree allele
    #[arg(long, default_value_t = 90)]
    min_tree_sites: usize,
}

fn main() {
    if let Err(error) = run() {
        eprintln!("rustroo: {error}");
        std::process::exit(1);
    }
}

fn run() -> Result<()> {
    let cli = Cli::parse();
    if let Some(Command::Validate(args)) = cli.command {
        let summary = validate(&args.input, &args.output_dir, &ValidationOptions {
            threads: args.threads, tree_mode: args.tree_mode,
            bootstrap_replicates: args.bootstrap_replicates, max_tree_alleles: args.max_tree_alleles,
            max_alignment_cells: args.max_alignment_cells, min_tree_sites: args.min_tree_sites,
        })?;
        println!("validated_families\t{}", summary.families);
        println!("unassessed_alignments\t{}", summary.unassessed_alignments);
        println!("diagnostic_trees\t{}", summary.trees_written);
        println!("validation\t{}", args.output_dir.display());
        return Ok(());
    }
    let output_dir = cli.output_dir.as_deref().ok_or_else(|| Error::new("missing output directory"))?;
    let params = Params {
        precluster_identity: cli.precluster_identity,
        family_identity: cli.family_identity,
        length_ratio: cli.length_ratio,
        min_cds_bp: cli.min_cds_bp,
        refind: !cli.no_refind,
        threads: cli.threads,
        audit_members: cli.audit_members,
    };
    let paths = collect_inputs(&cli.input)?;
    let parser_threads = if params.threads == 0 {
        std::thread::available_parallelism()
            .map(|count| count.get())
            .unwrap_or(1)
    } else {
        params.threads
    };
    eprintln!(
        "rustroo: reading {} files on {parser_threads} threads; clustering stays on 1 thread",
        paths.len()
    );
    let samples = if params.threads == 0 {
        paths
            .par_iter()
            .map(|path| parse_gff_file(path))
            .collect::<Result<Vec<_>>>()?
    } else {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(params.threads)
            .build()?;
        pool.install(|| {
            paths
                .par_iter()
                .map(|path| parse_gff_file(path))
                .collect::<Result<Vec<_>>>()
        })?
    };
    let pangenome = build(samples, &params)?;
    write_pangenome(&pangenome, output_dir)?;
    println!("genes\t{}", pangenome.genes.len());
    println!("samples\t{}", pangenome.sample_names.len());
    println!("reference\t{}", output_dir.join("genes.fa").display());
    println!("presence\t{}", output_dir.join("presence.tsv").display());
    println!("lengths\t{}", output_dir.join("lengths.tsv").display());
    println!(
        "report\t{}",
        output_dir.join("reference_report.tsv").display()
    );
    if cli.audit_members {
        println!("members\t{}", output_dir.join("members.tsv").display());
    }
    let broken = pangenome
        .genes
        .iter()
        .filter(|gene| !gene.intact)
        .count();
    if broken > 0 {
        eprintln!(
            "rustroo: {broken} genes had no intact copy; their reference is the least-damaged allele"
        );
    }
    Ok(())
}

fn collect_inputs(inputs: &[PathBuf]) -> Result<Vec<PathBuf>> {
    let mut paths = Vec::new();
    for input in inputs {
        if input.is_dir() {
            let mut found = Vec::new();
            for entry in fs::read_dir(input)? {
                let path = entry?.path();
                if path.is_file() && is_gff(&path) {
                    found.push(path);
                }
            }
            found.sort();
            paths.extend(found);
        } else if input
            .extension()
            .is_some_and(|extension| extension.eq_ignore_ascii_case("txt"))
        {
            let text = fs::read_to_string(input)?;
            for line in text.lines() {
                let line = line.trim();
                if line.is_empty() || line.starts_with('#') {
                    continue;
                }
                paths.push(PathBuf::from(line));
            }
        } else if input.is_file() {
            paths.push(input.clone());
        } else {
            return Err(Error::new(format!("input not found: {}", input.display())));
        }
    }
    paths.sort();
    paths.dedup();
    if paths.is_empty() {
        Err(Error::new("no GFF files found"))
    } else {
        Ok(paths)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn build_cli_keeps_existing_required_arguments_and_defaults() {
        let cli = Cli::try_parse_from(["rustroo", "-i", "gff", "-o", "out", "--audit-members"]).unwrap();
        assert!(cli.command.is_none());
        assert_eq!(cli.output_dir, Some(PathBuf::from("out")));
        assert_eq!(cli.family_identity, 0.70);
        assert!(cli.audit_members);
        assert!(Cli::try_parse_from(["rustroo", "-i", "gff"]).is_err());
    }

    #[test]
    fn standalone_validation_does_not_require_build_inputs() {
        let cli = Cli::try_parse_from(["rustroo", "validate", "-i", "run", "-o", "validation", "--threads", "2"]).unwrap();
        let Some(Command::Validate(args)) = cli.command else { panic!("validation subcommand"); };
        assert_eq!(args.input, PathBuf::from("run"));
        assert_eq!(args.tree_mode, TreeMode::Uncertain);
        assert_eq!(args.threads, 2);
        assert_eq!(args.bootstrap_replicates, 100);
        assert!(Cli::try_parse_from(["rustroo", "-i", "gff", "-o", "out", "validate", "-i", "run", "-o", "validation"]).is_err());
    }
}
