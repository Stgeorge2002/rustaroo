use std::fs;
use std::path::PathBuf;
use std::time::Instant;

use clap::{Args, Parser, Subcommand};
use rayon::prelude::*;

use rustroo::log;
use rustroo::{build, is_gff, parse_gff_file, validate, write_pangenome, Error, Params, Quality, Result, TreeMode, ValidationOptions};

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

    /// Worker threads. 0 uses the default pool. Homology decisions stay in longest-first order.
    #[arg(long, default_value_t = 0)]
    threads: usize,

    /// Also write members.tsv with original CDS IDs, coordinates and physical loci
    #[arg(long)]
    audit_members: bool,

    /// `high` drops a gene with no intact copy and an unconfirmed short gene.
    /// `strict` also drops a shorter or near-identical gene already covered by another.
    #[arg(long, value_enum, default_value = "standard")]
    quality: Quality,
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
    let timer = RunTimer::start();
    if let Err(error) = run() {
        log::record(&format!("rustroo: {error}"));
        timer.report();
        log::close();
        eprintln!("rustroo: {error}");
        std::process::exit(1);
    }
    timer.report();
    log::close();
}

/// Wall clock and CPU minutes for one command.
///
/// CPU time is this process's user and system time, summed across every
/// worker thread, plus children that have already been waited for. That is
/// the same core-time GNU time reports as user plus system, expressed in
/// minutes. Wall time is the actual elapsed time, so a parallel run can
/// spend more CPU minutes than minutes on the clock.
struct RunTimer {
    started: Instant,
    cpu_start: Option<f64>,
}

impl RunTimer {
    fn start() -> Self {
        Self {
            started: Instant::now(),
            cpu_start: process_cpu_seconds(),
        }
    }

    fn report(self) {
        let actual_min = format!("{:.3}", self.started.elapsed().as_secs_f64() / 60.0);
        let cpu_min = match (self.cpu_start, process_cpu_seconds()) {
            (Some(start), Some(end)) => format!("{:.3}", ((end - start).max(0.0)) / 60.0),
            _ => "na".to_string(),
        };
        println!("cpu_min\t{cpu_min}");
        println!("actual_min\t{actual_min}");
        log::record(&format!("result\tcpu_min\t{cpu_min}"));
        log::record(&format!("result\tactual_min\t{actual_min}"));
    }
}

/// Seconds of CPU already charged to this process.
///
/// `/proc/self/stat` is the aggregate of every thread. A per-thread clock
/// would under-count a Rayon run. The tick rate comes from `sysconf`, because
/// the proc fields are in ticks and that rate is not always 100.
fn process_cpu_seconds() -> Option<f64> {
    #[cfg(target_os = "linux")]
    {
        let text = fs::read_to_string("/proc/self/stat").ok()?;
        cpu_seconds_from_stat(&text, clock_ticks()?)
    }
    #[cfg(not(target_os = "linux"))]
    {
        None
    }
}

#[cfg(target_os = "linux")]
fn clock_ticks() -> Option<f64> {
    // glibc and musl both number `_SC_CLK_TCK` as 2 on Linux. The value
    // returned is the tick rate, not this name.
    const SC_CLK_TCK: i32 = 2;
    extern "C" {
        fn sysconf(name: i32) -> i64;
    }
    let ticks = unsafe { sysconf(SC_CLK_TCK) };
    (ticks > 0).then_some(ticks as f64)
}

/// User, system, and waited-for child ticks from one `/proc/self/stat` line.
///
/// The command name sits in parentheses and may itself contain spaces or
/// parentheses, so the fields after the last `)` are the stable ones.
#[cfg(any(target_os = "linux", test))]
fn cpu_seconds_from_stat(text: &str, ticks_per_second: f64) -> Option<f64> {
    if !ticks_per_second.is_finite() || ticks_per_second <= 0.0 {
        return None;
    }
    let rest = text.rsplit_once(')')?.1;
    let mut fields = rest.split_whitespace();
    // Field 3 is state. User time is field 14, then system, child user, child system.
    let utime: f64 = fields.nth(11)?.parse().ok()?;
    let stime: f64 = fields.next()?.parse().ok()?;
    let cutime: f64 = fields.next()?.parse().ok()?;
    let cstime: f64 = fields.next()?.parse().ok()?;
    let ticks = utime + stime + cutime + cstime;
    ticks.is_finite().then_some(ticks / ticks_per_second)
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
        quality: cli.quality,
    };
    let paths = collect_inputs(&cli.input)?;
    log::open(&output_dir.join("rustroo.log"))?;
    log::record(&format!("param\tprecluster_identity\t{}", params.precluster_identity));
    log::record(&format!("param\tfamily_identity\t{}", params.family_identity));
    log::record(&format!("param\tlength_ratio\t{}", params.length_ratio));
    log::record(&format!("param\tmin_cds_bp\t{}", params.min_cds_bp));
    log::record(&format!("param\trefind\t{}", params.refind));
    log::record(&format!("param\tthreads\t{}", params.threads));
    log::record(&format!("param\taudit_members\t{}", params.audit_members));
    log::record(&format!("param\tquality\t{}", quality_name(params.quality)));
    for path in &paths {
        log::record(&format!("input\t{}", path.display()));
    }
    let worker_threads = if params.threads == 0 {
        std::thread::available_parallelism()
            .map(|count| count.get())
            .unwrap_or(1)
    } else {
        params.threads
    };
    log::write_progress(&format!(
        "rustroo: reading {} files on {worker_threads} threads",
        paths.len()
    ));
    let compute = || -> Result<_> {
        let samples = paths
            .par_iter()
            .map(|path| parse_gff_file(path))
            .collect::<Result<Vec<_>>>()?;
        let pangenome = build(samples, &params)?;
        write_pangenome(&pangenome, output_dir)?;
        Ok(pangenome)
    };
    let pangenome = if params.threads == 0 {
        compute()?
    } else {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(params.threads)
            .build()?;
        pool.install(compute)?
    };
    println!("genes\t{}", pangenome.genes.len());
    println!("samples\t{}", pangenome.sample_names.len());
    println!("reference\t{}", output_dir.join("genes.fa").display());
    println!("presence\t{}", output_dir.join("presence.tsv").display());
    println!("lengths\t{}", output_dir.join("lengths.tsv").display());
    println!(
        "report\t{}",
        output_dir.join("reference_report.tsv").display()
    );
    println!("log\t{}", output_dir.join("rustroo.log").display());
    log::record(&format!("result\tgenes\t{}", pangenome.genes.len()));
    log::record(&format!("result\tsamples\t{}", pangenome.sample_names.len()));
    if cli.audit_members {
        println!("members\t{}", output_dir.join("members.tsv").display());
    }
    let broken = pangenome
        .genes
        .iter()
        .filter(|gene| !gene.intact)
        .count();
    if broken > 0 {
        log::write_progress(&format!(
            "rustroo: {broken} genes had no intact copy; their reference is the least-damaged allele"
        ));
    }
    log::flush();
    Ok(())
}

fn quality_name(quality: Quality) -> &'static str {
    match quality {
        Quality::Standard => "standard",
        Quality::High => "high",
        Quality::Strict => "strict",
    }
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
        assert_eq!(cli.quality, Quality::Standard);
        assert_eq!(
            Cli::try_parse_from(["rustroo", "-i", "gff", "-o", "out", "--quality", "high"]).unwrap().quality,
            Quality::High
        );
        assert_eq!(
            Cli::try_parse_from(["rustroo", "-i", "gff", "-o", "out", "--quality", "strict"]).unwrap().quality,
            Quality::Strict
        );
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

    #[test]
    fn stat_cpu_includes_this_process_and_waited_children() {
        let stat = "42 (rust roo) S 1 42 42 0 -1 0 1 2 3 4 10 20 30 40 5";
        let seconds = cpu_seconds_from_stat(stat, 100.0).unwrap();
        assert!((seconds - 1.0).abs() < 1e-12);
    }

    #[test]
    fn stat_cpu_ignores_parentheses_inside_the_command_name() {
        let stat = "7 (na)me) R 1 7 7 0 -1 0 0 0 0 0 50 0 0 0";
        let seconds = cpu_seconds_from_stat(stat, 100.0).unwrap();
        assert!((seconds - 0.5).abs() < 1e-12);
    }

    #[test]
    fn stat_cpu_rejects_a_truncated_or_untimed_line() {
        assert!(cpu_seconds_from_stat("1 (rustroo) S 1", 100.0).is_none());
        assert!(cpu_seconds_from_stat("1 (rustroo) S 1 1 1 0 -1 0 0 0 0 0 1 1 1 1", 0.0).is_none());
    }
}
