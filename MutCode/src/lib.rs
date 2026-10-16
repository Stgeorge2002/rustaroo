mod codec;
mod fasta;
mod io;
mod vcf;

use std::borrow::Cow;
use std::collections::HashSet;
use std::fs::{self, File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

pub use codec::{apply, Edit};
use fasta::write_record;
pub use fasta::Reference;
use io::{atomic_write, digest, file_digest, identifier, reader, staging_path};
use vcf::SampleInfo;

pub type Result<T> = std::result::Result<T, Box<dyn std::error::Error + Send + Sync>>;

fn fail<T>(message: impl Into<String>) -> Result<T> {
    Err(message.into().into())
}

const HEADER: &str =
    "sample\tgene\treference_length\treference_sha256\ttarget_length\ttarget_sha256\tcodes";
const ROW_BATCH_SIZE: usize = 512;

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct ReferenceInfo {
    gene: String,
    length: usize,
    sha256: String,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Metadata {
    format: String,
    target: String,
    coordinates: String,
    haplotype: Option<usize>,
    references: Vec<ReferenceInfo>,
    samples: Vec<SampleInfo>,
}

pub struct EncodeOptions {
    pub reference: PathBuf,
    pub vcfs: Vec<PathBuf>,
    pub snippy_dir: Option<PathBuf>,
    pub output: PathBuf,
    pub sample: Option<String>,
    pub haplotype: Option<usize>,
    pub threads: usize,
    pub check_consensus: Option<PathBuf>,
}

#[derive(Debug, PartialEq, Eq)]
pub struct Summary {
    pub samples: usize,
    pub genes: usize,
    pub rows: usize,
    pub variant_records: usize,
    pub applied_records: usize,
}

pub fn encode(options: &EncodeOptions) -> Result<Summary> {
    if options.threads == 0
        || options
            .haplotype
            .is_some_and(|value| !matches!(value, 1 | 2))
    {
        return fail("threads must be positive; haplotype must be 1 or 2");
    }
    if options.snippy_dir.is_some() != options.vcfs.is_empty() {
        return fail("supply exactly one of --snippy-dir or --vcf");
    }
    if options.output.exists() {
        return fail(format!(
            "output already exists: {}",
            options.output.display()
        ));
    }
    let reference_fingerprint = file_digest(&options.reference)?;
    let reference = Reference::load(&options.reference)?;
    let paths = match &options.snippy_dir {
        Some(directory) => vcf::final_vcfs(directory)?,
        None => options.vcfs.clone(),
    };
    let jobs = vcf::discover(&paths, &reference, options.sample.as_deref())?;
    if options.check_consensus.is_some() && jobs.len() != 1 {
        return fail("--check-consensus requires exactly one selected sample");
    }
    let consensus_hash = options
        .check_consensus
        .as_deref()
        .map(file_digest)
        .transpose()?;
    let consensus = options
        .check_consensus
        .as_deref()
        .map(Reference::load)
        .transpose()?;
    if let Some(consensus) = &consensus {
        if consensus.index.len() != reference.index.len()
            || reference
                .index
                .keys()
                .any(|id| !consensus.index.contains_key(id))
        {
            return fail("verification consensus must contain exactly the reference gene IDs");
        }
    }
    let metadata = Metadata {
        format: "mutcode-v1".into(),
        target: "called-consensus;unlisted=reference;FILTER=PASS-or-dot".into(),
        coordinates: "zero-based-original-reference;edit=[start,delete_length,inserted_sequence]"
            .into(),
        haplotype: options.haplotype,
        references: reference
            .sequences
            .iter()
            .map(|sequence| ReferenceInfo {
                gene: sequence.id.clone(),
                length: sequence.dna.len(),
                sha256: sequence.sha256.clone(),
            })
            .collect(),
        samples: jobs.iter().map(|job| job.info.clone()).collect(),
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(options.threads)
        .build()?;
    eprintln!(
        "mutcode: encoding with {} workers",
        pool.current_num_threads()
    );
    let mut records = 0;
    let mut applied = 0;
    atomic_write(&options.output, |writer| {
        writeln!(writer, "#mutcode\t{}", serde_json::to_string(&metadata)?)?;
        writeln!(writer, "{HEADER}")?;
        // Bound memory to one sample batch, preserving deterministic sample order.
        for batch in jobs.chunks(options.threads) {
            let calls = pool.install(|| {
                batch
                    .par_iter()
                    .map(|job| vcf::read_calls(job, &reference, options.haplotype))
                    .collect::<Result<Vec<_>>>()
            })?;
            for (job, calls) in batch.iter().zip(calls) {
                records += calls.records;
                applied += calls.applied;
                for (sequences, edits) in reference
                    .sequences
                    .chunks(ROW_BATCH_SIZE)
                    .zip(calls.genes.chunks(ROW_BATCH_SIZE))
                {
                    let rows = pool.install(|| {
                        sequences
                            .par_iter()
                            .zip(edits.par_iter())
                            .map(|(sequence, edits)| {
                                encode_row(sequence, edits, &job.info.sample, consensus.as_ref())
                            })
                            .collect::<Result<Vec<_>>>()
                    })?;
                    for row in rows {
                        writer.write_all(row.as_bytes())?;
                    }
                }
                eprintln!("mutcode: encoded {}", job.info.sample);
            }
        }
        if file_digest(&options.reference)? != reference_fingerprint {
            return fail("reference file changed during encoding");
        }
        if let (Some(path), Some(hash)) = (&options.check_consensus, &consensus_hash) {
            if &file_digest(path)? != hash {
                return fail("verification consensus changed during encoding");
            }
        }
        Ok(())
    })?;
    Ok(Summary {
        samples: jobs.len(),
        genes: reference.sequences.len(),
        rows: jobs
            .len()
            .checked_mul(reference.sequences.len())
            .ok_or("row count overflow")?,
        variant_records: records,
        applied_records: applied,
    })
}

fn encode_row(
    sequence: &fasta::Sequence,
    edits: &[Edit],
    sample: &str,
    consensus: Option<&Reference>,
) -> Result<String> {
    let target = if edits.is_empty() {
        Cow::Borrowed(sequence.dna.as_slice())
    } else {
        Cow::Owned(apply(&sequence.dna, edits)?)
    };
    let codes = if edits.is_empty() {
        "[]".to_string()
    } else {
        let codes = serde_json::to_string(edits)?;
        let decoded: Vec<Edit> = serde_json::from_str(&codes)?;
        if apply(&sequence.dna, &decoded)?.as_slice() != target.as_ref() {
            return fail("internal error: serialized codes fail byte-for-byte reconstruction");
        }
        codes
    };
    if let Some(consensus) = consensus {
        let expected = &consensus.sequences[consensus.index[&sequence.id]].dna;
        if expected.as_slice() != target.as_ref() {
            return fail(format!(
                "consensus mismatch for sample {sample}, gene {}",
                sequence.id
            ));
        }
    }
    let target_hash = if edits.is_empty() {
        sequence.sha256.clone()
    } else {
        digest(target.as_ref())
    };
    Ok(format!(
        "{sample}\t{}\t{}\t{}\t{}\t{}\t{codes}\n",
        sequence.id,
        sequence.dna.len(),
        sequence.sha256,
        target.len(),
        target_hash
    ))
}

pub struct DecodeOptions {
    pub reference: PathBuf,
    pub codes: PathBuf,
    pub output: Option<PathBuf>,
    pub sample: Option<String>,
    pub threads: usize,
}

pub fn decode(options: &DecodeOptions) -> Result<Summary> {
    if options.threads == 0 {
        return fail("threads must be positive");
    }
    if let Some(path) = &options.output {
        if path.exists() {
            return fail(format!("output already exists: {}", path.display()));
        }
    }
    let reference = Reference::load(&options.reference)?;
    let mut input = reader(&options.codes)?;
    let mut line = String::new();
    input.read_line(&mut line)?;
    let json = line
        .trim_end_matches(['\r', '\n'])
        .strip_prefix("#mutcode\t")
        .ok_or("missing MutCode format header")?;
    let metadata: Metadata = serde_json::from_str(json)?;
    if metadata.format != "mutcode-v1"
        || metadata.target != "called-consensus;unlisted=reference;FILTER=PASS-or-dot"
        || metadata.coordinates
            != "zero-based-original-reference;edit=[start,delete_length,inserted_sequence]"
        || metadata
            .haplotype
            .is_some_and(|value| !matches!(value, 1 | 2))
    {
        return fail("unsupported MutCode format or reconstruction policy");
    }
    if metadata.references.len() != reference.sequences.len() || metadata.samples.is_empty() {
        return fail("reference count mismatch or empty sample list");
    }
    let mut references = Vec::new();
    let mut ids = HashSet::new();
    for info in &metadata.references {
        let sequence = reference
            .index
            .get(&info.gene)
            .map(|index| &reference.sequences[*index])
            .ok_or_else(|| format!("missing reference: {}", info.gene))?;
        if sequence.dna.len() != info.length
            || sequence.sha256 != info.sha256
            || !ids.insert(&info.gene)
        {
            return fail(format!(
                "reference hash, length or uniqueness mismatch: {}",
                info.gene
            ));
        }
        references.push(sequence);
    }
    let mut samples = HashSet::new();
    for info in &metadata.samples {
        identifier(&info.sample)?;
        if !samples.insert(&info.sample)
            || info.vcf_sha256.len() != 64
            || !info.vcf_sha256.bytes().all(|byte| byte.is_ascii_hexdigit())
        {
            return fail("duplicate sample or invalid VCF fingerprint in codes");
        }
    }
    if options
        .sample
        .as_ref()
        .is_some_and(|sample| !samples.contains(sample))
    {
        return fail("requested sample not present in codes");
    }
    line.clear();
    input.read_line(&mut line)?;
    if line.trim_end_matches(['\r', '\n']) != HEADER {
        return fail("unexpected MutCode table columns");
    }
    let staging = options.output.as_deref().map(staging_path).transpose()?;
    if let Some(path) = &staging {
        fs::create_dir(path)?;
    }
    let result = decode_rows(
        options,
        &metadata,
        &references,
        &mut *input,
        staging.as_deref(),
    );
    match result {
        Ok(summary) => {
            if let (Some(staging), Some(output)) = (&staging, &options.output) {
                // Reserve a destination rather than overwriting a competing output.
                fs::create_dir(output)?;
                for entry in fs::read_dir(staging)? {
                    let path = entry?.path();
                    fs::rename(
                        &path,
                        output.join(path.file_name().ok_or("missing output filename")?),
                    )?;
                }
                fs::remove_dir(staging)?;
            }
            Ok(summary)
        }
        Err(error) => {
            if let Some(path) = &staging {
                // Keep explicitly named failure artifacts for inspection.
                eprintln!(
                    "mutcode: decoding failed; partial output retained at {}",
                    path.display()
                );
            }
            Err(error)
        }
    }
}

fn decode_rows(
    options: &DecodeOptions,
    metadata: &Metadata,
    references: &[&fasta::Sequence],
    input: &mut dyn std::io::BufRead,
    staging: Option<&Path>,
) -> Result<Summary> {
    let mut row = String::new();
    let mut written = 0;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(options.threads)
        .build()?;
    eprintln!(
        "mutcode: decoding/verifying with {} workers",
        pool.current_num_threads()
    );
    let mut map = match staging {
        Some(path) => Some(BufWriter::new(File::create(path.join("samples.tsv"))?)),
        None => None,
    };
    if let Some(map) = &mut map {
        writeln!(map, "sample\tfasta")?;
    }
    for (sample_index, info) in metadata.samples.iter().enumerate() {
        let selected = options
            .sample
            .as_ref()
            .is_none_or(|sample| *sample == info.sample);
        let mut output = if selected {
            if let Some(path) = staging {
                let filename = format!("s{:06}.fa", sample_index + 1);
                if let Some(map) = &mut map {
                    writeln!(map, "{}\t{filename}", info.sample)?;
                }
                Some(BufWriter::new(
                    OpenOptions::new()
                        .write(true)
                        .create_new(true)
                        .open(path.join(filename))?,
                ))
            } else {
                None
            }
        } else {
            None
        };
        for batch in references.chunks(ROW_BATCH_SIZE) {
            let mut rows = Vec::with_capacity(batch.len());
            for _ in batch {
                row.clear();
                if input.read_line(&mut row)? == 0 {
                    return fail("truncated codes table: missing sample/gene row");
                }
                rows.push(std::mem::take(&mut row));
            }
            let targets = pool.install(|| {
                rows.par_iter()
                    .zip(batch.par_iter())
                    .map(|(row, sequence)| decode_row(row, &info.sample, sequence))
                    .collect::<Result<Vec<_>>>()
            })?;
            if let Some(writer) = &mut output {
                for (sequence, target) in batch.iter().zip(targets) {
                    write_record(writer, &sequence.id, &target)?;
                }
            }
            written += batch.len();
        }

        fn decode_row<'a>(
            row: &str,
            sample: &str,
            sequence: &'a fasta::Sequence,
        ) -> Result<Cow<'a, [u8]>> {
            let fields: Vec<_> = row.trim_end_matches(['\r', '\n']).split('\t').collect();
            if fields.len() != 7 || fields[0] != sample || fields[1] != sequence.id {
                return fail("missing, duplicated, unexpected or reordered sample/gene row");
            }
            if fields[2].parse::<usize>()? != sequence.dna.len() || fields[3] != sequence.sha256 {
                return fail("row reference length or hash mismatch");
            }
            let edits: Vec<Edit> = serde_json::from_str(fields[6])?;
            let target = if edits.is_empty() {
                Cow::Borrowed(sequence.dna.as_slice())
            } else {
                Cow::Owned(apply(&sequence.dna, &edits)?)
            };
            let target_hash = if edits.is_empty() {
                sequence.sha256.clone()
            } else {
                digest(target.as_ref())
            };
            if fields[4].parse::<usize>()? != target.len() || fields[5] != target_hash {
                return fail(format!(
                    "target reconstruction hash/length mismatch: {sample} {}",
                    sequence.id
                ));
            }
            Ok(target)
        }
        if let Some(mut writer) = output {
            writer.flush()?;
            writer.get_ref().sync_all()?;
        }
    }
    row.clear();
    if input.read_line(&mut row)? != 0 {
        return fail("unexpected trailing rows in codes table");
    }
    if let Some(mut map) = map {
        map.flush()?;
        map.get_ref().sync_all()?;
    }
    Ok(Summary {
        samples: metadata.samples.len(),
        genes: references.len(),
        rows: written,
        variant_records: 0,
        applied_records: 0,
    })
}
