//! Bakta and Prokka GFF3.
//!
//! The CDS is taken from the contig FASTA, not from the `translation`
//! attribute, so the nucleotide reference and the protein used for clustering
//! are the same sequence. Minus-strand features are reverse-complemented into
//! the 5' to 3' coding sequence. A `##FASTA` block is the usual place the
//! contig lives; a sibling FASTA is accepted when the GFF is features only.

use crate::fxhash::HashMap;
use std::fs::File;
use std::io::Read;
use std::path::{Path, PathBuf};

use flate2::read::GzDecoder;

use crate::error::{Error, Result};
use crate::genetic_code::{reverse_complement, translate, uppercase_dna};
use crate::types::{Contig, ParsedGene, ParsedSample};

pub fn parse_gff_file(path: &Path) -> Result<ParsedSample> {
    let text = read_maybe_gzip(path)?;
    let name = sample_name(path);
    let sibling = sibling_fasta(path).and_then(|fasta| read_maybe_gzip(&fasta).ok());
    parse_gff(&name, &text, sibling.as_deref())
        .map_err(|error| Error::new(format!("{}: {error}", path.display())))
}

pub fn sample_name(path: &Path) -> String {
    let mut name = path
        .file_name()
        .map(|value| value.to_string_lossy().into_owned())
        .unwrap_or_else(|| "sample".to_string());
    if let Some(stripped) = name.strip_suffix(".gz") {
        name = stripped.to_string();
    }
    for extension in [".gff3", ".gff", ".fna", ".fasta", ".fa"] {
        if let Some(stripped) = name.strip_suffix(extension) {
            name = stripped.to_string();
            break;
        }
    }
    if name.is_empty() {
        "sample".to_string()
    } else {
        name
    }
}

/// Lean parser: identical `ParsedSample` to `parse_gff_reference` (kept below
/// and compared in tests), without a `Vec` per line, a `HashMap` per CDS, or a
/// `String` per contig piece.
pub fn parse_gff(name: &str, text: &str, sibling_fasta: Option<&str>) -> Result<ParsedSample> {
    let (features, embedded) = split_fasta(text);
    let fasta_text = embedded.or(sibling_fasta).ok_or_else(|| {
        Error::new("no contig FASTA. Expected a ##FASTA section or a sibling .fa/.fna/.fasta file")
    })?;
    let contigs = parse_fasta_fast(fasta_text);
    if contigs.is_empty() {
        return Err(Error::new("FASTA section contained no sequences"));
    }
    let mut by_name: HashMap<&str, usize> = HashMap::with_capacity_and_hasher(contigs.len(), Default::default());
    for (index, contig) in contigs.iter().enumerate() {
        by_name.insert(contig.name.as_str(), index);
    }

    let mut parts: Vec<CdsPart> = Vec::new();
    for line in features.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        // The first nine tab fields, exactly as `split('\t').collect()` then
        // indexing 0..=8 would see them; anything past column 9 is ignored.
        let mut columns = [""; 9];
        let mut count = 0usize;
        for field in line.split('\t') {
            columns[count] = field;
            count += 1;
            if count == 9 {
                break;
            }
        }
        if count < 9 || columns[2] != "CDS" {
            continue;
        }
        let Some(&contig_index) = by_name.get(columns[0]) else {
            continue;
        };
        let Ok(start) = columns[3].parse::<u32>() else {
            continue;
        };
        let Ok(end) = columns[4].parse::<u32>() else {
            continue;
        };
        if start == 0 || end < start {
            continue;
        }
        let strand = if columns[6] == "-" { -1 } else { 1 };
        let phase = columns[7].parse::<usize>().unwrap_or(0).min(2);
        let attrs = needed_attributes(columns[8]);
        let decoded = |value: Option<&str>| value.map(decode_value);
        parts.push(CdsPart {
            contig: contig_index,
            start,
            end,
            strand,
            phase,
            id: decoded(attrs.id).unwrap_or_default(),
            symbol: decoded(attrs.gene)
                .filter(|value| !value.is_empty())
                .or_else(|| decoded(attrs.name))
                .unwrap_or_default(),
            product: decoded(attrs.product).unwrap_or_default(),
        });
    }

    let genes = assemble_genes_with(&contigs, parts, gene_from_parts_fast)?;
    Ok(ParsedSample {
        name: name.to_string(),
        contigs,
        genes,
    })
}

/// The original parser, unchanged; the reference for equivalence tests.
#[allow(dead_code)]
pub(crate) fn parse_gff_reference(name: &str, text: &str, sibling_fasta: Option<&str>) -> Result<ParsedSample> {
    let (features, embedded) = split_fasta(text);
    let fasta_text = embedded.or(sibling_fasta).ok_or_else(|| {
        Error::new("no contig FASTA. Expected a ##FASTA section or a sibling .fa/.fna/.fasta file")
    })?;
    let contigs = parse_fasta(fasta_text);
    if contigs.is_empty() {
        return Err(Error::new("FASTA section contained no sequences"));
    }
    let mut by_name: HashMap<&str, usize> = HashMap::default();
    for (index, contig) in contigs.iter().enumerate() {
        by_name.insert(contig.name.as_str(), index);
    }

    let mut parts: Vec<CdsPart> = Vec::new();
    for line in features.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() < 9 || columns[2] != "CDS" {
            continue;
        }
        let Some(&contig_index) = by_name.get(columns[0]) else {
            continue;
        };
        let Ok(start) = columns[3].parse::<u32>() else {
            continue;
        };
        let Ok(end) = columns[4].parse::<u32>() else {
            continue;
        };
        if start == 0 || end < start {
            continue;
        }
        let strand = if columns[6] == "-" { -1 } else { 1 };
        let phase = columns[7].parse::<usize>().unwrap_or(0).min(2);
        let attrs = attributes(columns[8]);
        parts.push(CdsPart {
            contig: contig_index,
            start,
            end,
            strand,
            phase,
            id: attrs.get("ID").cloned().unwrap_or_default(),
            symbol: attrs
                .get("gene")
                .cloned()
                .filter(|value| !value.is_empty())
                .or_else(|| attrs.get("Name").cloned())
                .unwrap_or_default(),
            product: attrs.get("product").cloned().unwrap_or_default(),
        });
    }

    let genes = assemble_genes(&contigs, parts)?;
    Ok(ParsedSample {
        name: name.to_string(),
        contigs,
        genes,
    })
}

struct CdsPart {
    contig: usize,
    start: u32,
    end: u32,
    strand: i8,
    phase: usize,
    id: String,
    symbol: String,
    product: String,
}

fn assemble_genes(contigs: &[Contig], parts: Vec<CdsPart>) -> Result<Vec<ParsedGene>> {
    assemble_genes_with(contigs, parts, gene_from_parts)
}

fn assemble_genes_with(
    contigs: &[Contig],
    mut parts: Vec<CdsPart>,
    gene_from_parts: fn(&[Contig], &[CdsPart]) -> Option<ParsedGene>,
) -> Result<Vec<ParsedGene>> {
    parts.sort_by(|a, b| {
        a.contig
            .cmp(&b.contig)
            .then(a.id.cmp(&b.id))
            .then(a.start.cmp(&b.start))
    });
    let mut genes = Vec::new();
    let mut index = 0;
    while index < parts.len() {
        let mut group_end = index + 1;
        if !parts[index].id.is_empty() {
            while group_end < parts.len()
                && parts[group_end].contig == parts[index].contig
                && parts[group_end].id == parts[index].id
                && parts[group_end].strand == parts[index].strand
            {
                group_end += 1;
            }
        }
        let group = &parts[index..group_end];
        if let Some(gene) = gene_from_parts(contigs, group) {
            genes.push(gene);
        }
        index = group_end;
    }
    if genes.is_empty() {
        return Err(Error::new("no CDS features with sequence were found"));
    }
    Ok(genes)
}

#[allow(dead_code)]
fn gene_from_parts(contigs: &[Contig], parts: &[CdsPart]) -> Option<ParsedGene> {
    let strand = parts[0].strand;
    let mut ordered: Vec<&CdsPart> = parts.iter().collect();
    if strand < 0 {
        ordered.sort_by(|a, b| b.start.cmp(&a.start));
    } else {
        ordered.sort_by(|a, b| a.start.cmp(&b.start));
    }
    let mut dna = String::new();
    for (piece_index, part) in ordered.iter().enumerate() {
        let contig = &contigs[part.contig].sequence;
        let start = part.start.saturating_sub(1) as usize;
        let end = (part.end as usize).min(contig.len());
        if start >= end || end > contig.len() {
            return None;
        }
        let mut piece = contig[start..end].to_string();
        if strand < 0 {
            piece = reverse_complement(&piece);
        }
        if piece_index == 0 && part.phase > 0 && part.phase < piece.len() {
            piece = piece[part.phase..].to_string();
        }
        dna.push_str(&piece);
    }
    if dna.len() % 3 != 0 {
        dna.truncate(dna.len() - dna.len() % 3);
    }
    if dna.len() < 6 {
        return None;
    }
    let protein = translate(&dna);
    if protein.is_empty() {
        return None;
    }
    let first = ordered[0];
    let start = parts
        .iter()
        .map(|part| part.start)
        .min()
        .unwrap_or(first.start);
    let end = parts.iter().map(|part| part.end).max().unwrap_or(first.end);
    let symbol = parts
        .iter()
        .find(|part| !part.symbol.is_empty())
        .map(|part| part.symbol.clone())
        .unwrap_or_default();
    let product = parts
        .iter()
        .find(|part| !part.product.is_empty())
        .map(|part| part.product.clone())
        .unwrap_or_default();
    Some(ParsedGene {
        annotation_id: first.id.clone(),
        contig: first.contig,
        start,
        end,
        strand,
        dna,
        protein,
        symbol,
        product,
    })
}

/// Same bytes as `gene_from_parts`, built in one buffer: a minus-strand piece
/// is reverse-complemented straight from the contig slice, and the phase
/// offset is skipped instead of re-allocating the piece.
fn gene_from_parts_fast(contigs: &[Contig], parts: &[CdsPart]) -> Option<ParsedGene> {
    let strand = parts[0].strand;
    let mut ordered: Vec<&CdsPart> = parts.iter().collect();
    if strand < 0 {
        ordered.sort_by(|a, b| b.start.cmp(&a.start));
    } else {
        ordered.sort_by(|a, b| a.start.cmp(&b.start));
    }
    let mut dna: Vec<u8> = Vec::new();
    for (piece_index, part) in ordered.iter().enumerate() {
        let contig = &contigs[part.contig].sequence;
        let start = part.start.saturating_sub(1) as usize;
        let end = (part.end as usize).min(contig.len());
        if start >= end || end > contig.len() {
            return None;
        }
        // `str` slicing keeps the original's char-boundary behaviour.
        let piece = contig[start..end].as_bytes();
        let skip = if piece_index == 0 && part.phase > 0 && part.phase < piece.len() {
            part.phase
        } else {
            0
        };
        if strand < 0 {
            crate::genetic_code::extend_reverse_complement(&mut dna, piece, skip);
        } else {
            dna.extend_from_slice(&piece[skip..]);
        }
    }
    if dna.len() % 3 != 0 {
        dna.truncate(dna.len() - dna.len() % 3);
    }
    if dna.len() < 6 {
        return None;
    }
    let protein = crate::genetic_code::translate_bytes(&dna);
    if protein.is_empty() {
        return None;
    }
    // Contig bytes are ASCII after `uppercase_dna`, and so is every
    // complement, so this cannot fail; the lossy path is never taken.
    let dna = String::from_utf8(dna)
        .unwrap_or_else(|error| String::from_utf8_lossy(error.as_bytes()).into_owned());
    let first = ordered[0];
    let start = parts
        .iter()
        .map(|part| part.start)
        .min()
        .unwrap_or(first.start);
    let end = parts.iter().map(|part| part.end).max().unwrap_or(first.end);
    let symbol = parts
        .iter()
        .find(|part| !part.symbol.is_empty())
        .map(|part| part.symbol.clone())
        .unwrap_or_default();
    let product = parts
        .iter()
        .find(|part| !part.product.is_empty())
        .map(|part| part.product.clone())
        .unwrap_or_default();
    Some(ParsedGene {
        annotation_id: first.id.clone(),
        contig: first.contig,
        start,
        end,
        strand,
        dna,
        protein,
        symbol,
        product,
    })
}

/// Raw (undecoded) values of the four attributes the parser reads. A repeated
/// key keeps its last value, as `HashMap::insert` in `attributes` does.
#[derive(Default)]
struct NeededAttributes<'a> {
    id: Option<&'a str>,
    gene: Option<&'a str>,
    name: Option<&'a str>,
    product: Option<&'a str>,
}

fn needed_attributes(field: &str) -> NeededAttributes<'_> {
    let mut out = NeededAttributes::default();
    for part in field.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        let Some((key, value)) = part.split_once('=') else {
            continue;
        };
        match key.trim() {
            "ID" => out.id = Some(value),
            "gene" => out.gene = Some(value),
            "Name" => out.name = Some(value),
            "product" => out.product = Some(value),
            _ => {}
        }
    }
    out
}

/// The value transformation `attributes` applies, with a fast path: without a
/// `%`, `percent_decode` copies the bytes of an already valid `str`.
fn decode_value(value: &str) -> String {
    let value = value.trim().trim_matches('"');
    if value.as_bytes().contains(&b'%') {
        percent_decode(value)
    } else {
        value.to_string()
    }
}

/// `parse_fasta` with the case/U->N mapping applied while copying each line,
/// instead of concatenating and then re-collecting through `char`s.
fn parse_fasta_fast(text: &str) -> Vec<Contig> {
    let mut contigs = Vec::new();
    let mut name = String::new();
    let mut sequence: Vec<u8> = Vec::new();
    let flush = |name: &mut String, sequence: &mut Vec<u8>, contigs: &mut Vec<Contig>| {
        if name.is_empty() {
            return;
        }
        let bytes = std::mem::take(sequence);
        contigs.push(Contig {
            name: std::mem::take(name),
            // Every byte went through `uppercase_dna_byte`, so this is ASCII.
            sequence: String::from_utf8(bytes)
                .unwrap_or_else(|error| String::from_utf8_lossy(error.as_bytes()).into_owned()),
        });
    };
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            flush(&mut name, &mut sequence, &mut contigs);
            name = header
                .split_whitespace()
                .next()
                .unwrap_or("contig")
                .to_string();
        } else if !name.is_empty() {
            sequence.extend(line.bytes().map(crate::genetic_code::uppercase_dna_byte));
        }
    }
    flush(&mut name, &mut sequence, &mut contigs);
    contigs
}

#[allow(dead_code)]
fn attributes(field: &str) -> HashMap<String, String> {
    let mut map = HashMap::default();
    for part in field.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        let Some((key, value)) = part.split_once('=') else {
            continue;
        };
        let value = percent_decode(value.trim().trim_matches('"'));
        map.insert(key.trim().to_string(), value);
    }
    map
}

fn percent_decode(value: &str) -> String {
    let bytes = value.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut index = 0;
    while index < bytes.len() {
        if bytes[index] == b'%' && index + 2 < bytes.len() {
            if let Ok(byte) = u8::from_str_radix(
                std::str::from_utf8(&bytes[index + 1..index + 3]).unwrap_or(""),
                16,
            ) {
                out.push(byte);
                index += 3;
                continue;
            }
        }
        out.push(bytes[index]);
        index += 1;
    }
    String::from_utf8_lossy(&out).into_owned()
}

fn split_fasta(text: &str) -> (&str, Option<&str>) {
    for marker in ["\n##FASTA", "\n##fasta"] {
        if let Some(position) = text.find(marker) {
            let features = &text[..position];
            let fasta = text[position + marker.len()..].trim_start_matches(['\r', '\n']);
            return (features, Some(fasta));
        }
    }
    if text.starts_with("##FASTA") || text.starts_with("##fasta") {
        let rest = text
            .find('\n')
            .map(|index| &text[index + 1..])
            .unwrap_or("");
        return ("", Some(rest));
    }
    (text, None)
}

#[allow(dead_code)]
fn parse_fasta(text: &str) -> Vec<Contig> {
    let mut contigs = Vec::new();
    let mut name = String::new();
    let mut sequence = String::new();
    let flush = |name: &mut String, sequence: &mut String, contigs: &mut Vec<Contig>| {
        if name.is_empty() {
            return;
        }
        contigs.push(Contig {
            name: std::mem::take(name),
            sequence: uppercase_dna(sequence),
        });
        sequence.clear();
    };
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            flush(&mut name, &mut sequence, &mut contigs);
            name = header
                .split_whitespace()
                .next()
                .unwrap_or("contig")
                .to_string();
        } else if !name.is_empty() {
            sequence.push_str(line);
        }
    }
    flush(&mut name, &mut sequence, &mut contigs);
    contigs
}

fn sibling_fasta(path: &Path) -> Option<PathBuf> {
    let file_name = path.file_name()?.to_string_lossy();
    let mut stem = file_name.as_ref().to_string();
    if let Some(stripped) = stem.strip_suffix(".gz") {
        stem = stripped.to_string();
    }
    for extension in [".gff3", ".gff"] {
        if let Some(stripped) = stem.strip_suffix(extension) {
            stem = stripped.to_string();
            break;
        }
    }
    let parent = path.parent().unwrap_or_else(|| Path::new("."));
    for extension in [".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz"] {
        let candidate = parent.join(format!("{stem}{extension}"));
        if candidate.is_file() {
            return Some(candidate);
        }
    }
    None
}

fn read_maybe_gzip(path: &Path) -> Result<String> {
    let mut file = File::open(path)?;
    let mut magic = [0u8; 2];
    let read = file.read(&mut magic)?;
    let mut restored = Vec::from(&magic[..read]);
    file.read_to_end(&mut restored)?;
    if read >= 2 && magic[0] == 0x1f && magic[1] == 0x8b {
        let mut decoder = GzDecoder::new(restored.as_slice());
        let mut text = String::new();
        decoder.read_to_string(&mut text)?;
        Ok(text)
    } else {
        // Same text as `from_utf8_lossy(..).into_owned()`, without a copy
        // when the file is already valid UTF-8.
        Ok(match String::from_utf8(restored) {
            Ok(text) => text,
            Err(error) => String::from_utf8_lossy(error.as_bytes()).into_owned(),
        })
    }
}

pub fn is_gff(path: &Path) -> bool {
    let name = path
        .file_name()
        .map(|value| value.to_string_lossy().to_ascii_lowercase())
        .unwrap_or_default();
    name.ends_with(".gff")
        || name.ends_with(".gff3")
        || name.ends_with(".gff.gz")
        || name.ends_with(".gff3.gz")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genetic_code::{is_intact, reverse_complement};

    #[test]
    fn parses_embedded_fasta_and_reverse_strand() {
        let forward = crate::genetic_code::back_translate(b"MKTLLYVW");
        let contig = reverse_complement(&forward);
        let gff = format!(
            "##gff-version 3\ncontig1\tBakta\tCDS\t1\t{}\t.\t-\t0\tID=g1;gene=tesA;product=Test%20protein\n##FASTA\n>contig1\n{contig}\n",
            forward.len()
        );
        let sample = parse_gff("sample", &gff, None).unwrap();
        assert_eq!(sample.genes.len(), 1);
        assert_eq!(sample.genes[0].dna, forward);
        assert_eq!(sample.genes[0].symbol, "tesA");
        assert_eq!(sample.genes[0].product, "Test protein");
        assert_eq!(sample.genes[0].annotation_id, "g1");
        assert!(is_intact(&sample.genes[0].dna));
    }

    #[test]
    fn preserves_decoded_multipart_ids_and_missing_ids() {
        let gff = "##gff-version 3\n\
ctg\tBakta\tCDS\t1\t6\t.\t+\t0\tID=cds%3A1\n\
ctg\tBakta\tCDS\t7\t12\t.\t+\t0\tID=cds%3A1\n\
ctg\tBakta\tCDS\t13\t24\t.\t+\t0\tgene=unnamed\n\
##FASTA\n>ctg\nATGAAAAAATAAATGAAAAAATAA\n";
        let parsed = parse_gff("one", gff, None).unwrap();
        let multipart = parsed
            .genes
            .iter()
            .find(|gene| gene.annotation_id == "cds:1")
            .unwrap();
        assert_eq!((multipart.start, multipart.end), (1, 12));
        assert_eq!(multipart.dna, "ATGAAAAAATAA");
        assert_eq!(
            parsed
                .genes
                .iter()
                .filter(|gene| gene.annotation_id.is_empty())
                .count(),
            1
        );
    }

    fn same(name: &str, text: &str, sibling: Option<&str>) {
        let fast = format!("{:?}", parse_gff(name, text, sibling));
        let reference = format!("{:?}", parse_gff_reference(name, text, sibling));
        assert_eq!(fast, reference, "input:\n{text}");
    }

    #[test]
    fn perf_lean_parser_matches_reference_on_edge_cases() {
        let fasta = "##FASTA\n>ctg extra words\nATGAAACCCGGGTTTAAATAGATGCCCTGA\n>c2\nacgtuNNRYKMatgaaatga\n  ATGTTTTAA \r\n>c2\nATGGGGTAA\n>\nATG\n>c3\nATGé€AAATAA\n";
        let lines = [
            "ctg\tB\tCDS\t1\t30\t.\t+\t0\tID=a;ID=b;gene=x;gene=;Name=n1;Name=n2;product=p%20q",
            "ctg\tB\tCDS\t1\t30\t.\t-\t1\tID= \"quoted\" ;gene=\"g\";product=%zz%4",
            "ctg\tB\tCDS\t+1\t+30\t.\t-\t2\tID=%41%42%4;product=%E2%82%AC;Name=%C3",
            "ctg\tB\tCDS\t4\t29\t.\t+\t3\tID=m;gene=",
            "ctg\tB\tCDS\t4\t29\t.\t+\t.\tID=m2\textra\tcolumns",
            "ctg\tB\tCDS\t4\t29\t.\t+\tx\tID=m3;=v;k;;  ;\u{2003}gene\u{2003}=\u{2003}uni\u{2003}",
            "ctg\tB\tCDS\t0\t29\t.\t+\t0\tID=zero",
            "ctg\tB\tCDS\t20\t10\t.\t+\t0\tID=back",
            "ctg\tB\tCDS\t1\t300\t.\t+\t0\tID=long",
            "ctg\tB\tgene\t1\t30\t.\t+\t0\tID=gene",
            "ctg\tB\tCDS\t1\t30\t.\t+\t0",
            "nope\tB\tCDS\t1\t30\t.\t+\t0\tID=nocontig",
            "c2\tB\tCDS\t1\t9\t.\t+\t0\tID=multi;gene=a",
            "c2\tB\tCDS\t10\t20\t.\t+\t2\tID=multi;product=second",
            "c2\tB\tCDS\t3\t15\t.\t-\t1\tID=minus",
            "c2\tB\tCDS\t16\t21\t.\t-\t0\tID=minus;gene=late",
            "c2\tB\tCDS\t16\t21\t.\t+\t0\tID=minus",
            "c2\tB\tCDS\t1\t21\t.\t?\t0\tgene=noid",
            "c2\tB\tCDS\t2\t21\t.\t.\t0\tgene=noid",
            "c3\tB\tCDS\t1\t12\t.\t+\t0\tID=utf",
            "c3\tB\tCDS\t1\t5\t.\t+\t0\tID=tiny",
            "#ctg\tB\tCDS\t1\t30\t.\t+\t0\tID=comment",
            "",
        ];
        let mut body = String::from("##gff-version 3\n");
        for line in lines {
            body.push_str(line);
            body.push('\n');
        }
        same("s", &format!("{body}{fasta}"), None);
        same("s", &format!("{}{fasta}", body.replace('\n', "\r\n")), None);
        same("s", &body, Some(fasta.trim_start_matches("##FASTA\n")));
        same("s", &format!("{body}##fasta\n>ctg\nATGAAATAA\n"), None);
        same("s", "##FASTA\n>ctg\nATGAAACCCTAA\n", None);
        same("s", &body, None);
        same("s", &format!("{body}##FASTA\n"), None);
        same("s", "ctg\tB\tCDS\t1\t3\t.\t+\t0\tID=a\n##FASTA\n>ctg\nATG\n", None);
    }

    #[test]
    fn perf_lean_parser_matches_reference_on_random_inputs() {
        let mut state = 0x853c49e6748fea9bu64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let keys = ["ID", "gene", "Name", "product", "locus_tag", " ID ", "Gene"];
        let values = ["a", "", "\"q\"", "x%20y", "%", "%4", "%zz", "é", "b c", "dup"];
        let bases = b"ACGTNacgtnURYX";
        for _ in 0..400 {
            let mut fasta = String::from("##FASTA\n");
            let contigs = 1 + next() % 3;
            for c in 0..contigs {
                fasta.push_str(&format!(">c{c} desc\n"));
                let len = 20 + next() % 200;
                for i in 0..len {
                    fasta.push(bases[(next() % bases.len() as u64) as usize] as char);
                    if i % 60 == 59 {
                        fasta.push('\n');
                    }
                }
                fasta.push('\n');
            }
            let mut body = String::new();
            for _ in 0..(next() % 12) {
                let start = next() % 230;
                let end = start + next() % 120;
                let strand = ["+", "-", ".", "?"][(next() % 4) as usize];
                let phase = ["0", "1", "2", "3", ".", "x"][(next() % 6) as usize];
                let mut attrs = String::new();
                for _ in 0..(next() % 6) {
                    attrs.push_str(keys[(next() % keys.len() as u64) as usize]);
                    attrs.push('=');
                    attrs.push_str(values[(next() % values.len() as u64) as usize]);
                    attrs.push(';');
                }
                let contig = next() % (contigs + 1);
                body.push_str(&format!(
                    "c{contig}\tB\tCDS\t{start}\t{end}\t.\t{strand}\t{phase}\t{attrs}\n"
                ));
            }
            same("r", &format!("{body}{fasta}"), None);
        }
    }

    /// Every GFF in a real collection parses to the same `ParsedSample`.
    /// `cargo test --release -- --ignored perf_lean_parser_matches_reference_on_corpus`
    /// with `RUSTROO_GFF_CORPUS=<dir>` (default /workspace/strep_test/gffs).
    #[test]
    #[ignore]
    fn perf_lean_parser_matches_reference_on_corpus() {
        let dir = std::env::var("RUSTROO_GFF_CORPUS")
            .unwrap_or_else(|_| "/workspace/strep_test/gffs".to_string());
        let mut paths: Vec<_> = std::fs::read_dir(&dir)
            .expect("corpus directory")
            .map(|entry| entry.unwrap().path())
            .filter(|path| is_gff(path))
            .collect();
        paths.sort();
        assert!(!paths.is_empty());
        for path in &paths {
            let text = read_maybe_gzip(path).unwrap();
            let sibling = sibling_fasta(path).and_then(|fasta| read_maybe_gzip(&fasta).ok());
            let name = sample_name(path);
            let fast = parse_gff(&name, &text, sibling.as_deref()).unwrap();
            let reference = parse_gff_reference(&name, &text, sibling.as_deref()).unwrap();
            assert_eq!(format!("{fast:?}"), format!("{reference:?}"), "{}", path.display());
        }
        eprintln!("corpus files compared: {}", paths.len());
    }
}
