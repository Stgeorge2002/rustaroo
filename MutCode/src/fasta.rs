use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

use crate::io::{digest, identifier, is_dna, reader};
use crate::{fail, Result};

#[derive(Debug, Clone)]
pub struct Sequence {
    pub id: String,
    pub dna: Vec<u8>,
    pub sha256: String,
}

pub struct Reference {
    pub sequences: Vec<Sequence>,
    pub index: HashMap<String, usize>,
}

impl Reference {
    pub fn load(path: &Path) -> Result<Self> {
        let mut sequences = Vec::new();
        let mut index = HashMap::new();
        let mut id = None;
        let mut dna = Vec::new();
        let mut line = String::new();
        let mut input = reader(path)?;
        loop {
            line.clear();
            if input.read_line(&mut line)? == 0 {
                break;
            }
            let text = line.trim();
            if let Some(header) = text.strip_prefix('>') {
                if let Some(previous) = id.take() {
                    push(
                        &mut sequences,
                        &mut index,
                        previous,
                        std::mem::take(&mut dna),
                    )?;
                }
                let name = header
                    .split_whitespace()
                    .next()
                    .ok_or("empty FASTA header")?;
                identifier(name)?;
                id = Some(name.to_string());
            } else if !text.is_empty() {
                if id.is_none() {
                    return fail("sequence before first FASTA header");
                }
                if !is_dna(text.as_bytes()) {
                    return fail(format!(
                        "non-IUPAC DNA in FASTA record {}",
                        id.as_deref().unwrap()
                    ));
                }
                dna.extend_from_slice(text.as_bytes());
            }
        }
        if let Some(previous) = id {
            push(&mut sequences, &mut index, previous, dna)?;
        }
        if sequences.is_empty() {
            return fail("empty FASTA");
        }
        Ok(Self { sequences, index })
    }
}

fn push(
    sequences: &mut Vec<Sequence>,
    index: &mut HashMap<String, usize>,
    id: String,
    dna: Vec<u8>,
) -> Result<()> {
    if dna.is_empty() || index.contains_key(&id) {
        return fail(format!("empty or duplicate FASTA record: {id}"));
    }
    index.insert(id.clone(), sequences.len());
    sequences.push(Sequence {
        id,
        sha256: digest(&dna),
        dna,
    });
    Ok(())
}

pub fn write_record(writer: &mut dyn Write, id: &str, dna: &[u8]) -> Result<()> {
    writeln!(writer, ">{id}")?;
    for chunk in dna.chunks(80) {
        writer.write_all(chunk)?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}
