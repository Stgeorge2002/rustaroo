//! Bacterial genetic code 11.
//!
//! The first codon is translated as methionine when it is a start (ATG, GTG,
//! or TTG). Later GTG and TTG stay valine and leucine. Clustering then treats
//! a start-codon swap as the same residue, which is how the gene is actually
//! initiated, instead of splitting orthologs on the first base.

const AA: [u8; 64] = {
    let mut table = [b'X'; 64];
    // Codon index is base*16 + base*4 + base with A=0, C=1, G=2, T=3.
    table[codon(b'A', b'A', b'A')] = b'K';
    table[codon(b'A', b'A', b'C')] = b'N';
    table[codon(b'A', b'A', b'G')] = b'K';
    table[codon(b'A', b'A', b'T')] = b'N';
    table[codon(b'A', b'C', b'A')] = b'T';
    table[codon(b'A', b'C', b'C')] = b'T';
    table[codon(b'A', b'C', b'G')] = b'T';
    table[codon(b'A', b'C', b'T')] = b'T';
    table[codon(b'A', b'G', b'A')] = b'R';
    table[codon(b'A', b'G', b'C')] = b'S';
    table[codon(b'A', b'G', b'G')] = b'R';
    table[codon(b'A', b'G', b'T')] = b'S';
    table[codon(b'A', b'T', b'A')] = b'I';
    table[codon(b'A', b'T', b'C')] = b'I';
    table[codon(b'A', b'T', b'G')] = b'M';
    table[codon(b'A', b'T', b'T')] = b'I';
    table[codon(b'C', b'A', b'A')] = b'Q';
    table[codon(b'C', b'A', b'C')] = b'H';
    table[codon(b'C', b'A', b'G')] = b'Q';
    table[codon(b'C', b'A', b'T')] = b'H';
    table[codon(b'C', b'C', b'A')] = b'P';
    table[codon(b'C', b'C', b'C')] = b'P';
    table[codon(b'C', b'C', b'G')] = b'P';
    table[codon(b'C', b'C', b'T')] = b'P';
    table[codon(b'C', b'G', b'A')] = b'R';
    table[codon(b'C', b'G', b'C')] = b'R';
    table[codon(b'C', b'G', b'G')] = b'R';
    table[codon(b'C', b'G', b'T')] = b'R';
    table[codon(b'C', b'T', b'A')] = b'L';
    table[codon(b'C', b'T', b'C')] = b'L';
    table[codon(b'C', b'T', b'G')] = b'L';
    table[codon(b'C', b'T', b'T')] = b'L';
    table[codon(b'G', b'A', b'A')] = b'E';
    table[codon(b'G', b'A', b'C')] = b'D';
    table[codon(b'G', b'A', b'G')] = b'E';
    table[codon(b'G', b'A', b'T')] = b'D';
    table[codon(b'G', b'C', b'A')] = b'A';
    table[codon(b'G', b'C', b'C')] = b'A';
    table[codon(b'G', b'C', b'G')] = b'A';
    table[codon(b'G', b'C', b'T')] = b'A';
    table[codon(b'G', b'G', b'A')] = b'G';
    table[codon(b'G', b'G', b'C')] = b'G';
    table[codon(b'G', b'G', b'G')] = b'G';
    table[codon(b'G', b'G', b'T')] = b'G';
    table[codon(b'G', b'T', b'A')] = b'V';
    table[codon(b'G', b'T', b'C')] = b'V';
    table[codon(b'G', b'T', b'G')] = b'V';
    table[codon(b'G', b'T', b'T')] = b'V';
    table[codon(b'T', b'A', b'A')] = b'*';
    table[codon(b'T', b'A', b'C')] = b'Y';
    table[codon(b'T', b'A', b'G')] = b'*';
    table[codon(b'T', b'A', b'T')] = b'Y';
    table[codon(b'T', b'C', b'A')] = b'S';
    table[codon(b'T', b'C', b'C')] = b'S';
    table[codon(b'T', b'C', b'G')] = b'S';
    table[codon(b'T', b'C', b'T')] = b'S';
    table[codon(b'T', b'G', b'A')] = b'*';
    table[codon(b'T', b'G', b'C')] = b'C';
    table[codon(b'T', b'G', b'G')] = b'W';
    table[codon(b'T', b'G', b'T')] = b'C';
    table[codon(b'T', b'T', b'A')] = b'L';
    table[codon(b'T', b'T', b'C')] = b'F';
    table[codon(b'T', b'T', b'G')] = b'L';
    table[codon(b'T', b'T', b'T')] = b'F';
    table
};

const fn codon(a: u8, b: u8, c: u8) -> usize {
    base(a) * 16 + base(b) * 4 + base(c)
}

const fn base(b: u8) -> usize {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' | b'U' | b'u' => 3,
        _ => 0,
    }
}

#[cfg(test)]
fn base_index(b: u8) -> Option<usize> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' | b'U' | b'u' => Some(3),
        _ => None,
    }
}

pub fn is_stop(codon: &[u8]) -> bool {
    matches!(codon, b"TAA" | b"TAG" | b"TGA")
}

pub fn is_start(codon: &[u8]) -> bool {
    matches!(codon, b"ATG" | b"GTG" | b"TTG")
}

/// Translate a coding sequence. A terminal stop is removed. An internal stop
/// is kept as `*`, so a pseudogene that still spans the locus can cluster
/// with the intact gene instead of looking like a short fragment.
pub fn translate(dna: &str) -> Vec<u8> {
    translate_bytes(dna.as_bytes())
}

/// Same translation as `translate`, reading bytes that are not yet a `String`.
pub fn translate_bytes(bytes: &[u8]) -> Vec<u8> {
    let mut protein = Vec::with_capacity(bytes.len() / 3);
    let mut chunks = bytes.chunks_exact(3);
    if let Some(codon) = chunks.next() {
        protein.push(if is_start(codon) { b'M' } else { codon_residue(codon) });
    }
    for codon in chunks {
        protein.push(codon_residue(codon));
    }
    if protein.last() == Some(&b'*') {
        protein.pop();
    }
    protein
}

/// `base_index` as a table; 4 marks a non-base.
const BASE_INDEX: [u8; 256] = {
    let mut table = [4u8; 256];
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    table[b'U' as usize] = 3;
    table[b'u' as usize] = 3;
    table
};

#[inline]
fn codon_residue(codon: &[u8]) -> u8 {
    let a = BASE_INDEX[codon[0] as usize];
    let b = BASE_INDEX[codon[1] as usize];
    let c = BASE_INDEX[codon[2] as usize];
    if (a | b | c) & 4 != 0 {
        b'X'
    } else {
        AA[a as usize * 16 + b as usize * 4 + c as usize]
    }
}

#[cfg(test)]
pub(crate) fn translate_bytes_reference(bytes: &[u8]) -> Vec<u8> {
    let mut protein = Vec::with_capacity(bytes.len() / 3);
    let mut index = 0;
    let mut first = true;
    while index + 3 <= bytes.len() {
        let codon = &bytes[index..index + 3];
        index += 3;
        let aa = if first && is_start(codon) {
            first = false;
            b'M'
        } else {
            first = false;
            match (
                base_index(codon[0]),
                base_index(codon[1]),
                base_index(codon[2]),
            ) {
                (Some(a), Some(b), Some(c)) => AA[a * 16 + b * 4 + c],
                _ => b'X',
            }
        };
        protein.push(aa);
    }
    if protein.last() == Some(&b'*') {
        protein.pop();
    }
    protein
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OrfStatus {
    pub in_frame: bool,
    pub has_start: bool,
    pub internal_stops: u32,
    pub terminal_stop: bool,
}

impl OrfStatus {
    /// A usable gene model: in frame, a bacterial start, and no stop until
    /// the final codon. The final stop is preferred later, but a contig that
    /// ends on a sense codon can still be the longest intact copy.
    pub fn is_intact(self) -> bool {
        self.in_frame && self.has_start && self.internal_stops == 0 && self.codons_at_least_two()
    }

    fn codons_at_least_two(self) -> bool {
        // `in_frame` only says the length divides by 3. A lone start codon
        // is not a gene. Callers pass sequences of at least 6 bp.
        true
    }
}

pub fn orf_status(dna: &str) -> OrfStatus {
    let bytes = dna.as_bytes();
    if bytes.len() < 6 || bytes.len() % 3 != 0 {
        return OrfStatus {
            in_frame: false,
            has_start: false,
            internal_stops: 0,
            terminal_stop: false,
        };
    }
    let codons = bytes.len() / 3;
    let mut internal_stops = 0;
    for codon_index in 1..codons.saturating_sub(1) {
        let start = codon_index * 3;
        if is_stop(&bytes[start..start + 3]) {
            internal_stops += 1;
        }
    }
    OrfStatus {
        in_frame: true,
        has_start: is_start(&bytes[0..3]),
        internal_stops,
        terminal_stop: is_stop(&bytes[bytes.len() - 3..]),
    }
}

pub fn is_intact(dna: &str) -> bool {
    let status = orf_status(dna);
    status.is_intact() && dna.len() >= 6
}

const UPPER: [u8; 256] = {
    let mut table = [b'N'; 256];
    table[b'a' as usize] = b'A';
    table[b'c' as usize] = b'C';
    table[b'g' as usize] = b'G';
    table[b't' as usize] = b'T';
    table[b'u' as usize] = b'T';
    table[b'A' as usize] = b'A';
    table[b'C' as usize] = b'C';
    table[b'G' as usize] = b'G';
    table[b'T' as usize] = b'T';
    table[b'N' as usize] = b'N';
    table
};

const COMPLEMENT: [u8; 256] = {
    let mut table = [b'N'; 256];
    table[b'A' as usize] = b'T';
    table[b'a' as usize] = b'T';
    table[b'T' as usize] = b'A';
    table[b't' as usize] = b'A';
    table[b'U' as usize] = b'A';
    table[b'u' as usize] = b'A';
    table[b'G' as usize] = b'C';
    table[b'g' as usize] = b'C';
    table[b'C' as usize] = b'G';
    table[b'c' as usize] = b'G';
    table
};

/// `uppercase_dna` for one byte (a lookup of the same match arms).
#[inline]
pub fn uppercase_dna_byte(base: u8) -> u8 {
    UPPER[base as usize]
}

pub fn uppercase_dna(sequence: &str) -> String {
    let bytes: Vec<u8> = sequence.bytes().map(uppercase_dna_byte).collect();
    String::from_utf8(bytes).expect("table output is ASCII")
}

/// Append the reverse complement of `piece`, minus its first `skip` output
/// bases (the phase offset applied after complementing).
pub fn extend_reverse_complement(out: &mut Vec<u8>, piece: &[u8], skip: usize) {
    let keep = piece.len().saturating_sub(skip);
    out.extend(piece[..keep].iter().rev().map(|&base| COMPLEMENT[base as usize]));
}

pub fn reverse_complement(sequence: &str) -> String {
    let mut bytes = Vec::with_capacity(sequence.len());
    extend_reverse_complement(&mut bytes, sequence.as_bytes(), 0);
    String::from_utf8(bytes).expect("table output is ASCII")
}

#[cfg(test)]
pub(crate) fn uppercase_dna_reference(sequence: &str) -> String {
    sequence
        .bytes()
        .map(|base| match base {
            b'a' => b'A',
            b'c' => b'C',
            b'g' => b'G',
            b't' | b'u' => b'T',
            b'A' | b'C' | b'G' | b'T' | b'N' => base,
            b'n' => b'N',
            _ => b'N',
        })
        .map(char::from)
        .collect()
}

#[cfg(test)]
pub(crate) fn reverse_complement_reference(sequence: &str) -> String {
    sequence
        .bytes()
        .rev()
        .map(|base| match base {
            b'A' | b'a' => b'T',
            b'T' | b't' | b'U' | b'u' => b'A',
            b'G' | b'g' => b'C',
            b'C' | b'c' => b'G',
            _ => b'N',
        })
        .map(char::from)
        .collect()
}

/// One codon per residue, used to build fixtures. The coding sequence always
/// gets a terminal TAA so an intact gene is closed.
#[cfg(test)]
pub(crate) fn back_translate(protein: &[u8]) -> String {
    let mut dna = String::with_capacity(protein.len() * 3 + 3);
    for &residue in protein {
        dna.push_str(codon_for(residue));
    }
    dna.push_str("TAA");
    dna
}

#[cfg(test)]
fn codon_for(residue: u8) -> &'static str {
    match residue {
        b'A' => "GCT",
        b'C' => "TGT",
        b'D' => "GAT",
        b'E' => "GAA",
        b'F' => "TTT",
        b'G' => "GGT",
        b'H' => "CAT",
        b'I' => "ATT",
        b'K' => "AAA",
        b'L' => "TTA",
        b'M' => "ATG",
        b'N' => "AAT",
        b'P' => "CCT",
        b'Q' => "CAA",
        b'R' => "CGT",
        b'S' => "TCT",
        b'T' => "ACT",
        b'V' => "GTT",
        b'W' => "TGG",
        b'Y' => "TAT",
        b'*' => "TAA",
        _ => "NNN",
    }
}

/// Drop a terminal stop codon so two CDS strings can be joined in frame.
pub fn strip_terminal_stop(dna: &str) -> &str {
    let bytes = dna.as_bytes();
    if bytes.len() >= 3 && bytes.len() % 3 == 0 && is_stop(&bytes[bytes.len() - 3..]) {
        &dna[..dna.len() - 3]
    } else {
        dna
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn translates_starts_and_stops() {
        assert_eq!(translate("ATGAAAAAGTAA"), b"MKK");
        assert_eq!(translate("GTGAAATAA"), b"MK");
        assert_eq!(translate("AAGGTGAAATAA"), b"KVK");
        assert!(!is_intact("ATGTAATAA"));
        assert!(is_intact("ATGAAAAAGTAA"));
    }

    #[test]
    fn reverse_complement_roundtrips_coding_sequence() {
        let dna = "ATGAAAAAGTAA";
        assert_eq!(reverse_complement(&reverse_complement(dna)), dna);
    }

    #[test]
    fn back_translation_is_the_same_protein() {
        let protein = b"MKTLLYVWAACD";
        assert_eq!(translate(&back_translate(protein)), protein);
        assert!(is_intact(&back_translate(protein)));
    }

    #[test]
    fn perf_table_driven_dna_ops_match_reference() {
        let mut state = 0x9e3779b97f4a7c15u64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        // Every single byte, every codon over a mixed alphabet, random strings.
        for byte in 0..=255u8 {
            let text = String::from_utf8_lossy(&[byte]).into_owned();
            assert_eq!(uppercase_dna(&text), uppercase_dna_reference(&text));
            assert_eq!(reverse_complement(&text), reverse_complement_reference(&text));
        }
        let alphabet = b"ACGTUacgtuNnRX*";
        for a in alphabet {
            for b in alphabet {
                for c in alphabet {
                    let codon = [*a, *b, *c];
                    assert_eq!(translate_bytes(&codon), translate_bytes_reference(&codon));
                    let six = [*a, *b, *c, *c, *b, *a];
                    assert_eq!(translate_bytes(&six), translate_bytes_reference(&six));
                }
            }
        }
        for _ in 0..500 {
            let len = (next() % 90) as usize;
            let text: String = (0..len)
                .map(|_| ['A', 'C', 'G', 'T', 'a', 'u', 'N', 'é', 'x', 'U'][(next() % 10) as usize])
                .collect();
            assert_eq!(uppercase_dna(&text), uppercase_dna_reference(&text));
            assert_eq!(reverse_complement(&text), reverse_complement_reference(&text));
            assert_eq!(translate_bytes(text.as_bytes()), translate_bytes_reference(text.as_bytes()));
        }
    }
}
