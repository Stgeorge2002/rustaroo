//! Contig bases held for the gap search, 2 bits per base.
//!
//! `uppercase_dna` only ever emits A, C, G, T and N, so a parsed contig packs
//! exactly: ACGT as 2-bit codes and N as a sorted list of runs. Anything else
//! (only possible for contigs built directly, e.g. in tests) is kept as the
//! original `String`, so every read returns exactly the original bytes.

pub struct PackedSeq {
    len: usize,
    repr: Repr,
}

enum Repr {
    Packed { words: Vec<u64>, n_runs: Vec<(u32, u32)> },
    Raw(String),
}

const CODE: [u8; 256] = {
    let mut table = [255u8; 256];
    table[b'A' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b'N' as usize] = 4;
    table
};
const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

impl PackedSeq {
    pub fn new(sequence: String) -> Self {
        let bytes = sequence.as_bytes();
        let len = bytes.len();
        if len > u32::MAX as usize || bytes.iter().any(|&b| CODE[b as usize] == 255) {
            return Self { len, repr: Repr::Raw(sequence) };
        }
        let mut words = vec![0u64; len.div_ceil(32)];
        let mut n_runs: Vec<(u32, u32)> = Vec::new();
        for (i, &b) in bytes.iter().enumerate() {
            let code = CODE[b as usize];
            if code == 4 {
                match n_runs.last_mut() {
                    Some(run) if run.1 as usize == i => run.1 += 1,
                    _ => n_runs.push((i as u32, i as u32 + 1)),
                }
            } else {
                words[i / 32] |= u64::from(code) << ((i % 32) * 2);
            }
        }
        n_runs.shrink_to_fit();
        Self { len, repr: Repr::Packed { words, n_runs } }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    /// Exactly `original[from..to].to_string()` (same panics on a bad range).
    pub fn slice(&self, from: usize, to: usize) -> String {
        match &self.repr {
            Repr::Raw(text) => text[from..to].to_string(),
            Repr::Packed { words, n_runs } => {
                assert!(from <= to && to <= self.len, "slice {from}..{to} out of range {}", self.len);
                let mut out: Vec<u8> = (from..to)
                    .map(|i| BASES[((words[i / 32] >> ((i % 32) * 2)) & 3) as usize])
                    .collect();
                let first = n_runs.partition_point(|run| (run.1 as usize) <= from);
                for &(start, end) in &n_runs[first..] {
                    let (start, end) = (start as usize, end as usize);
                    if start >= to {
                        break;
                    }
                    for i in start.max(from)..end.min(to) {
                        out[i - from] = b'N';
                    }
                }
                String::from_utf8(out).expect("decoded bases are ASCII")
            }
        }
    }

    #[cfg(test)]
    pub fn unpack(&self) -> String {
        self.slice(0, self.len)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn packed_round_trip_and_every_slice() {
        let mut state = 0x9e3779b97f4a7c15u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for trial in 0..300 {
            let len = (next() % 200) as usize + trial % 3;
            let alphabet: &[char] = if trial % 7 == 0 {
                &['A', 'C', 'G', 'T', 'N', 'a', 'X', 'é']
            } else {
                &['A', 'C', 'G', 'T', 'N', 'N', 'N']
            };
            let mut text = String::new();
            for _ in 0..len {
                text.push(alphabet[(next() % alphabet.len() as u64) as usize]);
            }
            let packed = PackedSeq::new(text.clone());
            assert_eq!(packed.len(), text.len());
            if text.bytes().all(|b| b"ACGTN".contains(&b)) {
                assert!(matches!(packed.repr, Repr::Packed { .. }));
            }
            if text.is_ascii() {
                assert_eq!(packed.unpack(), text);
                for _ in 0..50 {
                    let a = (next() % (text.len() as u64 + 1)) as usize;
                    let b = (next() % (text.len() as u64 + 1)) as usize;
                    let (a, b) = (a.min(b), a.max(b));
                    assert_eq!(packed.slice(a, b), text[a..b]);
                }
            }
        }
        let raw = PackedSeq::new("ACGTé".to_string());
        assert_eq!(raw.unpack(), "ACGTé");
    }
}
