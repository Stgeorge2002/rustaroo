/// A CDS after it has been given integer ids. DNA and protein live in pools so
/// a core gene carried by every sample is stored once.
#[derive(Debug, Clone)]
pub struct Call {
    pub sample: u32,
    pub contig: u32,
    pub start: u32,
    pub end: u32,
    pub strand: i8,
    pub dna: u32,
    pub protein: u32,
    pub cluster: u32,
    /// Homology context anchor: splits do not refine it; proven unions may coarsen it.
    pub context: u32,
    /// Previously supported physical bundle, independent of the current representative.
    pub locus: Option<u32>,
    pub source: Option<u32>,
    /// Found in a gap. Annotated calls stay false even when no audit id is stored.
    pub refound: bool,
    pub fusion_piece: bool,
    pub symbol: String,
    pub product: String,
}

#[derive(Debug, Default)]
pub struct Pools {
    pub dna: Vec<String>,
    pub protein: Vec<Vec<u8>>,
    /// FNV-1a of the stored sequence. A bucket holds every id with that hash;
    /// equality, not the hash, decides a hit. The sequence itself is stored
    /// only in `dna` / `protein`.
    dna_index: std::collections::HashMap<u64, Vec<u32>>,
    protein_index: std::collections::HashMap<u64, Vec<u32>>,
}

impl Pools {
    pub fn intern_dna(&mut self, sequence: String) -> u32 {
        let hash = content_hash(sequence.as_bytes());
        intern_at(&mut self.dna, &mut self.dna_index, hash, sequence)
    }

    pub fn intern_protein(&mut self, sequence: Vec<u8>) -> u32 {
        let hash = content_hash(&sequence);
        intern_at(&mut self.protein, &mut self.protein_index, hash, sequence)
    }
}

pub(crate) fn content_hash(bytes: &[u8]) -> u64 {
    let mut hash = 0xcbf29ce484222325u64;
    for &byte in bytes {
        hash ^= u64::from(byte);
        hash = hash.wrapping_mul(0x100000001b3);
    }
    hash
}

fn intern_at<T: PartialEq>(
    items: &mut Vec<T>,
    index: &mut std::collections::HashMap<u64, Vec<u32>>,
    hash: u64,
    sequence: T,
) -> u32 {
    if let Some(bucket) = index.get(&hash) {
        for &id in bucket {
            if items[id as usize] == sequence {
                return id;
            }
        }
    }
    let id = u32::try_from(items.len()).expect("sequence pool exceeds u32::MAX");
    items.push(sequence);
    index.entry(hash).or_default().push(id);
    id
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_sequences_reuse_one_slot() {
        let mut pools = Pools::default();
        let first = pools.intern_dna("ATGAAATAA".to_string());
        let again = pools.intern_dna("ATGAAATAA".to_string());
        let other = pools.intern_dna("ATGCCCTAA".to_string());
        assert_eq!(first, again);
        assert_ne!(first, other);
        assert_eq!(pools.dna.len(), 2);
        let protein = pools.intern_protein(b"MKK".to_vec());
        assert_eq!(protein, pools.intern_protein(b"MKK".to_vec()));
        assert_eq!(pools.protein.len(), 1);
    }
}
