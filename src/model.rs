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
    pub fusion_piece: bool,
    pub symbol: String,
    pub product: String,
}

#[derive(Debug, Default)]
pub struct Pools {
    pub dna: Vec<String>,
    pub protein: Vec<Vec<u8>>,
    dna_of: std::collections::HashMap<String, u32>,
    protein_of: std::collections::HashMap<Vec<u8>, u32>,
}

impl Pools {
    pub fn intern_dna(&mut self, sequence: String) -> u32 {
        if let Some(&id) = self.dna_of.get(&sequence) {
            return id;
        }
        let id = self.dna.len() as u32;
        self.dna_of.insert(sequence.clone(), id);
        self.dna.push(sequence);
        id
    }

    pub fn intern_protein(&mut self, sequence: Vec<u8>) -> u32 {
        if let Some(&id) = self.protein_of.get(&sequence) {
            return id;
        }
        let id = self.protein.len() as u32;
        self.protein_of.insert(sequence.clone(), id);
        self.protein.push(sequence);
        id
    }
}
