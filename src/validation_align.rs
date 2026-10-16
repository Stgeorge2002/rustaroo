use std::cell::RefCell;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Step {
    Both,
    ReferenceOnly,
    AlleleOnly,
}

#[derive(Clone, Debug)]
pub(crate) struct Alignment {
    pub steps: Vec<Step>,
    pub ambiguous: bool,
}

const NEG: i32 = i32::MIN / 4;
const OPEN: i32 = -5;
const EXTEND: i32 = -1;

fn best(values: &[i32]) -> (i32, usize, bool) {
    let score = *values.iter().max().expect("nonempty alignment choices");
    let index = values.iter().position(|&value| value == score).unwrap();
    (score, index, score > NEG && values.iter().filter(|&&value| value == score).count() > 1)
}

/// Exact global affine-gap alignment; unlike clustering, no heuristic fallback.
pub(crate) fn align(reference: &[u8], allele: &[u8], cap: usize) -> Result<Alignment, &'static str> {
    if reference.is_empty() || allele.is_empty() {
        return Err("empty_protein");
    }
    if reference == allele {
        return Ok(Alignment { steps: vec![Step::Both; reference.len()], ambiguous: false });
    }
    let n = reference.len();
    let m = allele.len();
    if n > i32::MAX as usize / 8 || m > i32::MAX as usize / 8 {
        return Err("alignment_score_overflow");
    }
    let cells = n.checked_add(1).and_then(|rows| m.checked_add(1).and_then(|cols| rows.checked_mul(cols)))
        .ok_or("alignment_size_overflow")?;
    if cap != 0 && cells > cap {
        return Err("alignment_work_limit");
    }
    // Score rows and the traceback are reused on this worker. The recurrence,
    // tie bits, and predecessor codes are unchanged; only the allocation is.
    AFFINE.with(|cell| cell.borrow_mut().align(reference, allele, n, m, cells))
}

struct AffineScratch {
    trace: Vec<u16>,
    previous: [Vec<i32>; 3],
    current: [Vec<i32>; 3],
}

impl AffineScratch {
    fn align(
        &mut self, reference: &[u8], allele: &[u8], n: usize, m: usize, cells: usize,
    ) -> Result<Alignment, &'static str> {
        let width = m + 1;
        for row in &mut self.previous {
            if row.len() < width {
                row.resize(width, NEG);
            }
            row[..width].fill(NEG);
        }
        for row in &mut self.current {
            if row.len() < width {
                row.resize(width, NEG);
            }
        }
        if self.trace.len() < cells {
            self.trace.resize(cells, 0);
        }
        self.previous[0][0] = 0;
        for j in 1..=m {
            self.previous[2][j] = OPEN + (j - 1) as i32 * EXTEND;
            self.trace[j] = if j == 1 { 0 } else { 2 << 4 };
        }
        for i in 1..=n {
            for row in &mut self.current {
                row[..width].fill(NEG);
            }
            self.current[1][0] = OPEN + (i - 1) as i32 * EXTEND;
            self.trace[i * (m + 1)] = if i == 1 { 0 } else { 1 << 2 };
            for j in 1..=m {
                let (score, state, tied) = best(&[
                    self.previous[0][j - 1],
                    self.previous[1][j - 1],
                    self.previous[2][j - 1],
                ]);
                let canonical = reference[i - 1] != b'X' && reference[i - 1] != b'*';
                self.current[0][j] = score + if reference[i - 1] == allele[j - 1] && canonical { 2 } else { -1 };
                let mut code = state as u16 | (u16::from(tied) << 6);
                let (score, state, tied) = best(&[self.previous[0][j] + OPEN, self.previous[1][j] + EXTEND]);
                self.current[1][j] = score;
                code |= (state as u16) << 2 | (u16::from(tied) << 7);
                let (score, state, tied) = best(&[self.current[0][j - 1] + OPEN, self.current[2][j - 1] + EXTEND]);
                self.current[2][j] = score;
                code |= (if state == 0 { 0 } else { 2 }) << 4 | (u16::from(tied) << 8);
                self.trace[i * (m + 1) + j] = code;
            }
            std::mem::swap(&mut self.previous, &mut self.current);
        }
        let (_, mut state, mut ambiguous) = best(&[self.previous[0][m], self.previous[1][m], self.previous[2][m]]);
        let (mut i, mut j) = (n, m);
        let mut steps = Vec::with_capacity(n.max(m));
        while i > 0 || j > 0 {
            let code = self.trace[i * (m + 1) + j];
            ambiguous |= code & (1 << (6 + state)) != 0;
            let predecessor = ((code >> (state * 2)) & 3) as usize;
            match state {
                0 => { i -= 1; j -= 1; steps.push(Step::Both); }
                1 => { i -= 1; steps.push(Step::ReferenceOnly); }
                2 => { j -= 1; steps.push(Step::AlleleOnly); }
                _ => unreachable!("alignment state"),
            }
            state = predecessor;
        }
        steps.reverse();
        Ok(Alignment { steps, ambiguous })
    }
}

thread_local! {
    static AFFINE: RefCell<AffineScratch> = RefCell::new(AffineScratch {
        trace: Vec::new(),
        previous: std::array::from_fn(|_| Vec::new()),
        current: std::array::from_fn(|_| Vec::new()),
    });
}

#[derive(Clone)]
pub(crate) struct Profile {
    pub positions: Vec<Option<usize>>,
    pub insertions: Vec<Vec<usize>>,
    pub ambiguous: bool,
}

impl Profile {
    pub fn new(reference_len: usize, alignment: &Alignment) -> Self {
        let mut positions = vec![None; reference_len];
        let mut insertions = vec![Vec::new(); reference_len + 1];
        let (mut r, mut a) = (0, 0);
        for step in &alignment.steps {
            match step {
                Step::Both => { positions[r] = Some(a); r += 1; a += 1; }
                Step::ReferenceOnly => { r += 1; }
                Step::AlleleOnly => { insertions[r].push(a); a += 1; }
            }
        }
        Self { positions, insertions, ambiguous: alignment.ambiguous }
    }

    pub fn codon_core(&self, dna: &[u8]) -> Vec<u8> {
        let mut row = Vec::with_capacity(self.positions.len() * 3);
        for position in &self.positions {
            if let Some(index) = position {
                row.extend_from_slice(&dna[index * 3..index * 3 + 3]);
            } else {
                row.extend_from_slice(b"---");
            }
        }
        row
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_identity_and_work_limit_are_explicit() {
        assert_eq!(align(b"MACD", b"MACD", 1).unwrap().steps, vec![Step::Both; 4]);
        assert_eq!(align(b"MACD", b"MAC", 1).unwrap_err(), "alignment_work_limit");
        assert_eq!(align(b"", b"MAC", 0).unwrap_err(), "empty_protein");
    }

    #[test]
    fn traceback_localizes_terminal_and_internal_losses() {
        let reference = b"MACDEFGHIKLMNPQRSTVWY";
        let suffix = align(reference, &reference[5..], 0).unwrap();
        assert_eq!(&suffix.steps[..5], &[Step::ReferenceOnly; 5]);
        let mut deletion = reference[..7].to_vec();
        deletion.extend_from_slice(&reference[12..]);
        let alignment = align(reference, &deletion, 0).unwrap();
        let profile = Profile::new(reference.len(), &alignment);
        assert!(profile.positions[7..12].iter().all(Option::is_none));
        assert_eq!(profile.positions.iter().flatten().count(), deletion.len());
    }

    #[test]
    fn insertions_are_not_projected_onto_reference_tree_sites() {
        let reference = b"MACDEFGHIKLMNPQRSTVWY";
        let mut allele = reference[..8].to_vec();
        allele.extend_from_slice(b"AAAAA");
        allele.extend_from_slice(&reference[8..]);
        let profile = Profile::new(reference.len(), &align(reference, &allele, 0).unwrap());
        assert_eq!(profile.insertions[8].len(), 5);
        let dna = crate::genetic_code::back_translate(&allele);
        assert_eq!(profile.codon_core(dna.as_bytes()).len(), reference.len() * 3);
    }

    #[test]
    fn repeat_placement_ties_are_reported() {
        assert!(align(b"MAAAAAAAAC", b"MAAAAAC", 0).unwrap().ambiguous);
    }

    #[test]
    fn traceback_scores_match_exhaustive_affine_paths_for_small_pairs() {
        fn exhaustive(a: &[u8], b: &[u8], state: usize) -> i32 {
            if a.is_empty() && b.is_empty() { return 0; }
            let mut score = NEG;
            if !a.is_empty() && !b.is_empty() {
                score = (if a[0] == b[0] { 2 } else { -1 }) + exhaustive(&a[1..], &b[1..], 0);
            }
            if !a.is_empty() && state != 2 {
                score = score.max((if state == 1 { EXTEND } else { OPEN }) + exhaustive(&a[1..], b, 1));
            }
            if !b.is_empty() && state != 1 {
                score = score.max((if state == 2 { EXTEND } else { OPEN }) + exhaustive(a, &b[1..], 2));
            }
            score
        }
        for a in [b"A".as_slice(), b"AC", b"ACA", b"CCAA"] {
            for b in [b"C".as_slice(), b"AA", b"CAC", b"AACC"] {
                let alignment = align(a, b, 0).unwrap();
                let (mut i, mut j, mut state, mut score) = (0, 0, 0, 0);
                for step in alignment.steps {
                    match step {
                        Step::Both => { score += if a[i] == b[j] { 2 } else { -1 }; i += 1; j += 1; state = 0; }
                        Step::ReferenceOnly => { score += if state == 1 { EXTEND } else { OPEN }; i += 1; state = 1; }
                        Step::AlleleOnly => { score += if state == 2 { EXTEND } else { OPEN }; j += 1; state = 2; }
                    }
                }
                assert_eq!((i, j), (a.len(), b.len()));
                assert_eq!(score, exhaustive(a, b, 0));
            }
        }
    }
}
