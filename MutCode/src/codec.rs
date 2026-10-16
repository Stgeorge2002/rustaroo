use serde::{Deserialize, Serialize};

use crate::io::is_dna;
use crate::{fail, Result};

/// Reference offset, reference bases deleted, literal inserted bases.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Edit(pub usize, pub usize, pub String);

pub fn replacement(start: usize, reference: &[u8], alternate: &str) -> Option<Edit> {
    let alt = alternate.as_bytes();
    let prefix = reference
        .iter()
        .zip(alt)
        .take_while(|(a, b)| a == b)
        .count();
    let suffix = reference[prefix..]
        .iter()
        .rev()
        .zip(alt[prefix..].iter().rev())
        .take_while(|(a, b)| a == b)
        .count();
    let deleted = reference.len() - prefix - suffix;
    let inserted = &alternate[prefix..alt.len() - suffix];
    if deleted == 0 && inserted.is_empty() {
        None
    } else {
        Some(Edit(start + prefix, deleted, inserted.to_string()))
    }
}

pub fn validate(reference_len: usize, edits: &[Edit]) -> Result<usize> {
    let mut cursor = 0;
    let mut previous_start = None;
    let mut target_len = reference_len;
    for Edit(start, deleted, inserted) in edits {
        let end = start
            .checked_add(*deleted)
            .ok_or("edit coordinate overflow")?;
        if *start < cursor
            || end > reference_len
            || previous_start.is_some_and(|previous| previous >= *start)
        {
            return fail("edits are overlapping, unsorted, duplicated or out of bounds");
        }
        if (*deleted == 0 && inserted.is_empty()) || !is_dna(inserted.as_bytes()) {
            return fail("empty edit or nonliteral/non-IUPAC inserted sequence");
        }
        target_len = target_len
            .checked_sub(*deleted)
            .and_then(|length| length.checked_add(inserted.len()))
            .ok_or("target length overflow")?;
        cursor = end;
        previous_start = Some(*start);
    }
    Ok(target_len)
}

pub fn apply(reference: &[u8], edits: &[Edit]) -> Result<Vec<u8>> {
    let mut sequence = Vec::with_capacity(validate(reference.len(), edits)?);
    let mut cursor = 0;
    for Edit(start, deleted, inserted) in edits {
        sequence.extend_from_slice(&reference[cursor..*start]);
        sequence.extend_from_slice(inserted.as_bytes());
        cursor = start + deleted;
    }
    sequence.extend_from_slice(&reference[cursor..]);
    Ok(sequence)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_literal_mutation_classes_and_boundaries() {
        for (reference, alternate, expected) in [
            ("A", "T", "T"),
            ("AC", "GT", "GT"),
            ("A", "ACT", "ACT"),
            ("ACT", "A", "A"),
            ("ACGT", "AGTTA", "AGTTA"),
            ("ACGT", "NCGT", "NCGT"),
            ("ACGT", "acgt", "acgt"),
            ("AAA", "AAAA", "AAAA"),
            ("ACGT", "ACGT", "ACGT"),
        ] {
            let edits: Vec<_> = replacement(0, reference.as_bytes(), alternate)
                .into_iter()
                .collect();
            assert_eq!(
                apply(reference.as_bytes(), &edits).unwrap(),
                expected.as_bytes()
            );
        }
        let edits = vec![Edit(0, 0, "N".into()), Edit(4, 0, "T".into())];
        assert_eq!(apply(b"ACGT", &edits).unwrap(), b"NACGTT");
        assert_eq!(apply(b"ACGT", &[Edit(0, 4, "".into())]).unwrap(), b"");
    }

    #[test]
    fn reject_invalid_edits() {
        for edits in [
            vec![Edit(3, 2, "".into())],
            vec![Edit(1, 2, "".into()), Edit(2, 0, "A".into())],
            vec![Edit(1, 0, "A".into()), Edit(1, 0, "T".into())],
            vec![Edit(0, 0, "".into())],
            vec![Edit(0, 1, "<DEL>".into())],
        ] {
            assert!(apply(b"ACGT", &edits).is_err());
        }
    }

    #[test]
    fn exhaustive_short_replacements_roundtrip() {
        let mut words = vec![String::new()];
        for length in 1..=4 {
            for bits in 0..1 << length {
                words.push(
                    (0..length)
                        .map(|bit| if bits & (1 << bit) == 0 { 'A' } else { 'N' })
                        .collect(),
                );
            }
        }
        for reference in &words {
            for alternate in &words {
                let edits: Vec<_> = replacement(0, reference.as_bytes(), alternate)
                    .into_iter()
                    .collect();
                assert_eq!(
                    apply(reference.as_bytes(), &edits).unwrap(),
                    alternate.as_bytes()
                );
            }
        }
    }
}
