use std::collections::{BTreeMap, BTreeSet};

use rayon::prelude::*;

use crate::{Error, Result};

#[derive(Clone)]
struct Node {
    children: Vec<(usize, f64)>,
}

pub(crate) struct Tree {
    nodes: Vec<Node>,
    root: usize,
    tips: usize,
}

pub(crate) fn neighbor_joining(distances: &[Vec<f64>]) -> Result<Tree> {
    let n = distances.len();
    if n < 3 || distances.iter().any(|row| row.len() != n) {
        return Err(Error::new(
            "neighbor joining requires a square matrix with at least three tips",
        ));
    }
    for i in 0..n {
        for j in 0..n {
            let d = distances[i][j];
            if !d.is_finite()
                || d < 0.0
                || (d - distances[j][i]).abs() > 1e-9
                || (i == j && d.abs() > 1e-9)
            {
                return Err(Error::new("invalid neighbor-joining distance matrix"));
            }
        }
    }
    let mut matrix = vec![vec![0.0; n * 2]; n * 2];
    for i in 0..n {
        matrix[i][..n].copy_from_slice(&distances[i]);
    }
    let mut nodes = vec![
        Node {
            children: Vec::new()
        };
        n
    ];
    let mut active: Vec<_> = (0..n).collect();
    while active.len() > 3 {
        let m = active.len();
        let sums: BTreeMap<_, _> = active
            .iter()
            .map(|&i| (i, active.iter().map(|&j| matrix[i][j]).sum::<f64>()))
            .collect();
        let mut choice = (f64::INFINITY, active[0], active[1]);
        for (position, &i) in active.iter().enumerate() {
            for &j in &active[position + 1..] {
                let q = (m - 2) as f64 * matrix[i][j] - sums[&i] - sums[&j];
                if q < choice.0 {
                    choice = (q, i, j);
                }
            }
        }
        let (_, i, j) = choice;
        let left = matrix[i][j] / 2.0 + (sums[&i] - sums[&j]) / (2.0 * (m - 2) as f64);
        let right = matrix[i][j] - left;
        let u = nodes.len();
        nodes.push(Node {
            children: vec![(i, left), (j, right)],
        });
        for &k in &active {
            if k != i && k != j {
                let distance = (matrix[i][k] + matrix[j][k] - matrix[i][j]) / 2.0;
                matrix[u][k] = distance;
                matrix[k][u] = distance;
            }
        }
        active.retain(|&k| k != i && k != j);
        active.push(u);
    }
    let (a, b, c) = (active[0], active[1], active[2]);
    let a_length = (matrix[a][b] + matrix[a][c] - matrix[b][c]) / 2.0;
    let root = nodes.len();
    nodes.push(Node {
        children: vec![
            (a, a_length),
            (b, matrix[a][b] - a_length),
            (c, matrix[a][c] - a_length),
        ],
    });
    Ok(Tree {
        nodes,
        root,
        tips: n,
    })
}

impl Tree {
    pub fn negative_branches(&self) -> usize {
        self.nodes
            .iter()
            .flat_map(|node| &node.children)
            .filter(|(_, length)| *length < -1e-10)
            .count()
    }

    pub fn splits(&self) -> BTreeMap<Vec<usize>, f64> {
        fn visit(
            tree: &Tree,
            index: usize,
            splits: &mut BTreeMap<Vec<usize>, f64>,
        ) -> BTreeSet<usize> {
            if index < tree.tips {
                return BTreeSet::from([index]);
            }
            let mut leaves = BTreeSet::new();
            for &(child, length) in &tree.nodes[index].children {
                let subset = visit(tree, child, splits);
                if subset.len() >= 2 && tree.tips - subset.len() >= 2 {
                    let a: Vec<_> = subset.iter().copied().collect();
                    let b: Vec<_> = (0..tree.tips).filter(|tip| !subset.contains(tip)).collect();
                    splits.insert(a.min(b), length);
                }
                leaves.extend(subset);
            }
            leaves
        }
        let mut splits = BTreeMap::new();
        visit(self, self.root, &mut splits);
        splits
    }

    pub fn newick(
        &self,
        labels: &[String],
        support: &BTreeMap<Vec<usize>, usize>,
        replicates: usize,
    ) -> String {
        fn render(
            tree: &Tree,
            index: usize,
            labels: &[String],
            support: &BTreeMap<Vec<usize>, usize>,
            replicates: usize,
        ) -> (String, BTreeSet<usize>) {
            if index < tree.tips {
                return (
                    format!("'{}'", labels[index].replace('\'', "''")),
                    BTreeSet::from([index]),
                );
            }
            let mut parts = Vec::new();
            let mut leaves = BTreeSet::new();
            for &(child, length) in &tree.nodes[index].children {
                let (text, subset) = render(tree, child, labels, support, replicates);
                parts.push(format!("{text}:{length:.8}"));
                leaves.extend(subset);
            }
            let mut text = format!("({})", parts.join(","));
            if index != tree.root
                && leaves.len() >= 2
                && tree.tips - leaves.len() >= 2
                && replicates > 0
            {
                let a: Vec<_> = leaves.iter().copied().collect();
                let b: Vec<_> = (0..tree.tips).filter(|tip| !leaves.contains(tip)).collect();
                let count = support.get(&a.min(b)).copied().unwrap_or(0);
                text.push_str(&format!("{:.1}", count as f64 * 100.0 / replicates as f64));
            }
            (text, leaves)
        }
        format!(
            "{};",
            render(self, self.root, labels, support, replicates).0
        )
    }
}

pub(crate) fn common_sites(rows: &[Vec<u8>]) -> Result<Vec<usize>> {
    let Some(first) = rows.first() else {
        return Err(Error::new("no alignment rows for a tree"));
    };
    if rows.iter().any(|row| row.len() != first.len()) {
        return Err(Error::new("unequal tree alignment row lengths"));
    }
    Ok((0..first.len())
        .filter(|&site| {
            rows.iter()
                .all(|row| matches!(row[site], b'A' | b'C' | b'G' | b'T'))
        })
        .collect())
}

pub(crate) fn variable_sites(rows: &[Vec<u8>], sites: &[usize]) -> usize {
    sites
        .iter()
        .filter(|&&site| rows.iter().skip(1).any(|row| row[site] != rows[0][site]))
        .count()
}

#[derive(Debug, PartialEq, Eq)]
pub(crate) enum DistanceError {
    NoSites,
    Saturated,
}

impl std::fmt::Display for DistanceError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::NoSites => "no_common_nucleotide_sites",
            Self::Saturated => "jc69_saturation",
        })
    }
}

pub(crate) fn distances(
    rows: &[Vec<u8>],
    sites: &[usize],
) -> std::result::Result<Vec<Vec<f64>>, DistanceError> {
    if sites.is_empty() {
        return Err(DistanceError::NoSites);
    }
    let n = rows.len();
    let mut result = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in i + 1..n {
            let mismatches = sites
                .iter()
                .filter(|&&site| rows[i][site] != rows[j][site])
                .count();
            let p = mismatches as f64 / sites.len() as f64;
            if p >= 0.75 {
                return Err(DistanceError::Saturated);
            }
            let d = -0.75 * (1.0 - 4.0 * p / 3.0).ln();
            result[i][j] = d;
            result[j][i] = d;
        }
    }
    Ok(result)
}

pub(crate) fn bootstrap(
    rows: &[Vec<u8>],
    sites: &[usize],
    replicates: usize,
) -> Result<(BTreeMap<Vec<usize>, usize>, usize, usize)> {
    let mut support = BTreeMap::new();
    if replicates == 0 {
        return Ok((support, 0, 0));
    }
    // The xorshift is a chain. Draw every replicate's columns first so the
    // stream matches a serial pass, then build the trees together.
    let mut state = 0x72a9_14bc_5e37_8d01u64;
    let mut sampled_sets = Vec::with_capacity(replicates);
    for _ in 0..replicates {
        let sampled: Vec<_> = (0..sites.len())
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                sites[(state % sites.len() as u64) as usize]
            })
            .collect();
        sampled_sets.push(sampled);
    }
    let outcomes: Vec<Result<Option<Vec<Vec<usize>>>>> = sampled_sets
        .par_iter()
        .map(|sampled| {
            if variable_sites(rows, sampled) == 0 {
                return Ok(None);
            }
            let matrix = match distances(rows, sampled) {
                Ok(matrix) => matrix,
                Err(DistanceError::Saturated) => return Ok(None),
                Err(DistanceError::NoSites) => {
                    return Err(Error::new("bootstrap has no common nucleotide sites"));
                }
            };
            let tree = neighbor_joining(&matrix)?;
            Ok(Some(
                tree.splits()
                    .into_iter()
                    .filter(|(_, length)| *length > 1e-10)
                    .map(|(split, _)| split)
                    .collect(),
            ))
        })
        .collect();
    let (mut usable, mut rejected) = (0, 0);
    for outcome in outcomes {
        match outcome? {
            None => rejected += 1,
            Some(splits) => {
                for split in splits {
                    *support.entry(split).or_default() += 1;
                }
                usable += 1;
            }
        }
    }
    Ok((support, usable, rejected))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn additive_four_tip_tree_recovers_the_known_split() {
        let matrix = vec![
            vec![0., 2., 6., 6.],
            vec![2., 0., 6., 6.],
            vec![6., 6., 0., 2.],
            vec![6., 6., 2., 0.],
        ];
        let tree = neighbor_joining(&matrix).unwrap();
        assert_eq!(
            tree.splits().keys().cloned().collect::<Vec<_>>(),
            vec![vec![0, 1]]
        );
        assert_eq!(tree.negative_branches(), 0);
        assert!(tree
            .newick(
                &["a".into(), "b".into(), "c".into(), "d".into()],
                &BTreeMap::new(),
                0
            )
            .ends_with(';'));
    }

    #[test]
    fn missing_sites_and_saturation_are_not_zero_distance_fallbacks() {
        let rows = vec![b"AC-N".to_vec(), b"ACTN".to_vec(), b"ACGN".to_vec()];
        assert_eq!(common_sites(&rows).unwrap(), vec![0, 1]);
        assert!(distances(&[b"AAAA".to_vec(), b"CCCC".to_vec()], &[0, 1, 2, 3]).is_err());
        assert!(neighbor_joining(&[vec![0., -1.], vec![-1., 0.]]).is_err());
    }

    #[test]
    fn bootstrap_is_deterministic_and_supports_a_known_split() {
        let rows = vec![
            b"AAAAAAAAAACCCCCCCCCC".to_vec(),
            b"AAAAAAAAAACCCCCCCCCC".to_vec(),
            b"AAAAAAAAAAGGGGGGGGGG".to_vec(),
            b"AAAAAAAAAAGGGGGGGGGG".to_vec(),
        ];
        let sites = common_sites(&rows).unwrap();
        let first = bootstrap(&rows, &sites, 20).unwrap();
        assert!(first.1 > 0);
        assert_eq!(first, bootstrap(&rows, &sites, 20).unwrap());
        assert_eq!(first.0.get(&vec![0, 1]), Some(&first.1));
    }

    #[test]
    fn nonadditive_distances_keep_negative_branches_visible() {
        let tree =
            neighbor_joining(&[vec![0., 1., 1.], vec![1., 0., 4.], vec![1., 4., 0.]]).unwrap();
        assert_eq!(tree.negative_branches(), 1);
        assert!(tree
            .newick(&["a".into(), "b".into(), "c".into()], &BTreeMap::new(), 0)
            .contains(":-1.00000000"));
    }

    #[test]
    fn rejected_bootstraps_do_not_inflate_conditional_support() {
        let matrix = vec![
            vec![0., 2., 6., 6.],
            vec![2., 0., 6., 6.],
            vec![6., 6., 0., 2.],
            vec![6., 6., 2., 0.],
        ];
        let tree = neighbor_joining(&matrix).unwrap();
        let support = BTreeMap::from([(vec![0, 1], 20)]);
        assert!(tree
            .newick(
                &["a".into(), "b".into(), "c".into(), "d".into()],
                &support,
                50
            )
            .contains(")40.0:"));
    }
}
