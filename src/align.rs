//! Pairwise protein identity.
//!
//! Homology uses a global alignment and scores identical residues over the
//! shorter sequence, which is the identity CD-HIT reports. Deciding whether a
//! long allele is a clean extension of a short one uses k-mer diagonals: a
//! terminal truncation sits on one diagonal, an internal insertion splits into
//! two. That test only has to be right for alleles that are already nearly the
//! same gene.

use std::cell::RefCell;

const MATCH: i32 = 2;
const MISMATCH: i32 = -1;
const GAP: i32 = -2;
const KMER: usize = 4;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Fit {
    pub identity: f64,
    pub coverage: f64,
    pub long_start: usize,
    pub long_end: usize,
    pub internal_gap: usize,
}

impl Fit {
    pub fn extra_prefix(self, long_len: usize) -> usize {
        self.long_start.min(long_len)
    }

    pub fn extra_suffix(self, long_len: usize) -> usize {
        long_len.saturating_sub(self.long_end)
    }
}

pub fn length_ratio(a: &[u8], b: &[u8]) -> f64 {
    let (short, long) = if a.len() < b.len() {
        (a.len(), b.len())
    } else {
        (b.len(), a.len())
    };
    if long == 0 {
        0.0
    } else {
        short as f64 / long as f64
    }
}

/// Best full-coverage ungapped identity within two residues of either terminus.
/// Two slices allow testing an observed fragment join without allocating it.
pub(crate) fn terminal_identity_parts(left: &[u8], right: &[u8], long: &[u8]) -> Option<f64> {
    let length = left.len() + right.len();
    if length == 0 || length > long.len() {
        return None;
    }
    if right.is_empty() && left == long {
        return Some(1.0);
    }
    let gap = long.len() - length;
    let prefix_end = gap.min(2);
    let suffix_start = gap.saturating_sub(2).max(prefix_end + 1);
    (0..=prefix_end)
        .chain(suffix_start..=gap)
        .map(|offset| {
            let matches = left
                .iter()
                .chain(right)
                .zip(&long[offset..offset + length])
                .filter(|(a, b)| a == b)
                .count();
            matches as f64 / length as f64
        })
        .max_by(f64::total_cmp)
}

/// Identity of a global alignment, identical residues divided by the shorter
/// length. Equal sequences never enter the matrix. Everything else is
/// computed in a band around the center diagonal, so the cost tracks how far
/// the two sequences actually drift apart in length rather than their full
/// product; `None` only when that drift is too large for a band to be exact
/// and affordable, in which case the k-mer diagonal fitter below stands in.
pub fn global_identity(a: &[u8], b: &[u8]) -> Option<f64> {
    if a.is_empty() || b.is_empty() {
        return None;
    }
    if a == b {
        return Some(1.0);
    }
    let (short, long) = if a.len() <= b.len() { (a, b) } else { (b, a) };
    let band = band_width(short.len(), long.len());
    let work = short.len().saturating_mul(band.saturating_mul(2) + 1);
    if work > BANDED_WORK_CAP {
        let fit = fit_short_in_long(short, long)?;
        return if fit.coverage < 0.8 {
            None
        } else {
            Some(fit.identity)
        };
    }
    let matches = banded_matches(short, long, band)?;
    Some(matches as f64 / short.len() as f64)
}

const BANDED_WORK_CAP: usize = 20_000_000;

/// PERF PROTOTYPE: same value as `global_identity` whenever that value could
/// satisfy `id + f64::EPSILON >= min_identity`; otherwise `None`. Callers that
/// discard below-threshold identities (cluster `choose_representative`) see
/// identical decisions. Uses (1) an exact residue-composition upper bound on
/// the number of identical aligned pairs and (2) an exact per-row upper bound
/// inside the banded DP (best matches so far in the row + remaining rows).
use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};
static CALLS: AtomicU64 = AtomicU64::new(0);
static COMP_PRUNED: AtomicU64 = AtomicU64::new(0);
static DP_RUN: AtomicU64 = AtomicU64::new(0);
static DP_ABORT: AtomicU64 = AtomicU64::new(0);
static CELLS_FULL: AtomicU64 = AtomicU64::new(0);
static CELLS_DONE: AtomicU64 = AtomicU64::new(0);
pub fn report_counters(label: &str) {
    if !crate::log::stages_enabled() {
        return;
    }
    eprintln!("COUNTERS\t{label}\tcalls={}\tcomp_pruned={}\tdp_run={}\tdp_abort={}\tcells_full={}\tcells_done={}",
        CALLS.swap(0, AtomicOrdering::Relaxed), COMP_PRUNED.swap(0, AtomicOrdering::Relaxed), DP_RUN.swap(0, AtomicOrdering::Relaxed),
        DP_ABORT.swap(0, AtomicOrdering::Relaxed), CELLS_FULL.swap(0, AtomicOrdering::Relaxed), CELLS_DONE.swap(0, AtomicOrdering::Relaxed));
}
pub fn global_identity_at_least(a: &[u8], b: &[u8], min_identity: f64) -> Option<f64> {
    #[cfg(feature = "counters")]
    CALLS.fetch_add(1, AtomicOrdering::Relaxed);
    if a.is_empty() || b.is_empty() {
        return None;
    }
    if a == b {
        return Some(1.0);
    }
    let (short, long) = if a.len() <= b.len() { (a, b) } else { (b, a) };
    let band = band_width(short.len(), long.len());
    let work = short.len().saturating_mul(band.saturating_mul(2) + 1);
    if work > BANDED_WORK_CAP {
        // Unbounded fallback path is kept exactly as before.
        return global_identity(a, b);
    }
    let n = short.len();
    let passes = |matches: usize| matches as f64 / n as f64 + f64::EPSILON >= min_identity;
    // Smallest match count that could still pass; predicate is monotone in k.
    let mut needed = n + 1;
    {
        let (mut lo, mut hi) = (0usize, n + 1);
        while lo < hi {
            let mid = (lo + hi) / 2;
            if mid <= n && passes(mid) { hi = mid; } else { lo = mid + 1; }
        }
        if lo <= n { needed = lo; }
    }
    if needed > n {
        return None;
    }
    if needed > 0 {
        let mut counts = [0u32; 256];
        for &r in short { counts[r as usize] += 1; }
        let mut other = [0u32; 256];
        for &r in long { other[r as usize] += 1; }
        let mut bound = 0usize;
        for k in 0..256 { bound += counts[k].min(other[k]) as usize; }
        if bound < needed {
            #[cfg(feature = "counters")]
            COMP_PRUNED.fetch_add(1, AtomicOrdering::Relaxed);
            return None;
        }
    }
    #[cfg(feature = "counters")]
    DP_RUN.fetch_add(1, AtomicOrdering::Relaxed);
    #[cfg(feature = "counters")]
    CELLS_FULL.fetch_add(work as u64, AtomicOrdering::Relaxed);
    let matches = BAND.with(|cell| cell.borrow_mut().matches_bounded(short, long, band, needed))?;
    Some(matches as f64 / short.len() as f64)
}

/// Half-width of the band around the center diagonal: the length gap itself
/// plus 10% slack either side for indels. Candidates reaching here already
/// passed the caller's length-ratio filter, so the true optimal alignment
/// stays inside this band; it is the same width the old large-sequence-only
/// fallback used, just applied to every pair instead of a rare few.
fn band_width(short_len: usize, long_len: usize) -> usize {
    short_len.abs_diff(long_len) + long_len / 10 + 16
}

/// Exact global-alignment match count, restricted to a band around the
/// center diagonal. Two rolling rows of width `O(band)` stand in for the old
/// full `(n+1) x (m+1)` score and traceback matrices: the match count that
/// used to come from a backward traceback is instead carried forward
/// alongside the score at each cell, so neither a full matrix nor a
/// traceback pass is needed. When `band` covers the whole matrix this is
/// exactly the old full computation.
///
/// The rows live in a thread-local buffer. A fresh `Vec` per pair dominated
/// the runtime once the band itself got small, and each worker only ever
/// aligns one pair at a time.
fn banded_matches(short: &[u8], long: &[u8], band: usize) -> Option<usize> {
    BAND.with(|cell| cell.borrow_mut().matches(short, long, band))
}

const NEG: i32 = i32::MIN / 4;

struct BandScratch {
    prev_score: Vec<i32>,
    prev_matches: Vec<u32>,
    curr_score: Vec<i32>,
    curr_matches: Vec<u32>,
}

impl BandScratch {
    fn matches(&mut self, short: &[u8], long: &[u8], band: usize) -> Option<usize> {
        let n = short.len();
        let m = long.len();
        let width = band.saturating_mul(2) + 1;
        self.ensure(width);
        // The buffer outlives this pair. A gap row must start from zero
        // matches; leftover counts from the previous pair are not a diagonal.
        self.prev_score[..width].fill(NEG);
        self.prev_matches[..width].fill(0);
        self.curr_score[..width].fill(NEG);
        self.curr_matches[..width].fill(0);

        let mut prev_lo = 0usize;
        let mut prev_hi = band.min(m);
        for j in prev_lo..=prev_hi {
            self.prev_score[j - prev_lo] = j as i32 * GAP;
        }

        for i in 1..=n {
            let center = (i * m) / n;
            let lo = center.saturating_sub(band);
            let hi = (center + band).min(m);
            for j in lo..=hi {
                let diag = if j == 0 {
                    (NEG, 0)
                } else {
                    let (ps, pm) = band_cell(
                        prev_lo,
                        prev_hi,
                        &self.prev_score,
                        &self.prev_matches,
                        j - 1,
                    );
                    if ps <= NEG {
                        (NEG, 0)
                    } else {
                        let is_match = short[i - 1] == long[j - 1];
                        (
                            ps + if is_match { MATCH } else { MISMATCH },
                            pm + u32::from(is_match),
                        )
                    }
                };
                let up = {
                    let (ps, pm) =
                        band_cell(prev_lo, prev_hi, &self.prev_score, &self.prev_matches, j);
                    if ps <= NEG {
                        (NEG, 0)
                    } else {
                        (ps + GAP, pm)
                    }
                };
                let left = if j == lo {
                    (NEG, 0)
                } else {
                    let (cs, cm) = (self.curr_score[j - 1 - lo], self.curr_matches[j - 1 - lo]);
                    if cs <= NEG {
                        (NEG, 0)
                    } else {
                        (cs + GAP, cm)
                    }
                };
                let (score, matches) = if diag.0 >= up.0 && diag.0 >= left.0 {
                    diag
                } else if up.0 >= left.0 {
                    up
                } else {
                    left
                };
                self.curr_score[j - lo] = score;
                self.curr_matches[j - lo] = matches;
            }
            std::mem::swap(&mut self.prev_score, &mut self.curr_score);
            std::mem::swap(&mut self.prev_matches, &mut self.curr_matches);
            prev_lo = lo;
            prev_hi = hi;
        }
        let (final_score, final_matches) =
            band_cell(prev_lo, prev_hi, &self.prev_score, &self.prev_matches, m);
        if final_score <= NEG {
            None
        } else {
            Some(final_matches as usize)
        }
    }

    /// Same recurrence as `matches`; returns `None` as soon as no path can
    /// still reach `needed` identical pairs.
    fn matches_bounded(&mut self, short: &[u8], long: &[u8], band: usize, needed: usize) -> Option<usize> {
        let n = short.len();
        let m = long.len();
        let width = band.saturating_mul(2) + 1;
        self.ensure(width);
        self.prev_score[..width].fill(NEG);
        self.prev_matches[..width].fill(0);
        self.curr_score[..width].fill(NEG);
        self.curr_matches[..width].fill(0);
        let mut prev_lo = 0usize;
        let mut prev_hi = band.min(m);
        for j in prev_lo..=prev_hi {
            self.prev_score[j - prev_lo] = j as i32 * GAP;
        }
        // Exact remaining-composition bound: residues still unconsumed in
        // short[i..] and long[lo..] can contribute at most sum_c min(counts).
        let mut cs = [0u32; 256];
        let mut cl = [0u32; 256];
        for &r in short { cs[r as usize] += 1; }
        for &r in long { cl[r as usize] += 1; }
        let mut comp: usize = (0..256).map(|c| cs[c].min(cl[c]) as usize).sum();
        let mut long_from = 0usize;
        for i in 1..=n {
            let center = (i * m) / n;
            let lo = center.saturating_sub(band);
            let hi = (center + band).min(m);
            {
                let c = short[i - 1] as usize;
                if cs[c] <= cl[c] { comp -= 1; }
                cs[c] -= 1;
                while long_from < lo {
                    let c = long[long_from] as usize;
                    if cl[c] <= cs[c] { comp -= 1; }
                    cl[c] -= 1;
                    long_from += 1;
                }
            }
            let mut row_best = 0u32;
            for j in lo..=hi {
                let diag = if j == 0 {
                    (NEG, 0)
                } else {
                    let (ps, pm) = band_cell(prev_lo, prev_hi, &self.prev_score, &self.prev_matches, j - 1);
                    if ps <= NEG {
                        (NEG, 0)
                    } else {
                        let is_match = short[i - 1] == long[j - 1];
                        (ps + if is_match { MATCH } else { MISMATCH }, pm + u32::from(is_match))
                    }
                };
                let up = {
                    let (ps, pm) = band_cell(prev_lo, prev_hi, &self.prev_score, &self.prev_matches, j);
                    if ps <= NEG { (NEG, 0) } else { (ps + GAP, pm) }
                };
                let left = if j == lo {
                    (NEG, 0)
                } else {
                    let (cs, cm) = (self.curr_score[j - 1 - lo], self.curr_matches[j - 1 - lo]);
                    if cs <= NEG { (NEG, 0) } else { (cs + GAP, cm) }
                };
                let (score, matches) = if diag.0 >= up.0 && diag.0 >= left.0 {
                    diag
                } else if up.0 >= left.0 {
                    up
                } else {
                    left
                };
                self.curr_score[j - lo] = score;
                self.curr_matches[j - lo] = matches;
                if score > NEG && matches > row_best {
                    row_best = matches;
                }
            }
            std::mem::swap(&mut self.prev_score, &mut self.curr_score);
            std::mem::swap(&mut self.prev_matches, &mut self.curr_matches);
            prev_lo = lo;
            prev_hi = hi;
            if (row_best as usize) + (n - i).min(comp) < needed {
                #[cfg(feature = "counters")]
                DP_ABORT.fetch_add(1, AtomicOrdering::Relaxed);
                #[cfg(feature = "counters")]
                CELLS_DONE.fetch_add((i * width) as u64, AtomicOrdering::Relaxed);
                return None;
            }
        }
        #[cfg(feature = "counters")]
        CELLS_DONE.fetch_add((n * width) as u64, AtomicOrdering::Relaxed);
        let (final_score, final_matches) = band_cell(prev_lo, prev_hi, &self.prev_score, &self.prev_matches, m);
        if final_score <= NEG { None } else { Some(final_matches as usize) }
    }

    fn ensure(&mut self, width: usize) {
        if self.prev_score.len() >= width {
            return;
        }
        self.prev_score.resize(width, NEG);
        self.prev_matches.resize(width, 0);
        self.curr_score.resize(width, NEG);
        self.curr_matches.resize(width, 0);
    }
}

fn band_cell(lo: usize, hi: usize, score: &[i32], matches: &[u32], j: usize) -> (i32, u32) {
    if j < lo || j > hi {
        (NEG, 0)
    } else {
        (score[j - lo], matches[j - lo])
    }
}

thread_local! {
    static BAND: RefCell<BandScratch> = RefCell::new(BandScratch {
        prev_score: Vec::new(),
        prev_matches: Vec::new(),
        curr_score: Vec::new(),
        curr_matches: Vec::new(),
    });
}

/// Place `short` inside `long`. The diagonal that actually matches the most
/// residues wins, not the one with the most raw k-mer hits: a tandem repeat
/// votes for several offsets, and the true offset is the contiguous one.
/// `internal_gap` is set only when that best diagonal still fails to explain
/// `short` and a second diagonal, more than five residues away, is needed.
/// A clean truncation of a repetitive gene stays one diagonal.
#[allow(dead_code)]
pub fn fit_short_in_long_reference(short: &[u8], long: &[u8]) -> Option<Fit> {
    if short.is_empty() || long.len() < short.len() || short.len() < KMER {
        return None;
    }
    let mut votes: Vec<(isize, u32)> = Vec::new();
    let mut positions: std::collections::HashMap<u32, Vec<usize>> =
        std::collections::HashMap::new();
    for j in 0..=long.len() - KMER {
        let key = pack(&long[j..j + KMER])?;
        positions.entry(key).or_default().push(j);
    }
    let mut kmers = 0u32;
    for i in 0..=short.len() - KMER {
        let Some(key) = pack(&short[i..i + KMER]) else {
            continue;
        };
        let Some(hits) = positions.get(&key) else {
            continue;
        };
        kmers += 1;
        for &j in hits {
            let diagonal = j as isize - i as isize;
            if let Some(slot) = votes.iter_mut().find(|(d, _)| *d == diagonal) {
                slot.1 += 1;
            } else {
                votes.push((diagonal, 1));
            }
        }
    }
    if kmers == 0 || votes.is_empty() {
        return None;
    }
    votes.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    let best_votes = votes[0].1;
    // A handful of competitive diagonals is enough. Scoring every offset of a
    // low-complexity protein would redo the alignment the k-mers already ranked.
    let contenders: Vec<isize> = votes
        .iter()
        .filter(|(_, votes)| votes.saturating_mul(5) >= best_votes)
        .map(|(diagonal, _)| *diagonal)
        .take(8)
        .collect();
    let mut chosen: Option<(isize, Fit)> = None;
    for diagonal in contenders {
        let placed = place(short, long, diagonal);
        let Some(placed) = placed else {
            continue;
        };
        let replace = match chosen {
            None => true,
            Some((_, previous)) => {
                placed.identity > previous.identity + 1e-12
                    || ((placed.identity - previous.identity).abs() <= 1e-12
                        && placed.coverage > previous.coverage + 1e-12)
            }
        };
        if replace {
            chosen = Some((diagonal, placed));
        }
    }
    let (best_diag, mut fit) = chosen?;
    // Two strong diagonals only mean an internal insertion when one diagonal
    // cannot already hold `short`. A second offset that also matches is a
    // repeat of the same peptide, and the longer allele is still the gene.
    if fit.identity < 0.95 || fit.coverage < 0.90 {
        if let Some((second, _)) = votes.iter().find(|(diagonal, votes)| {
            *diagonal != best_diag
                && (*diagonal - best_diag).unsigned_abs() > 5
                && votes.saturating_mul(5) >= best_votes
        }) {
            fit.internal_gap = (*second - best_diag).unsigned_abs();
        }
    }
    Some(fit)
}

struct FitScratch {
    keyed: Vec<(u32, u32)>,
    votes: Vec<u32>,
    touched: Vec<u32>,
    ranked: Vec<(isize, u32)>,
}

thread_local! {
    static FIT: RefCell<FitScratch> = RefCell::new(FitScratch { keyed: Vec::new(), votes: Vec::new(), touched: Vec::new(), ranked: Vec::new() });
}

/// The sorted `(4-mer, position)` table of one `long` sequence, built once
/// and reused for every `short` fitted into it. `None` when `long` has a
/// window that does not pack, which is exactly when `fit_short_in_long`
/// returns `None` for every `short` of a valid length.
pub struct LongIndex {
    keyed: Option<Vec<(u32, u32)>>,
}

impl LongIndex {
    pub fn new(long: &[u8]) -> Self {
        if long.len() < KMER {
            return Self { keyed: Some(Vec::new()) };
        }
        let mut keyed = Vec::with_capacity(long.len() - KMER + 1);
        for j in 0..=long.len() - KMER {
            let Some(key) = pack(&long[j..j + KMER]) else {
                return Self { keyed: None };
            };
            keyed.push((key, j as u32));
        }
        keyed.sort_unstable();
        Self { keyed: Some(keyed) }
    }
}

/// PERF PROTOTYPE: identical result to `fit_short_in_long_reference`.
/// A sorted (k-mer, position) array replaces the per-call `HashMap<u32, Vec>`
/// and a dense per-diagonal counter replaces the linear `votes.find`. The
/// final vote order is the same total order (votes desc, diagonal asc).
pub fn fit_short_in_long(short: &[u8], long: &[u8]) -> Option<Fit> {
    if short.is_empty() || long.len() < short.len() || short.len() < KMER {
        return None;
    }
    FIT.with(|cell| {
        let scratch = &mut *cell.borrow_mut();
        let mut keyed = std::mem::take(&mut scratch.keyed);
        keyed.clear();
        let mut packable = true;
        for j in 0..=long.len() - KMER {
            let Some(key) = pack(&long[j..j + KMER]) else {
                packable = false;
                break;
            };
            keyed.push((key, j as u32));
        }
        let result = if packable {
            keyed.sort_unstable();
            fit_with_keyed(scratch, short, long, &keyed)
        } else {
            None
        };
        scratch.keyed = keyed;
        result
    })
}

/// `fit_short_in_long(short, long)` with `long`'s table prebuilt; identical
/// result (the same early checks, then the same vote and placement logic on
/// the same sorted table).
pub fn fit_short_in_long_indexed(short: &[u8], long: &[u8], index: &LongIndex) -> Option<Fit> {
    if short.is_empty() || long.len() < short.len() || short.len() < KMER {
        return None;
    }
    let keyed = index.keyed.as_deref()?;
    FIT.with(|cell| fit_with_keyed(&mut cell.borrow_mut(), short, long, keyed))
}

fn fit_with_keyed(
    scratch: &mut FitScratch,
    short: &[u8],
    long: &[u8],
    keyed: &[(u32, u32)],
) -> Option<Fit> {
    {
        // diagonal = j - i, in [-(short.len()-KMER), long.len()-KMER]
        let offset = short.len() - KMER;
        let span = offset + long.len() - KMER + 1;
        if scratch.votes.len() < span {
            scratch.votes.resize(span, 0);
        }
        scratch.touched.clear();
        let mut kmers = 0u32;
        for i in 0..=short.len() - KMER {
            let Some(key) = pack(&short[i..i + KMER]) else { continue };
            let start = keyed.partition_point(|&(k, _)| k < key);
            if start >= keyed.len() || keyed[start].0 != key {
                continue;
            }
            kmers += 1;
            let mut p = start;
            while p < keyed.len() && keyed[p].0 == key {
                let j = keyed[p].1 as usize;
                let slot = j + offset - i;
                if scratch.votes[slot] == 0 {
                    scratch.touched.push(slot as u32);
                }
                scratch.votes[slot] += 1;
                p += 1;
            }
        }
        scratch.ranked.clear();
        for &slot in &scratch.touched {
            let slot = slot as usize;
            scratch.ranked.push((slot as isize - offset as isize, scratch.votes[slot]));
            scratch.votes[slot] = 0;
        }
        if kmers == 0 || scratch.ranked.is_empty() {
            return None;
        }
        let votes = &mut scratch.ranked;
        votes.sort_unstable_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        let best_votes = votes[0].1;
        let mut chosen: Option<(isize, Fit)> = None;
        // Same contender selection as the reference: filter then take(8).
        let contenders: Vec<isize> = votes
            .iter()
            .filter(|(_, v)| v.saturating_mul(5) >= best_votes)
            .map(|(d, _)| *d)
            .take(8)
            .collect();
        for diagonal in contenders {
            let Some(placed) = place(short, long, diagonal) else { continue };
            let replace = match chosen {
                None => true,
                Some((_, previous)) => {
                    placed.identity > previous.identity + 1e-12
                        || ((placed.identity - previous.identity).abs() <= 1e-12
                            && placed.coverage > previous.coverage + 1e-12)
                }
            };
            if replace {
                chosen = Some((diagonal, placed));
            }
        }
        let (best_diag, mut fit) = chosen?;
        if fit.identity < 0.95 || fit.coverage < 0.90 {
            if let Some((second, _)) = votes.iter().find(|(diagonal, v)| {
                *diagonal != best_diag
                    && (*diagonal - best_diag).unsigned_abs() > 5
                    && v.saturating_mul(5) >= best_votes
            }) {
                fit.internal_gap = (*second - best_diag).unsigned_abs();
            }
        }
        Some(fit)
    }
}

fn place(short: &[u8], long: &[u8], diagonal: isize) -> Option<Fit> {
    let mut matches = 0usize;
    let mut compared = 0usize;
    for i in 0..short.len() {
        let j = i as isize + diagonal;
        if j < 0 || j >= long.len() as isize {
            continue;
        }
        compared += 1;
        if short[i] == long[j as usize] {
            matches += 1;
        }
    }
    if compared == 0 {
        return None;
    }
    let long_start = usize::try_from(diagonal).unwrap_or(0);
    let long_end = long_start.saturating_add(compared).min(long.len());
    Some(Fit {
        identity: matches as f64 / short.len() as f64,
        coverage: compared as f64 / short.len() as f64,
        long_start,
        long_end,
        internal_gap: 0,
    })
}

fn pack(kmer: &[u8]) -> Option<u32> {
    let mut value = 0u32;
    for &residue in kmer {
        let code = match residue {
            b'A' => 0,
            b'C' => 1,
            b'D' => 2,
            b'E' => 3,
            b'F' => 4,
            b'G' => 5,
            b'H' => 6,
            b'I' => 7,
            b'K' => 8,
            b'L' => 9,
            b'M' => 10,
            b'N' => 11,
            b'P' => 12,
            b'Q' => 13,
            b'R' => 14,
            b'S' => 15,
            b'T' => 16,
            b'V' => 17,
            b'W' => 18,
            b'Y' => 19,
            _ => return None,
        };
        value = (value << 5) | code;
    }
    Some(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rng(seed: u64) -> impl FnMut() -> u64 {
        let mut state = seed | 1;
        move || { state ^= state << 13; state ^= state >> 7; state ^= state << 17; state }
    }

    fn random_pair(next: &mut impl FnMut() -> u64) -> (Vec<u8>, Vec<u8>) {
        let alphabet: &[u8] = if next() % 4 == 0 { b"ACDEK" } else { b"ACDEFGHIKLMNPQRSTVWYX*" };
        let len = 4 + (next() % 200) as usize;
        let a: Vec<u8> = (0..len).map(|_| alphabet[(next() % alphabet.len() as u64) as usize]).collect();
        let mut b = a.clone();
        match next() % 4 {
            0 => { let cut = (next() as usize) % b.len(); b.drain(..cut.min(b.len().saturating_sub(4))); }
            1 => { for _ in 0..(next() % 40) { let i = (next() as usize) % b.len(); b[i] = alphabet[(next() % alphabet.len() as u64) as usize]; } }
            2 => { let ins: Vec<u8> = (0..(next() % 30)).map(|_| alphabet[(next() % alphabet.len() as u64) as usize]).collect(); let at = (next() as usize) % b.len(); b.splice(at..at, ins); }
            _ => { b = (0..(4 + next() % 200) as usize).map(|_| alphabet[(next() % alphabet.len() as u64) as usize]).collect(); }
        }
        (a, b)
    }

    #[test]
    fn perf_indexed_fit_matches_reference_fit() {
        let mut next = rng(0x1234_5678);
        for _ in 0..2_000 {
            let (a, b) = random_pair(&mut next);
            let ia = LongIndex::new(&a);
            let ib = LongIndex::new(&b);
            assert_eq!(fit_short_in_long_indexed(&a, &b, &ib), fit_short_in_long_reference(&a, &b));
            assert_eq!(fit_short_in_long_indexed(&b, &a, &ia), fit_short_in_long_reference(&b, &a));
        }
    }

    #[test]
    fn perf_fast_fit_matches_reference_fit() {
        let mut next = rng(0x9e3779b97f4a7c15);
        for _ in 0..2_000 {
            let (a, b) = random_pair(&mut next);
            assert_eq!(fit_short_in_long(&a, &b), fit_short_in_long_reference(&a, &b), "{a:?} {b:?}");
            assert_eq!(fit_short_in_long(&b, &a), fit_short_in_long_reference(&b, &a), "{b:?} {a:?}");
        }
    }

    #[test]
    fn perf_bounded_identity_matches_whenever_threshold_reachable() {
        let mut next = rng(0x1234_5678_9abc_def1);
        for _ in 0..2_000 {
            let (a, b) = random_pair(&mut next);
            let exact = global_identity(&a, &b);
            for threshold in [0.98, 0.70, 0.40, 0.0] {
                let bounded = global_identity_at_least(&a, &b, threshold);
                match exact {
                    Some(id) if id + f64::EPSILON >= threshold => assert_eq!(bounded, Some(id)),
                    _ => assert!(bounded.is_none() || bounded == exact, "{a:?} {b:?} {threshold}"),
                }
            }
        }
    }

    #[test]
    fn identical_and_one_mismatch() {
        assert_eq!(global_identity(b"ACDEFG", b"ACDEFG"), Some(1.0));
        let id = global_identity(b"ACDEFGHIKL", b"ACDEAGHIKL").unwrap();
        assert!((id - 0.9).abs() < 1e-9);
    }

    #[test]
    fn insertion_keeps_identity_over_the_shorter_sequence() {
        let id = global_identity(b"ACDEFGHIKL", b"ACDEFFGHIKL").unwrap();
        assert!(id > 0.9, "{id}");
    }

    #[test]
    fn prefix_is_a_clean_extension() {
        let long = b"MKTLLYVWAACDEFGHIKLMNPQR";
        let short = &long[..12];
        let fit = fit_short_in_long(short, long).unwrap();
        assert!(fit.identity > 0.95);
        assert_eq!(fit.internal_gap, 0);
        assert_eq!(fit.extra_prefix(long.len()), 0);
        assert!(fit.extra_suffix(long.len()) > 0);
    }

    #[test]
    fn internal_insertion_splits_the_diagonal() {
        let head = b"MKTLLYVWAACDEFGHIKLM";
        let tail = b"NPQRSTVWYACDEFGHIKLM";
        let junk = b"WWWWWWWWWWWWWWWWWWWW";
        let mut long = head.to_vec();
        long.extend_from_slice(junk);
        long.extend_from_slice(tail);
        let mut short = head.to_vec();
        short.extend_from_slice(tail);
        let fit = fit_short_in_long(&short, &long).unwrap();
        assert!(
            fit.internal_gap > 5,
            "gap {} identity {}",
            fit.internal_gap,
            fit.identity
        );
    }

    /// The banded forward pass must agree exactly with a brute-force full
    /// matrix for pairs within the band, since it is the same recurrence and
    /// tie-break rule, just computed without a traceback. This guards the
    /// `global_matches` -> `banded_matches` rewrite against any drift.
    #[test]
    fn banded_identity_matches_brute_force_full_alignment() {
        fn brute_force_identity(a: &[u8], b: &[u8]) -> f64 {
            let (short, long) = if a.len() <= b.len() { (a, b) } else { (b, a) };
            let n = short.len();
            let m = long.len();
            let mut score = vec![vec![0i32; m + 1]; n + 1];
            let mut matches = vec![vec![0u32; m + 1]; n + 1];
            for i in 0..=n {
                score[i][0] = i as i32 * GAP;
            }
            for j in 0..=m {
                score[0][j] = j as i32 * GAP;
            }
            for i in 1..=n {
                for j in 1..=m {
                    let is_match = short[i - 1] == long[j - 1];
                    let diag = (
                        score[i - 1][j - 1] + if is_match { MATCH } else { MISMATCH },
                        matches[i - 1][j - 1] + is_match as u32,
                    );
                    let up = (score[i - 1][j] + GAP, matches[i - 1][j]);
                    let left = (score[i][j - 1] + GAP, matches[i][j - 1]);
                    let (best_score, best_matches) = if diag.0 >= up.0 && diag.0 >= left.0 {
                        diag
                    } else if up.0 >= left.0 {
                        up
                    } else {
                        left
                    };
                    score[i][j] = best_score;
                    matches[i][j] = best_matches;
                }
            }
            matches[n][m] as f64 / n as f64
        }

        let mut state = 0x1234_5678_9abc_def0u64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let alphabet = b"ACDEFGHIKLMNPQRSTVWY";
        for _ in 0..80 {
            let base_len = 20 + (next() % 60) as usize;
            let drift = (next() % 5) as usize;
            let a: Vec<u8> = (0..base_len)
                .map(|_| alphabet[(next() % 20) as usize])
                .collect();
            let mut b = a.clone();
            for _ in 0..(next() % 5) {
                let idx = (next() as usize) % b.len();
                b[idx] = alphabet[(next() % 20) as usize];
            }
            if drift > 0 && b.len() > drift {
                if next() % 2 == 0 {
                    b.truncate(b.len() - drift);
                } else {
                    for _ in 0..drift {
                        b.push(alphabet[(next() % 20) as usize]);
                    }
                }
            }
            let expected = brute_force_identity(&a, &b);
            let actual = global_identity(&a, &b).unwrap();
            assert!(
                (expected - actual).abs() < 1e-9,
                "a={a:?} b={b:?} expected={expected} actual={actual}"
            );
        }
    }
}
