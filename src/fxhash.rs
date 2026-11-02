//! A fast, fixed (unseeded) hasher for the pipeline's internal maps.
//!
//! The standard `RandomState` is seeded randomly per process, so no output
//! of the original program can depend on map iteration order (it would differ
//! between runs). Replacing it with a fixed hash therefore cannot change any
//! output; it only removes SipHash's cost. Every map iteration that feeds an
//! output is still followed by a sort or a total-order selection.
//!
//! Same multiply-add scheme as `rustc-hash` 2.x (`FxHasher`).

use std::hash::{BuildHasherDefault, Hasher};

pub type FxBuildHasher = BuildHasherDefault<FxHasher>;
pub type HashMap<K, V> = std::collections::HashMap<K, V, FxBuildHasher>;
pub type HashSet<T> = std::collections::HashSet<T, FxBuildHasher>;

const K: u64 = 0xf1357aea2e62a9c5;

#[derive(Clone, Copy)]
pub struct FxHasher {
    hash: u64,
}

impl Default for FxHasher {
    #[inline]
    fn default() -> Self {
        Self { hash: seed() }
    }
}

/// Always 0 in a normal build. The `hashseed` feature (verification builds
/// only) reads `RUSTROO_HASH_SEED`, so runs with different seeds visit every
/// map in a different order; identical outputs across seeds show that no
/// output depends on map iteration order.
#[cfg(not(feature = "hashseed"))]
#[inline]
fn seed() -> u64 {
    0
}

#[cfg(feature = "hashseed")]
fn seed() -> u64 {
    static SEED: std::sync::OnceLock<u64> = std::sync::OnceLock::new();
    *SEED.get_or_init(|| {
        std::env::var("RUSTROO_HASH_SEED")
            .ok()
            .and_then(|value| value.parse().ok())
            .unwrap_or(0)
    })
}

impl FxHasher {
    #[inline]
    fn add(&mut self, value: u64) {
        self.hash = self.hash.wrapping_add(value).wrapping_mul(K);
    }
}

impl Hasher for FxHasher {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        let mut chunks = bytes.chunks_exact(8);
        for chunk in &mut chunks {
            self.add(u64::from_le_bytes(chunk.try_into().unwrap()));
        }
        let rest = chunks.remainder();
        if !rest.is_empty() {
            let mut tail = [0u8; 8];
            tail[..rest.len()].copy_from_slice(rest);
            self.add(u64::from_le_bytes(tail) ^ ((rest.len() as u64) << 59));
        }
    }

    #[inline]
    fn write_u8(&mut self, value: u8) {
        self.add(u64::from(value));
    }

    #[inline]
    fn write_u16(&mut self, value: u16) {
        self.add(u64::from(value));
    }

    #[inline]
    fn write_u32(&mut self, value: u32) {
        self.add(u64::from(value));
    }

    #[inline]
    fn write_u64(&mut self, value: u64) {
        self.add(value);
    }

    #[inline]
    fn write_usize(&mut self, value: usize) {
        self.add(value as u64);
    }

    #[inline]
    fn finish(&self) -> u64 {
        self.hash.rotate_left(26)
    }
}
