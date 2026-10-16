//! Run record written beside the pangenome, separate from the terminal.
//!
//! The terminal keeps the short stage lines. This file keeps those lines too,
//! then the decisions the stage lines only count: every input, every removed
//! family, and every reference. Writing is skipped entirely when no run has
//! opened a file, so tests pay nothing.

use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::Mutex;

use crate::error::Result;

static LOG: Mutex<Option<BufWriter<File>>> = Mutex::new(None);

pub fn open(path: &Path) -> Result<()> {
    if let Some(parent) = path.parent().filter(|parent| !parent.as_os_str().is_empty()) {
        fs::create_dir_all(parent)?;
    }
    let mut out = BufWriter::new(File::create(path)?);
    writeln!(out, "# rustroo log")?;
    *lock() = Some(out);
    Ok(())
}

pub fn close() {
    if let Some(mut out) = lock().take() {
        let _ = out.flush();
    }
}

pub fn enabled() -> bool {
    lock().is_some()
}

/// A line for the file only. The terminal is unchanged.
pub fn record(message: &str) {
    if let Some(out) = lock().as_mut() {
        let _ = writeln!(out, "{message}");
    }
}

/// The same short line the terminal already shows, copied into the file.
pub fn write_progress(message: &str) {
    eprintln!("{message}");
    if let Some(out) = lock().as_mut() {
        let _ = writeln!(out, "{message}");
        let _ = out.flush();
    }
}

pub fn flush() {
    if let Some(out) = lock().as_mut() {
        let _ = out.flush();
    }
}

fn lock() -> std::sync::MutexGuard<'static, Option<BufWriter<File>>> {
    LOG.lock().unwrap_or_else(|poisoned| poisoned.into_inner())
}

/// Terminal progress that is also the next line of the run log.
macro_rules! progress {
    ($($arg:tt)*) => {{
        $crate::log::write_progress(&format!($($arg)*));
    }};
}

pub(crate) use progress;
