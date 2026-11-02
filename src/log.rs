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
    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
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

/// PERF-INSTRUMENTATION (rustaroo_perf copy only): stage wall time + RSS.
pub fn stages_enabled() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("RUSTROO_STAGES").is_some())
}

pub fn stage(label: &str) {
    if !stages_enabled() {
        return;
    }
    use std::sync::OnceLock;
    use std::time::Instant;
    static START: OnceLock<std::sync::Mutex<(Instant, Instant)>> = OnceLock::new();
    let cell = START.get_or_init(|| std::sync::Mutex::new((Instant::now(), Instant::now())));
    let mut guard = cell.lock().unwrap();
    let now = Instant::now();
    let since_last = now.duration_since(guard.1).as_secs_f64();
    let total = now.duration_since(guard.0).as_secs_f64();
    guard.1 = now;
    let status = std::fs::read_to_string("/proc/self/status").unwrap_or_default();
    let field = |name: &str| status.lines().find(|l| l.starts_with(name)).map(|l| l.split_whitespace().nth(1).unwrap_or("0").parse::<u64>().unwrap_or(0)).unwrap_or(0);
    let cpu = std::fs::read_to_string("/proc/self/stat").ok().and_then(|t| { let r = t.rsplit_once(')')?.1.to_string(); let f: Vec<&str> = r.split_whitespace().collect(); Some((f[11].parse::<f64>().ok()? + f[12].parse::<f64>().ok()?) / 100.0) }).unwrap_or(0.0);
    eprintln!("STAGE\t{label}\t{since_last:.3}\t{total:.3}\tcpu_total={cpu:.2}\trss_mb={}\thwm_mb={}", field("VmRSS:") / 1024, field("VmHWM:") / 1024);
}
