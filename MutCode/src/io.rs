use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use sha2::{Digest, Sha256};

use crate::{fail, Result};

pub fn reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let mut file = BufReader::new(File::open(path)?);
    let gzip = file.fill_buf()?.starts_with(&[0x1f, 0x8b]);
    if gzip {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(file))
    }
}

pub fn digest(sequence: &[u8]) -> String {
    format!("{:x}", Sha256::digest(sequence))
}

pub fn file_digest(path: &Path) -> Result<String> {
    let mut file = File::open(path)?;
    let mut hasher = Sha256::new();
    let mut buffer = [0; 65536];
    loop {
        let count = file.read(&mut buffer)?;
        if count == 0 {
            break;
        }
        hasher.update(&buffer[..count]);
    }
    Ok(format!("{:x}", hasher.finalize()))
}

pub fn is_dna(sequence: &[u8]) -> bool {
    sequence.iter().all(|base| {
        matches!(
            base.to_ascii_uppercase(),
            b'A' | b'C'
                | b'G'
                | b'T'
                | b'R'
                | b'Y'
                | b'S'
                | b'W'
                | b'K'
                | b'M'
                | b'B'
                | b'D'
                | b'H'
                | b'V'
                | b'N'
        )
    })
}

pub fn identifier(value: &str) -> Result<()> {
    if value.is_empty() || value.chars().any(char::is_control) {
        return fail("empty identifier or control character in identifier");
    }
    Ok(())
}

// A same-directory staging path keeps publication on the same filesystem.
pub fn staging_path(path: &Path) -> Result<PathBuf> {
    let name = path.file_name().ok_or("output path has no filename")?;
    let stamp = SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos();
    Ok(path.with_file_name(format!(
        ".{}.{}.{}.partial",
        name.to_string_lossy(),
        std::process::id(),
        stamp
    )))
}

pub fn atomic_write(path: &Path, produce: impl FnOnce(&mut dyn Write) -> Result<()>) -> Result<()> {
    if path.exists() {
        return fail(format!("output already exists: {}", path.display()));
    }
    let temp = staging_path(path)?;
    let result = (|| {
        let mut file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temp)?;
        if path.extension().is_some_and(|ext| ext == "gz") {
            let mut writer = GzEncoder::new(&mut file, Compression::fast());
            produce(&mut writer)?;
            writer.finish()?;
        } else {
            let mut writer = BufWriter::new(&mut file);
            produce(&mut writer)?;
            writer.flush()?;
        }
        file.sync_all()?;
        // Unlike rename, hard_link refuses to overwrite an existing destination.
        fs::hard_link(&temp, path)?;
        Ok(())
    })();
    if temp.exists() {
        if let Err(error) = fs::remove_file(&temp) {
            return fail(format!(
                "temporary-file cleanup failed at {}: {error}; operation result: {result:?}",
                temp.display()
            ));
        }
    }
    result
}
