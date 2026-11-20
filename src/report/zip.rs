use anyhow::{Context, Result};
use std::fs::{self, File};
use std::io::{Read, Write};
use std::path::Path;
use zip::write::SimpleFileOptions;
use zip::{CompressionMethod, ZipWriter};

/// Fixed archive timestamp so repeated runs produce byte-identical bundles.
pub fn deterministic_time() -> zip::DateTime {
    zip::DateTime::from_date_and_time(1980, 1, 1, 0, 0, 0).expect("valid fixed timestamp")
}

pub fn options() -> SimpleFileOptions {
    SimpleFileOptions::default()
        .compression_method(CompressionMethod::Deflated)
        .last_modified_time(deterministic_time())
}

/// Bundle the report directory, mirroring FastQC's `<sample>_fastqc.zip` layout.
pub fn write_zip(out_dir: &Path, sample_name: &str, extra_files: &[String]) -> Result<()> {
    let source = out_dir.join(format!("{sample_name}_fastqc"));
    write_zip_from(out_dir, &source, sample_name, extra_files)
}

pub(crate) fn write_zip_from(
    out_dir: &Path,
    source: &Path,
    sample_name: &str,
    extra_files: &[String],
) -> Result<()> {
    let root = format!("{}_fastqc", sample_name);
    let zip_name = format!("{}_fastqc.zip", sample_name);
    let zip_path = out_dir.join(&zip_name);
    let mut temporary = tempfile::NamedTempFile::new_in(out_dir).with_context(|| {
        format!(
            "failed to create a temporary archive in {}",
            out_dir.display()
        )
    })?;
    let mut zip = ZipWriter::new(temporary.as_file_mut());
    write_zip_entries(&mut zip, source, &root, extra_files)?;
    zip.finish().context("failed to finalize zip")?;
    temporary
        .persist(&zip_path)
        .with_context(|| format!("failed to publish {}", zip_path.display()))?;
    Ok(())
}

fn write_zip_entries(
    zip: &mut ZipWriter<&mut File>,
    source: &Path,
    root: &str,
    extra_files: &[String],
) -> Result<()> {
    let options = options();

    zip.add_directory(format!("{}/", root), options)
        .with_context(|| "failed to add directory entry to zip")?;

    let base = ["fastqc_data.txt", "summary.txt", "fastqc_report.html"];
    let mut names: Vec<String> = base.iter().map(|s| (*s).to_string()).collect();
    names.extend(extra_files.iter().cloned());
    // Fixed order keeps the archive reproducible.
    names.sort();
    names.dedup();

    let mut added_dirs: Vec<String> = Vec::new();
    for name in names {
        let src_path = source.join(&name);
        if !src_path.is_file() {
            continue;
        }
        // Create parent directory entries for nested files (latex/…).
        let normalized = name.replace('\\', "/");
        if let Some(parent) = normalized.rsplit_once('/').map(|(p, _)| p) {
            let mut acc = String::new();
            for part in parent.split('/') {
                acc.push_str(part);
                acc.push('/');
                if !added_dirs.contains(&acc) {
                    zip.add_directory(format!("{}/{}", root, acc), options)?;
                    added_dirs.push(acc.clone());
                }
            }
        }
        let zip_path = format!("{}/{}", root, normalized);
        add_file(zip, &src_path, &zip_path, options)
            .with_context(|| format!("failed to add {} to zip", normalized))?;
    }
    Ok(())
}

fn add_file(
    zip: &mut ZipWriter<&mut File>,
    src_path: &Path,
    zip_path: &str,
    options: SimpleFileOptions,
) -> Result<()> {
    let mut file =
        File::open(src_path).with_context(|| format!("failed to open {}", src_path.display()))?;
    zip.start_file(zip_path, options)?;
    let mut buf = [0u8; 64 * 1024];
    loop {
        let n = file.read(&mut buf)?;
        if n == 0 {
            break;
        }
        zip.write_all(&buf[..n])?;
    }
    Ok(())
}

/// Relative paths of every file under `dir`, using `/` separators.
pub fn list_files_recursive(dir: &Path, prefix: &str) -> Result<Vec<String>> {
    let mut out = Vec::new();
    if !dir.is_dir() {
        return Ok(out);
    }
    let mut entries: Vec<_> = fs::read_dir(dir)
        .with_context(|| format!("failed to list {}", dir.display()))?
        .collect::<std::result::Result<Vec<_>, _>>()?;
    entries.sort_by_key(|e| e.file_name());
    for entry in entries {
        let name = entry.file_name().to_string_lossy().into_owned();
        let rel = if prefix.is_empty() {
            name.clone()
        } else {
            format!("{}/{}", prefix, name)
        };
        let path = entry.path();
        if path.is_dir() {
            out.extend(list_files_recursive(&path, &rel)?);
        } else {
            out.push(rel);
        }
    }
    Ok(out)
}
