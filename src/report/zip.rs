use anyhow::{Context, Result};
use std::fs::{self, File};
use std::io::{Read, Write};
use std::path::Path;
use zip::write::SimpleFileOptions;
use zip::{CompressionMethod, ZipWriter};

pub fn write_zip(out_dir: &Path, sample_name: &str) -> Result<()> {
    let root = format!("{}_fastqc", sample_name);
    let zip_name = format!("{}_fastqc.zip", sample_name);
    let zip_path = out_dir.join(&zip_name);
    let tmp_path = out_dir.join(format!("{}.tmp", zip_name));

    let file = File::create(&tmp_path)
        .with_context(|| format!("failed to create {}", tmp_path.display()))?;
    let mut zip = ZipWriter::new(file);
    let result = write_zip_entries(&mut zip, out_dir, &root);

    match result.and_then(|_| zip.finish().with_context(|| "failed to finalize zip")) {
        Ok(_) => {
            fs::rename(&tmp_path, &zip_path)
                .with_context(|| format!("failed to move zip to {}", zip_path.display()))?;
            Ok(())
        }
        Err(e) => {
            let _ = fs::remove_file(&tmp_path);
            Err(e)
        }
    }
}

fn write_zip_entries(zip: &mut ZipWriter<File>, out_dir: &Path, root: &str) -> Result<()> {
    let options = SimpleFileOptions::default()
        .compression_method(CompressionMethod::Deflated)
        .last_modified_time(zip::DateTime::from_date_and_time(1980, 1, 1, 0, 0, 0).unwrap());

    zip.add_directory(format!("{}/", root), options)
        .with_context(|| "failed to add directory entry to zip")?;

    let files = ["fastqc_data.txt", "summary.txt", "fastqc_report.html"];

    for name in files {
        let src_path = out_dir.join(root).join(name);
        let zip_path = format!("{}/{}", root, name);
        add_file(zip, &src_path, &zip_path, options)
            .with_context(|| format!("failed to add {} to zip", name))?;
    }
    Ok(())
}

fn add_file(
    zip: &mut ZipWriter<File>,
    src_path: &Path,
    zip_path: &str,
    options: SimpleFileOptions,
) -> Result<()> {
    let mut file =
        File::open(src_path).with_context(|| format!("failed to open {}", src_path.display()))?;
    zip.start_file(zip_path, options)?;
    let mut buf = [0u8; 8192];
    loop {
        let n = file.read(&mut buf)?;
        if n == 0 {
            break;
        }
        zip.write_all(&buf[..n])?;
    }
    Ok(())
}
