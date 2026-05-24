use crate::core::metrics::FinalMetrics;
use crate::core::model::{FinalizeContext, Mode, Status};
use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

pub fn write(path: &Path, metrics: &FinalMetrics, ctx: &FinalizeContext) -> Result<()> {
    let mut w = BufWriter::new(File::create(path).with_context(|| "create summary.txt failed")?);

    let file = &ctx.file_name;

    let line = |w: &mut BufWriter<File>, status: Status, name: &str| -> Result<()> {
        writeln!(w, "{}\t{}\t{}", status.as_str_upper(), name, file)?;
        Ok(())
    };

    line(&mut w, metrics.statuses.basic, "Basic Statistics")?;

    match ctx.mode {
        Mode::Short => {
            line(
                &mut w,
                metrics.statuses.per_base_qual,
                "Per base sequence quality",
            )?;
            line(
                &mut w,
                metrics.statuses.per_seq_qual,
                "Per sequence quality scores",
            )?;
            line(
                &mut w,
                metrics.statuses.per_base_content,
                "Per base sequence content",
            )?;
            line(
                &mut w,
                metrics.statuses.per_seq_gc,
                "Per sequence GC content",
            )?;
            line(&mut w, metrics.statuses.per_base_n, "Per base N content")?;
            line(
                &mut w,
                metrics.statuses.length_dist,
                "Sequence Length Distribution",
            )?;
            line(
                &mut w,
                metrics.statuses.duplication,
                "Sequence Duplication Levels",
            )?;
            line(
                &mut w,
                metrics.statuses.overrepresented,
                "Overrepresented sequences",
            )?;
            line(&mut w, metrics.statuses.adapter_content, "Adapter Content")?;
            #[cfg(not(feature = "no-kmer"))]
            line(&mut w, metrics.statuses.kmer_content, "Kmer Content")?;
        }
        Mode::Long => {
            line(
                &mut w,
                metrics.statuses.length_dist,
                "Sequence Length Distribution",
            )?;
            line(
                &mut w,
                metrics.statuses.per_seq_qual,
                "Per sequence quality scores",
            )?;
            line(
                &mut w,
                metrics.statuses.per_seq_gc,
                "Per sequence GC content",
            )?;
            line(&mut w, metrics.statuses.per_seq_n, "Per sequence N content")?;
            line(&mut w, metrics.statuses.adapter_content, "Adapter Content")?;
        }
    }

    Ok(())
}
