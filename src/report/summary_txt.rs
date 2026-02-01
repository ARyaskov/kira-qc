use crate::core::engine::RunOutput;
use crate::core::model::Mode;
use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

pub fn write(path: &Path, output: &RunOutput) -> Result<()> {
    let metrics = output.agg.finalize(&output.ctx);
    let mut w = BufWriter::new(File::create(path).with_context(|| "create summary.txt failed")?);

    let file = &output.ctx.file_name;

    writeln!(
        w,
        "{}\t{}\t{}",
        metrics.statuses.basic.as_str_upper(),
        "Basic Statistics",
        file
    )?;

    match output.ctx.mode {
        Mode::Short => {
            writeln!(
                w,
                "{}\t{}\t{}",
                metrics.statuses.per_base_qual.as_str_upper(),
                "Per base sequence quality",
                file
            )?;
            writeln!(
                w,
                "{}\t{}\t{}",
                metrics.statuses.per_seq_qual.as_str_upper(),
                "Per sequence quality scores",
                file
            )?;
            writeln!(
                w,
                "{}\t{}\t{}",
                metrics.statuses.per_base_content.as_str_upper(),
                "Per base sequence content",
                file
            )?;
            writeln!(
                w,
                "{}\t{}\t{}",
                metrics.statuses.per_seq_gc.as_str_upper(),
                "Per sequence GC content",
                file
            )?;
            writeln!(
                w,
                "{}\t{}\t{}",
                metrics.statuses.per_base_n.as_str_upper(),
                "Per base N content",
                file
            )?;
            writeln!(
                w,
                "{}\t{}\t{}",
                metrics.statuses.length_dist.as_str_upper(),
                "Sequence Length Distribution",
                file
            )?;
            writeln!(
                w,
                "{}\t{}\t{}",
                metrics.statuses.duplication.as_str_upper(),
                "Sequence Duplication Levels",
                file
            )?;
            writeln!(
                w,
                "{}\t{}\t{}",
                metrics.statuses.overrepresented.as_str_upper(),
                "Overrepresented sequences",
                file
            )?;
            writeln!(
                w,
                "{}\t{}\t{}",
                metrics.statuses.adapter_content.as_str_upper(),
                "Adapter Content",
                file
            )?;
            #[cfg(not(feature = "no-kmer"))]
            writeln!(
                w,
                "{}\t{}\t{}",
                metrics.statuses.kmer_content.as_str_upper(),
                "Kmer Content",
                file
            )?;
        }
        Mode::Long => {
            writeln!(
                w,
                "{}\t{}\t{}",
                metrics.statuses.length_dist.as_str_upper(),
                "Sequence Length Distribution",
                file
            )?;
            writeln!(
                w,
                "{}\t{}\t{}",
                metrics.statuses.per_seq_qual.as_str_upper(),
                "Per sequence quality scores",
                file
            )?;
            writeln!(
                w,
                "{}\t{}\t{}",
                metrics.statuses.per_seq_gc.as_str_upper(),
                "Per sequence GC content",
                file
            )?;
            writeln!(
                w,
                "{}\t{}\t{}",
                metrics.statuses.per_seq_n.as_str_upper(),
                "Per sequence N content",
                file
            )?;
            writeln!(
                w,
                "{}\t{}\t{}",
                metrics.statuses.adapter_content.as_str_upper(),
                "Adapter Content",
                file
            )?;
        }
    }

    Ok(())
}
