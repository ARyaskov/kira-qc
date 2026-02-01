use crate::core::engine::RunOutput;
use crate::core::model::Mode;
use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

pub fn write(path: &Path, output: &RunOutput) -> Result<()> {
    let metrics = output.agg.finalize(&output.ctx);
    let mut w =
        BufWriter::new(File::create(path).with_context(|| "create fastqc_data.txt failed")?);

    write_basic(&mut w, &metrics, &output.ctx.file_name)?;
    match output.ctx.mode {
        Mode::Short => {
            write_per_base_quality(&mut w, &metrics)?;
            write_per_seq_quality(&mut w, &metrics)?;
            write_per_base_content(&mut w, &metrics)?;
            write_per_seq_gc(&mut w, &metrics)?;
            write_per_base_n(&mut w, &metrics)?;
            write_length_dist_short(&mut w, &metrics)?;
            write_duplication(&mut w, &metrics)?;
            write_overrep(&mut w, &metrics)?;
            write_adapter_content_short(&mut w, &metrics)?;
            #[cfg(not(feature = "no-kmer"))]
            write_kmer_content(&mut w, &metrics)?;
        }
        Mode::Long => {
            write_length_dist_long(&mut w, &metrics)?;
            write_per_seq_quality(&mut w, &metrics)?;
            write_per_seq_gc(&mut w, &metrics)?;
            write_per_seq_n(&mut w, &metrics)?;
            write_adapter_content_long(&mut w, &metrics)?;
        }
    }

    Ok(())
}

fn write_basic(
    w: &mut dyn Write,
    metrics: &crate::core::metrics::FinalMetrics,
    file_name: &str,
) -> Result<()> {
    writeln!(
        w,
        ">>Basic Statistics\t{}",
        metrics.statuses.basic.as_str_lower()
    )?;
    writeln!(w, "#Measure\tValue")?;
    writeln!(w, "Filename\t{}", file_name)?;
    writeln!(w, "File type\t{}", metrics.basic.file_type)?;
    writeln!(w, "Encoding\t{}", metrics.basic.encoding)?;
    writeln!(w, "Total Sequences\t{}", metrics.basic.total_sequences)?;
    writeln!(
        w,
        "Filtered Sequences\t{}",
        metrics.basic.filtered_sequences
    )?;
    if metrics.basic.min_len == metrics.basic.max_len {
        writeln!(w, "Sequence length\t{}", metrics.basic.min_len)?;
    } else {
        writeln!(
            w,
            "Sequence length\t{}-{}",
            metrics.basic.min_len, metrics.basic.max_len
        )?;
    }
    writeln!(w, "%GC\t{}", metrics.basic.gc_percent)?;
    writeln!(w, ">>END_MODULE")?;
    Ok(())
}

fn write_per_base_quality(
    w: &mut dyn Write,
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<()> {
    writeln!(
        w,
        ">>Per base sequence quality\t{}",
        metrics.statuses.per_base_qual.as_str_lower()
    )?;
    writeln!(
        w,
        "#Base\tMean\tMedian\tLower Quartile\tUpper Quartile\t10th Percentile\t90th Percentile"
    )?;
    for row in &metrics.per_base_qual {
        writeln!(
            w,
            "{}\t{:.1}\t{}\t{}\t{}\t{}\t{}",
            row.base,
            row.mean,
            row.median,
            row.lower_quartile,
            row.upper_quartile,
            row.p10,
            row.p90
        )?;
    }
    writeln!(w, ">>END_MODULE")?;
    Ok(())
}

fn write_per_seq_quality(
    w: &mut dyn Write,
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<()> {
    writeln!(
        w,
        ">>Per sequence quality scores\t{}",
        metrics.statuses.per_seq_qual.as_str_lower()
    )?;
    writeln!(w, "#Quality\tCount")?;
    for row in &metrics.per_seq_qual {
        writeln!(w, "{}\t{}", row.mean_q, row.count)?;
    }
    writeln!(w, ">>END_MODULE")?;
    Ok(())
}

fn write_per_base_content(
    w: &mut dyn Write,
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<()> {
    writeln!(
        w,
        ">>Per base sequence content\t{}",
        metrics.statuses.per_base_content.as_str_lower()
    )?;
    writeln!(w, "#Base\tG\tA\tT\tC")?;
    for row in &metrics.per_base_content {
        writeln!(
            w,
            "{}\t{:.1}\t{:.1}\t{:.1}\t{:.1}",
            row.base, row.g, row.a, row.t, row.c
        )?;
    }
    writeln!(w, ">>END_MODULE")?;
    Ok(())
}

fn write_per_seq_gc(w: &mut dyn Write, metrics: &crate::core::metrics::FinalMetrics) -> Result<()> {
    writeln!(
        w,
        ">>Per sequence GC content\t{}",
        metrics.statuses.per_seq_gc.as_str_lower()
    )?;
    writeln!(w, "#GC Content\tCount")?;
    for row in &metrics.per_seq_gc {
        writeln!(w, "{}\t{}", row.gc, row.count)?;
    }
    writeln!(w, ">>END_MODULE")?;
    Ok(())
}

fn write_per_base_n(w: &mut dyn Write, metrics: &crate::core::metrics::FinalMetrics) -> Result<()> {
    writeln!(
        w,
        ">>Per base N content\t{}",
        metrics.statuses.per_base_n.as_str_lower()
    )?;
    writeln!(w, "#Base\tN-Count")?;
    for row in &metrics.per_base_n {
        writeln!(w, "{}\t{:.1}", row.base, row.n_percent)?;
    }
    writeln!(w, ">>END_MODULE")?;
    Ok(())
}

fn write_per_seq_n(w: &mut dyn Write, metrics: &crate::core::metrics::FinalMetrics) -> Result<()> {
    writeln!(
        w,
        ">>Per sequence N content\t{}",
        metrics.statuses.per_seq_n.as_str_lower()
    )?;
    writeln!(w, "#N%\tCount")?;
    for row in &metrics.per_seq_n {
        writeln!(w, "{}\t{}", row.n_percent, row.count)?;
    }
    writeln!(w, ">>END_MODULE")?;
    Ok(())
}

fn write_length_dist_short(
    w: &mut dyn Write,
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<()> {
    writeln!(
        w,
        ">>Sequence Length Distribution\t{}",
        metrics.statuses.length_dist.as_str_lower()
    )?;
    writeln!(w, "#Length\tCount")?;
    for row in &metrics.length_dist {
        writeln!(w, "{}\t{}", row.length, row.count)?;
    }
    writeln!(w, ">>END_MODULE")?;
    Ok(())
}

fn write_length_dist_long(
    w: &mut dyn Write,
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<()> {
    writeln!(
        w,
        ">>Sequence Length Distribution\t{}",
        metrics.statuses.length_dist.as_str_lower()
    )?;
    if let Some(ref ll) = metrics.long_length {
        writeln!(w, "#Metric\tValue")?;
        writeln!(w, "Min\t{}", ll.min)?;
        writeln!(w, "Max\t{}", ll.max)?;
        writeln!(w, "Mean\t{:.1}", ll.mean)?;
        writeln!(w, "N50\t{}", ll.n50)?;
        writeln!(w, "N90\t{}", ll.n90)?;
        writeln!(w, "#Length\tCount")?;
        for i in 0..ll.bins.len() {
            writeln!(w, "{}\t{}", ll.labels[i], ll.bins[i])?;
        }
    }
    writeln!(w, ">>END_MODULE")?;
    Ok(())
}

fn write_duplication(
    w: &mut dyn Write,
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<()> {
    writeln!(
        w,
        ">>Sequence Duplication Levels\t{}",
        metrics.statuses.duplication.as_str_lower()
    )?;
    writeln!(w, "#Duplication Level\tRelative Count")?;
    for row in &metrics.duplication {
        writeln!(w, "{}\t{:.2}", row.level.as_str(), row.relative)?;
    }
    writeln!(w, ">>END_MODULE")?;
    Ok(())
}

fn write_overrep(w: &mut dyn Write, metrics: &crate::core::metrics::FinalMetrics) -> Result<()> {
    writeln!(
        w,
        ">>Overrepresented sequences\t{}",
        metrics.statuses.overrepresented.as_str_lower()
    )?;
    writeln!(w, "#Sequence\tCount\tPercentage\tPossible Source")?;
    for row in &metrics.overrepresented {
        writeln!(
            w,
            "{}\t{}\t{:.2}\t{}",
            row.sequence, row.count, row.percent, row.source
        )?;
    }
    writeln!(w, ">>END_MODULE")?;
    Ok(())
}

fn write_adapter_content_short(
    w: &mut dyn Write,
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<()> {
    writeln!(
        w,
        ">>Adapter Content\t{}",
        metrics.statuses.adapter_content.as_str_lower()
    )?;
    write!(w, "#Position")?;
    for name in crate::core::metrics::ADAPTERS {
        write!(w, "\t{}", name)?;
    }
    writeln!(w)?;
    for row in &metrics.adapter_content {
        write!(w, "{}", row.position)?;
        for v in row.values.iter() {
            write!(w, "\t{:.1}", v)?;
        }
        writeln!(w)?;
    }
    writeln!(w, ">>END_MODULE")?;
    Ok(())
}

fn write_adapter_content_long(
    w: &mut dyn Write,
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<()> {
    writeln!(
        w,
        ">>Adapter Content\t{}",
        metrics.statuses.adapter_content.as_str_lower()
    )?;
    write!(w, "#Adapter")?;
    for name in crate::core::metrics::ADAPTERS {
        write!(w, "\t{}", name)?;
    }
    writeln!(w)?;
    if let Some(row) = metrics.adapter_content.first() {
        write!(w, "Any")?;
        for v in row.values.iter() {
            write!(w, "\t{:.1}", v)?;
        }
        writeln!(w)?;
    }
    writeln!(w, ">>END_MODULE")?;
    Ok(())
}

#[cfg(not(feature = "no-kmer"))]
fn write_kmer_content(
    w: &mut dyn Write,
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<()> {
    writeln!(
        w,
        ">>Kmer Content\t{}",
        metrics.statuses.kmer_content.as_str_lower()
    )?;
    writeln!(
        w,
        "#Sequence\tCount\tPValue\tObs/Exp Max\tMax Obs/Exp Position"
    )?;
    for row in &metrics.kmer_rows {
        writeln!(
            w,
            "{}\t{}\t{:.2e}\t{:.2}\t{}",
            row.sequence, row.count, row.p_value, row.obs_exp, row.max_pos
        )?;
    }
    writeln!(w, ">>END_MODULE")?;
    Ok(())
}
