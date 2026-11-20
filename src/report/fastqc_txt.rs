use crate::core::metrics::{BasicStats, FinalMetrics};
use crate::core::model::FinalizeContext;
use crate::report::{ModuleId, active_modules};
use crate::{FASTQC_FORMAT_VERSION, KIRA_QC_VERSION};
use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

pub fn write(path: &Path, metrics: &FinalMetrics, ctx: &FinalizeContext) -> Result<()> {
    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let mut w = BufWriter::new(file);
    write_to(&mut w, metrics, ctx)
}

pub fn write_to(w: &mut dyn Write, metrics: &FinalMetrics, ctx: &FinalizeContext) -> Result<()> {
    // MultiQC keys off this header line to identify the file and its format.
    writeln!(w, "##FastQC\t{}", FASTQC_FORMAT_VERSION)?;
    writeln!(w, "##kira-qc\t{}", KIRA_QC_VERSION)?;

    for module in active_modules(metrics, ctx) {
        writeln!(
            w,
            ">>{}\t{}",
            module.title(),
            module.status(metrics).as_str_lower()
        )?;
        write_body(w, module, metrics, ctx)?;
        writeln!(w, ">>END_MODULE")?;
    }
    Ok(())
}

fn write_body(
    w: &mut dyn Write,
    module: ModuleId,
    metrics: &FinalMetrics,
    ctx: &FinalizeContext,
) -> Result<()> {
    match module {
        ModuleId::Basic => write_basic(w, metrics, ctx),
        ModuleId::PerBaseQual => write_per_base_quality(w, metrics),
        ModuleId::PerTile => write_per_tile(w, metrics),
        ModuleId::PerSeqQual => write_per_seq_quality(w, metrics),
        ModuleId::PerBaseContent => write_per_base_content(w, metrics),
        ModuleId::PerSeqGc => write_per_seq_gc(w, metrics),
        ModuleId::PerBaseN => write_per_base_n(w, metrics),
        ModuleId::PerSeqN => write_per_seq_n(w, metrics),
        ModuleId::LengthDist => write_length_dist(w, metrics, ctx),
        ModuleId::Duplication => write_duplication(w, metrics, ctx),
        ModuleId::Overrepresented => write_overrep(w, metrics),
        ModuleId::AdapterContent => write_adapter_content(w, metrics, ctx),
        ModuleId::KmerContent => write_kmer_content(w, metrics),
    }
}

fn write_basic(w: &mut dyn Write, metrics: &FinalMetrics, ctx: &FinalizeContext) -> Result<()> {
    let b = &metrics.basic;
    writeln!(w, "#Measure\tValue")?;
    writeln!(w, "Filename\t{}", ctx.file_name)?;
    writeln!(w, "File type\t{}", b.file_type)?;
    writeln!(w, "Encoding\t{}", b.encoding)?;
    writeln!(w, "Total Sequences\t{}", b.total_sequences)?;
    // Added in FastQC 0.12 and parsed by MultiQC.
    writeln!(
        w,
        "Total Bases\t{}",
        BasicStats::format_bases(b.total_bases)
    )?;
    writeln!(
        w,
        "Sequences flagged as poor quality\t{}",
        b.poor_quality_sequences
    )?;
    writeln!(w, "Sequence length\t{}", b.length_range())?;
    writeln!(w, "%GC\t{}", b.gc_percent)?;
    Ok(())
}

fn write_per_base_quality(w: &mut dyn Write, metrics: &FinalMetrics) -> Result<()> {
    writeln!(
        w,
        "#Base\tMean\tMedian\tLower Quartile\tUpper Quartile\t10th Percentile\t90th Percentile"
    )?;
    for row in &metrics.per_base_qual {
        writeln!(
            w,
            "{}\t{:.1}\t{}\t{}\t{}\t{}\t{}",
            row.group.label(),
            row.mean,
            row.median,
            row.lower_quartile,
            row.upper_quartile,
            row.p10,
            row.p90
        )?;
    }
    Ok(())
}

fn write_per_tile(w: &mut dyn Write, metrics: &FinalMetrics) -> Result<()> {
    writeln!(w, "#Tile\tBase\tMean")?;
    for row in &metrics.per_tile {
        writeln!(
            w,
            "{}\t{}\t{:.5}",
            row.tile,
            row.group.label(),
            row.mean_deviation
        )?;
    }
    Ok(())
}

fn write_per_seq_quality(w: &mut dyn Write, metrics: &FinalMetrics) -> Result<()> {
    writeln!(w, "#Quality\tCount")?;
    for row in &metrics.per_seq_qual {
        writeln!(w, "{}\t{}", row.mean_q, row.count)?;
    }
    Ok(())
}

fn write_per_base_content(w: &mut dyn Write, metrics: &FinalMetrics) -> Result<()> {
    writeln!(w, "#Base\tG\tA\tT\tC")?;
    for row in &metrics.per_base_content {
        writeln!(
            w,
            "{}\t{:.1}\t{:.1}\t{:.1}\t{:.1}",
            row.group.label(),
            row.g,
            row.a,
            row.t,
            row.c
        )?;
    }
    Ok(())
}

fn write_per_seq_gc(w: &mut dyn Write, metrics: &FinalMetrics) -> Result<()> {
    writeln!(w, "#GC Content\tCount")?;
    for row in &metrics.per_seq_gc {
        writeln!(w, "{}\t{}", row.gc, row.count)?;
    }
    Ok(())
}

fn write_per_base_n(w: &mut dyn Write, metrics: &FinalMetrics) -> Result<()> {
    writeln!(w, "#Base\tN-Count")?;
    for row in &metrics.per_base_n {
        writeln!(w, "{}\t{:.1}", row.group.label(), row.n_percent)?;
    }
    Ok(())
}

fn write_per_seq_n(w: &mut dyn Write, metrics: &FinalMetrics) -> Result<()> {
    writeln!(w, "#N%\tCount")?;
    for row in &metrics.per_seq_n {
        writeln!(w, "{}\t{}", row.n_percent, row.count)?;
    }
    Ok(())
}

fn write_length_dist(
    w: &mut dyn Write,
    metrics: &FinalMetrics,
    _ctx: &FinalizeContext,
) -> Result<()> {
    if let Some(ref ll) = metrics.long_length {
        // Long-read summary: the headline numbers first, then the log-scaled
        // histogram. N50/N90 are exact, not bin midpoints.
        writeln!(w, "#Metric\tValue")?;
        writeln!(w, "Min\t{}", ll.min)?;
        writeln!(w, "Max\t{}", ll.max)?;
        writeln!(w, "Mean\t{:.1}", ll.mean)?;
        writeln!(w, "Median\t{}", ll.median)?;
        writeln!(w, "N50\t{}", ll.n50)?;
        writeln!(w, "N90\t{}", ll.n90)?;
        writeln!(w, "Total Bases\t{}", ll.total_bases)?;
        writeln!(w, "#Length\tCount")?;
        for (label, count) in &ll.bins {
            writeln!(w, "{}\t{}", label, count)?;
        }
    } else {
        writeln!(w, "#Length\tCount")?;
        for row in &metrics.length_dist {
            writeln!(w, "{}\t{}", row.length, row.count)?;
        }
    }
    Ok(())
}

fn write_duplication(
    w: &mut dyn Write,
    metrics: &FinalMetrics,
    ctx: &FinalizeContext,
) -> Result<()> {
    let Some(ref summary) = metrics.duplication_summary else {
        writeln!(
            w,
            "#Duplication Level\tPercentage of deduplicated\tPercentage of total"
        )?;
        return Ok(());
    };

    // MultiQC reads this exact line for the "% Dups" column in General Stats,
    // and the three-column table below for the duplication plot.
    writeln!(
        w,
        "#Total Deduplicated Percentage\t{:.4}",
        summary.total_deduplicated_percent
    )?;
    if !ctx.dup_window_complete {
        writeln!(
            w,
            "#Estimated from the first {} reads",
            ctx.dup_window_reads
        )?;
    }
    writeln!(
        w,
        "#Duplication Level\tPercentage of deduplicated\tPercentage of total"
    )?;
    for row in &metrics.duplication {
        writeln!(
            w,
            "{}\t{:.6}\t{:.6}",
            row.level, row.percent_of_deduplicated, row.percent_of_total
        )?;
    }
    Ok(())
}

fn write_overrep(w: &mut dyn Write, metrics: &FinalMetrics) -> Result<()> {
    writeln!(w, "#Sequence\tCount\tPercentage\tPossible Source")?;
    for row in &metrics.overrepresented {
        writeln!(
            w,
            "{}\t{}\t{}\t{}",
            row.sequence, row.count, row.percent, row.source
        )?;
    }
    Ok(())
}

fn write_adapter_content(
    w: &mut dyn Write,
    metrics: &FinalMetrics,
    ctx: &FinalizeContext,
) -> Result<()> {
    let is_long = metrics.long_length.is_some() || ctx.mode == crate::core::model::Mode::Long;
    // Column headers are adapter *names*, as FastQC writes them; the previous
    // implementation emitted raw sequences, which no parser expects.
    write!(w, "#{}", if is_long { "Adapter" } else { "Position" })?;
    for name in &metrics.adapter_names {
        write!(w, "\t{}", name)?;
    }
    writeln!(w)?;

    if is_long {
        if let Some(row) = metrics.adapter_content.first() {
            write!(w, "Any")?;
            for v in &row.values {
                write!(w, "\t{:.1}", v)?;
            }
            writeln!(w)?;
        }
        return Ok(());
    }

    for row in &metrics.adapter_content {
        write!(w, "{}", row.group.label())?;
        for v in &row.values {
            write!(w, "\t{:.1}", v)?;
        }
        writeln!(w)?;
    }
    Ok(())
}

fn write_kmer_content(w: &mut dyn Write, metrics: &FinalMetrics) -> Result<()> {
    writeln!(
        w,
        "#Sequence\tCount\tPValue\tObs/Exp Max\tMax Obs/Exp Position"
    )?;
    for row in &metrics.kmer_rows {
        writeln!(
            w,
            "{}\t{}\t{:.3e}\t{:.2}\t{}",
            row.sequence, row.count, row.p_value, row.obs_exp, row.max_pos
        )?;
    }
    Ok(())
}
