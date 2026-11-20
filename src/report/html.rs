use crate::core::metrics::FinalMetrics;
use crate::core::model::{FinalizeContext, Mode, Status};
use crate::report::{ModuleId, active_modules, fmt_int, format_timestamp, report_timestamp};
use crate::{KIRA_QC_VERSION, core::metrics::BasicStats};
use anyhow::{Context, Result};
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

const PLOT_W: f64 = 800.0;
const PLOT_H: f64 = 260.0;

/// Escape text that ends up inside HTML or SVG markup. Filenames, sample names,
/// contaminant labels and read sequences all reach the report unmodified, so
/// none of them may be interpolated raw.
pub fn esc(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&#39;"),
            c => out.push(c),
        }
    }
    out
}

pub fn write(path: &Path, metrics: &FinalMetrics, ctx: &FinalizeContext) -> Result<()> {
    let html = render(metrics, ctx)?;
    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let mut w = BufWriter::new(file);
    w.write_all(html.as_bytes())?;
    Ok(())
}

pub fn render(metrics: &FinalMetrics, ctx: &FinalizeContext) -> Result<String> {
    let mut html = String::with_capacity(256 * 1024);
    let ts = report_timestamp();
    let modules = active_modules(metrics, ctx);

    writeln!(html, "<!DOCTYPE html>")?;
    writeln!(html, "<html lang=\"en\">")?;
    writeln!(html, "<head>")?;
    writeln!(html, "<meta charset=\"utf-8\"/>")?;
    writeln!(
        html,
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\"/>"
    )?;
    writeln!(
        html,
        "<title>kira-qc report: {}</title>",
        esc(&ctx.sample_name)
    )?;
    write_style(&mut html)?;
    writeln!(html, "</head>")?;
    writeln!(html, "<body>")?;

    let mode_label = match ctx.mode {
        Mode::Short => "Short-read (Illumina)",
        Mode::Long => "Long-read (ONT / PacBio)",
    };

    writeln!(html, "<div class=\"page\">")?;
    writeln!(html, "<aside class=\"sidebar\">")?;
    writeln!(html, "<h2 id=\"summary\">Summary</h2>")?;
    writeln!(html, "<ul>")?;
    for m in &modules {
        writeln!(
            html,
            "<li>{} <a href=\"#{}\">{}</a></li>",
            status_icon_svg(m.status(metrics), 14),
            m.anchor(),
            esc(m.title())
        )?;
    }
    writeln!(html, "</ul>")?;
    writeln!(html, "</aside>")?;

    writeln!(html, "<main class=\"main\">")?;
    writeln!(html, "<h1>kira-qc report</h1>")?;
    writeln!(
        html,
        "<div class=\"meta\">File: {}<br/>Mode: {}<br/>kira-qc {}<br/>Generated: {}</div>",
        esc(&ctx.file_name),
        mode_label,
        esc(KIRA_QC_VERSION),
        esc(&format_timestamp(ts))
    )?;

    for m in &modules {
        section_header(&mut html, m.status(metrics), m.title(), m.anchor())?;
        section_body(&mut html, *m, metrics, ctx)?;
        section_footer(&mut html)?;
    }

    writeln!(html, "<div class=\"meta\">Produced by kira-qc</div>")?;
    writeln!(html, "</main>")?;
    writeln!(html, "</div>")?;
    writeln!(html, "</body></html>")?;
    Ok(html)
}

fn write_style(html: &mut String) -> Result<()> {
    writeln!(html, "<style>")?;
    writeln!(
        html,
        "body{{font-family:Arial,Helvetica,sans-serif;margin:0;background:#eee;color:#222;}}"
    )?;
    writeln!(
        html,
        ".page{{display:flex;align-items:flex-start;gap:16px;padding:16px;}}"
    )?;
    writeln!(
        html,
        ".sidebar{{width:260px;position:sticky;top:16px;align-self:flex-start;background:#f6f6f6;border:1px solid #ddd;border-radius:4px;padding:10px;}}"
    )?;
    writeln!(html, ".sidebar h2{{margin:4px 0 8px 0;font-size:16px;}}")?;
    writeln!(html, ".sidebar ul{{list-style:none;margin:0;padding:0;}}")?;
    writeln!(
        html,
        ".sidebar li{{display:flex;align-items:center;gap:8px;padding:4px 0;font-size:13px;}}"
    )?;
    writeln!(html, ".sidebar a{{color:#003366;text-decoration:none;}}")?;
    writeln!(html, ".sidebar a:hover{{text-decoration:underline;}}")?;
    writeln!(
        html,
        ".main{{flex:1;min-width:0;background:#fff;border:1px solid #ddd;border-radius:4px;box-shadow:0 1px 3px rgba(0,0,0,0.08);padding:16px 20px;}}"
    )?;
    writeln!(html, "h1{{margin:0 0 6px 0;font-size:22px;}}")?;
    writeln!(html, "h2{{margin:20px 0 6px 0;font-size:18px;}}")?;
    writeln!(
        html,
        ".meta{{color:#555;font-size:12px;margin-bottom:12px;}}"
    )?;
    writeln!(
        html,
        ".module{{padding:8px 0 14px 0;border-bottom:1px solid #eee;}}"
    )?;
    writeln!(html, ".module:last-child{{border-bottom:none;}}")?;
    writeln!(
        html,
        ".module h2{{display:flex;align-items:center;gap:8px;}}"
    )?;
    writeln!(html, ".plot{{margin:8px 0 6px 0;overflow-x:auto;}}")?;
    writeln!(
        html,
        ".desc{{color:#444;font-size:13px;max-width:1000px;margin:4px 0 10px 0;}}"
    )?;
    writeln!(
        html,
        ".note{{color:#666;font-size:12px;font-style:italic;margin:4px 0 8px 0;}}"
    )?;
    writeln!(html, ".tablewrap{{overflow-x:auto;max-height:520px;}}")?;
    writeln!(
        html,
        ".table{{border-collapse:collapse;width:100%;max-width:1000px;font-size:12px;}}"
    )?;
    writeln!(
        html,
        ".table th,.table td{{border:1px solid #ddd;padding:4px 6px;text-align:right;}}"
    )?;
    writeln!(
        html,
        ".table th:first-child,.table td:first-child{{text-align:left;}}"
    )?;
    writeln!(
        html,
        ".table td.seq{{font-family:monospace;text-align:left;}}"
    )?;
    writeln!(
        html,
        ".bs-table{{border-collapse:collapse;font-size:12px;width:420px;}}"
    )?;
    writeln!(
        html,
        ".bs-table th{{background:#3b6ea5;color:#fff;text-align:left;padding:4px 6px;border:1px solid #2f5a86;}}"
    )?;
    writeln!(
        html,
        ".bs-table td{{border:1px solid #ddd;padding:4px 6px;text-align:left;}}"
    )?;
    writeln!(html, "details{{margin:6px 0 0 0;}}")?;
    writeln!(
        html,
        ".back{{font-size:12px;margin-top:6px;display:inline-block;}}"
    )?;
    writeln!(
        html,
        "section:target{{outline:2px solid #99c;outline-offset:4px;border-radius:4px;}}"
    )?;
    writeln!(html, "svg{{background:#fafafa;border:1px solid #e5e5e5;}}")?;
    writeln!(html, "</style>")?;
    Ok(())
}

fn section_body(
    out: &mut String,
    module: ModuleId,
    metrics: &FinalMetrics,
    ctx: &FinalizeContext,
) -> Result<()> {
    match module {
        ModuleId::Basic => basic_stats(out, metrics, ctx),
        ModuleId::PerBaseQual => per_base_quality(out, metrics),
        ModuleId::PerTile => per_tile(out, metrics),
        ModuleId::PerSeqQual => per_seq_quality(out, metrics, ctx),
        ModuleId::PerBaseContent => per_base_content(out, metrics),
        ModuleId::PerSeqGc => per_seq_gc(out, metrics),
        ModuleId::PerBaseN => per_base_n(out, metrics),
        ModuleId::PerSeqN => per_seq_n(out, metrics),
        ModuleId::LengthDist => length_dist(out, metrics),
        ModuleId::Duplication => duplication(out, metrics, ctx),
        ModuleId::Overrepresented => overrep(out, metrics, ctx),
        ModuleId::AdapterContent => adapter_content(out, metrics, ctx),
        ModuleId::KmerContent => kmer_content(out, metrics),
    }
}

fn module_desc(out: &mut String, text: &str) -> Result<()> {
    writeln!(out, "<p class=\"desc\">{}</p>", esc(text))?;
    Ok(())
}

fn note(out: &mut String, text: &str) -> Result<()> {
    writeln!(out, "<p class=\"note\">{}</p>", esc(text))?;
    Ok(())
}

fn status_icon_svg(status: Status, size: u32) -> String {
    let (fill, mark) = match status {
        Status::Pass => ("#2e8b57", "M6 10 L10 14 L18 6"),
        Status::Warn => ("#e6a400", "M11 5 L11 13 M11 16 L11 18"),
        Status::Fail => ("#c00000", "M6 6 L18 18 M18 6 L6 18"),
    };
    format!(
        "<svg width=\"{s}\" height=\"{s}\" viewBox=\"0 0 24 24\" role=\"img\" aria-label=\"{a}\"><circle cx=\"12\" cy=\"12\" r=\"11\" fill=\"{f}\"/><path d=\"{p}\" stroke=\"#fff\" stroke-width=\"2\" fill=\"none\" stroke-linecap=\"round\" stroke-linejoin=\"round\"/></svg>",
        s = size,
        a = status.as_str_upper(),
        f = fill,
        p = mark
    )
}

fn section_header(out: &mut String, status: Status, title: &str, id: &str) -> Result<()> {
    writeln!(out, "<section id=\"{}\" class=\"module\">", id)?;
    writeln!(
        out,
        "<h2>{} {}</h2>",
        status_icon_svg(status, 16),
        esc(title)
    )?;
    Ok(())
}

fn section_footer(out: &mut String) -> Result<()> {
    writeln!(
        out,
        "<a class=\"back\" href=\"#summary\">Back to Summary</a>"
    )?;
    writeln!(out, "</section>")?;
    Ok(())
}

/// Collapsible data table.
fn table<F>(out: &mut String, summary: &str, f: F) -> Result<()>
where
    F: FnOnce(&mut String) -> Result<()>,
{
    writeln!(
        out,
        "<details><summary>{}</summary><div class=\"tablewrap\"><table class=\"table\">",
        esc(summary)
    )?;
    f(out)?;
    writeln!(out, "</table></div></details>")?;
    Ok(())
}

fn basic_stats(out: &mut String, metrics: &FinalMetrics, ctx: &FinalizeContext) -> Result<()> {
    module_desc(
        out,
        "Summary of input size, read counts, length range and GC%. A wide length range or unusual GC can indicate mixed libraries or contamination.",
    )?;
    let b = &metrics.basic;
    writeln!(out, "<table class=\"bs-table\">")?;
    writeln!(out, "<tr><th>Measure</th><th>Value</th></tr>")?;
    writeln!(
        out,
        "<tr><td>Filename</td><td>{}</td></tr>",
        esc(&ctx.file_name)
    )?;
    writeln!(
        out,
        "<tr><td>File type</td><td>{}</td></tr>",
        esc(b.file_type)
    )?;
    writeln!(
        out,
        "<tr><td>Encoding</td><td>{}</td></tr>",
        esc(b.encoding)
    )?;
    writeln!(
        out,
        "<tr><td>Total Sequences</td><td>{}</td></tr>",
        fmt_int(b.total_sequences)
    )?;
    writeln!(
        out,
        "<tr><td>Total Bases</td><td>{}</td></tr>",
        esc(&BasicStats::format_bases(b.total_bases))
    )?;
    writeln!(
        out,
        "<tr><td>Sequences flagged as poor quality</td><td>{}</td></tr>",
        fmt_int(b.poor_quality_sequences)
    )?;
    writeln!(
        out,
        "<tr><td>Sequence length</td><td>{}</td></tr>",
        esc(&b.length_range())
    )?;
    writeln!(out, "<tr><td>%GC</td><td>{}</td></tr>", b.gc_percent)?;
    writeln!(out, "</table>")?;
    Ok(())
}

fn per_base_quality(out: &mut String, metrics: &FinalMetrics) -> Result<()> {
    module_desc(
        out,
        "Quality score distribution at each base position. Systematic drops toward read ends usually reflect sequencing chemistry degradation or adapter read-through.",
    )?;
    let max_q = metrics
        .per_base_qual
        .iter()
        .map(|r| r.p90 as f64)
        .fold(40.0, f64::max);
    svg_boxplot(out, &metrics.per_base_qual, PLOT_W, PLOT_H, max_q)?;
    table(out, "Data", |o| {
        writeln!(
            o,
            "<tr><th>Base</th><th>Mean</th><th>Median</th><th>Lower Quartile</th><th>Upper Quartile</th><th>10th Percentile</th><th>90th Percentile</th></tr>"
        )?;
        for r in &metrics.per_base_qual {
            writeln!(
                o,
                "<tr><td>{}</td><td>{:.1}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>",
                esc(&r.group.label()),
                r.mean,
                r.median,
                r.lower_quartile,
                r.upper_quartile,
                r.p10,
                r.p90
            )?;
        }
        Ok(())
    })
}

fn per_tile(out: &mut String, metrics: &FinalMetrics) -> Result<()> {
    module_desc(
        out,
        "Deviation of each flowcell tile from the average quality at that position. Cold spots point at a physical problem on the flowcell rather than at the library.",
    )?;
    svg_tile_heatmap(out, metrics, PLOT_W, PLOT_H)?;
    table(out, "Data", |o| {
        writeln!(
            o,
            "<tr><th>Tile</th><th>Base</th><th>Mean deviation</th></tr>"
        )?;
        for r in &metrics.per_tile {
            writeln!(
                o,
                "<tr><td>{}</td><td>{}</td><td>{:.3}</td></tr>",
                r.tile,
                esc(&r.group.label()),
                r.mean_deviation
            )?;
        }
        Ok(())
    })
}

fn per_seq_quality(out: &mut String, metrics: &FinalMetrics, ctx: &FinalizeContext) -> Result<()> {
    module_desc(
        out,
        "Distribution of mean quality per read. A left-shifted distribution indicates a globally low-quality run or a mixture of good and bad reads.",
    )?;
    if ctx.mode == Mode::Long {
        note(
            out,
            "Long-read mode averages error probabilities before converting back to a Phred score, which is the convention used by NanoPlot and related tools.",
        )?;
    }
    let data: Vec<(f64, f64)> = metrics
        .per_seq_qual
        .iter()
        .map(|r| (r.mean_q as f64, r.count as f64))
        .collect();
    svg_bars(out, &data, PLOT_W, PLOT_H, "Mean Q", "Count")?;
    table(out, "Data", |o| {
        writeln!(o, "<tr><th>Quality</th><th>Count</th></tr>")?;
        for r in &metrics.per_seq_qual {
            writeln!(
                o,
                "<tr><td>{}</td><td>{}</td></tr>",
                r.mean_q,
                fmt_int(r.count)
            )?;
        }
        Ok(())
    })
}

fn per_base_content(out: &mut String, metrics: &FinalMetrics) -> Result<()> {
    module_desc(
        out,
        "Percentage of each base at each position. The module flags |A-T| and |G-C| imbalance, so a GC-rich genome is not reported as bias.",
    )?;
    legend_base_content(out)?;
    svg_base_content(out, &metrics.per_base_content, PLOT_W, PLOT_H)?;
    table(out, "Data", |o| {
        writeln!(
            o,
            "<tr><th>Base</th><th>G</th><th>A</th><th>T</th><th>C</th></tr>"
        )?;
        for r in &metrics.per_base_content {
            writeln!(
                o,
                "<tr><td>{}</td><td>{:.1}</td><td>{:.1}</td><td>{:.1}</td><td>{:.1}</td></tr>",
                esc(&r.group.label()),
                r.g,
                r.a,
                r.t,
                r.c
            )?;
        }
        Ok(())
    })
}

fn per_seq_gc(out: &mut String, metrics: &FinalMetrics) -> Result<()> {
    module_desc(
        out,
        "Distribution of GC% across reads, compared against a fitted normal distribution. Broad or multi-modal shapes suggest contamination or mixed libraries.",
    )?;
    let data: Vec<(f64, f64)> = metrics
        .per_seq_gc
        .iter()
        .map(|r| (r.gc as f64, r.count as f64))
        .collect();
    svg_bars_xbands(
        out,
        &data,
        PLOT_W,
        PLOT_H,
        0.0,
        100.0,
        &[(40.0, 60.0, "#cdeccf")],
        "GC%",
        "Count",
    )?;
    table(out, "Data", |o| {
        writeln!(o, "<tr><th>GC%</th><th>Count</th></tr>")?;
        for r in &metrics.per_seq_gc {
            writeln!(o, "<tr><td>{}</td><td>{}</td></tr>", r.gc, fmt_int(r.count))?;
        }
        Ok(())
    })
}

fn per_base_n(out: &mut String, metrics: &FinalMetrics) -> Result<()> {
    module_desc(
        out,
        "Proportion of Ns at each position. Spikes suggest base-calling failures or low-complexity regions.",
    )?;
    let data: Vec<(f64, f64)> = metrics
        .per_base_n
        .iter()
        .map(|r| (r.group.mid(), r.n_percent))
        .collect();
    let (y_min, y_max) = auto_range(data.iter().map(|(_, y)| *y), 0.0, 100.0);
    svg_line_ybands(
        out,
        &data,
        PLOT_W,
        PLOT_H,
        y_min,
        y_max,
        "#555",
        &[
            (0.0, 5.0, "#cdeccf"),
            (5.0, 20.0, "#ffe5b4"),
            (20.0, 100.0, "#f4c7c3"),
        ],
        "Position",
        "% N",
    )?;
    table(out, "Data", |o| {
        writeln!(o, "<tr><th>Base</th><th>N%</th></tr>")?;
        for r in &metrics.per_base_n {
            writeln!(
                o,
                "<tr><td>{}</td><td>{:.1}</td></tr>",
                esc(&r.group.label()),
                r.n_percent
            )?;
        }
        Ok(())
    })
}

fn per_seq_n(out: &mut String, metrics: &FinalMetrics) -> Result<()> {
    module_desc(
        out,
        "Distribution of N% per read. Excess high-N reads indicate poor base-calling or unresolved segments.",
    )?;
    let data: Vec<(f64, f64)> = metrics
        .per_seq_n
        .iter()
        .map(|r| (r.n_percent as f64, r.count as f64))
        .collect();
    svg_bars_xbands(
        out,
        &data,
        PLOT_W,
        PLOT_H,
        0.0,
        100.0,
        &[
            (0.0, 10.0, "#cdeccf"),
            (10.0, 20.0, "#ffe5b4"),
            (20.0, 100.0, "#f4c7c3"),
        ],
        "N%",
        "Count",
    )?;
    table(out, "Data", |o| {
        writeln!(o, "<tr><th>N%</th><th>Count</th></tr>")?;
        for r in &metrics.per_seq_n {
            writeln!(
                o,
                "<tr><td>{}</td><td>{}</td></tr>",
                r.n_percent,
                fmt_int(r.count)
            )?;
        }
        Ok(())
    })
}

fn length_dist(out: &mut String, metrics: &FinalMetrics) -> Result<()> {
    if let Some(ref ll) = metrics.long_length {
        module_desc(
            out,
            "Read length frequencies in log-scaled bins. N50 and N90 are computed from the full length distribution, not from bin midpoints.",
        )?;
        let labels: Vec<String> = ll.bins.iter().map(|(l, _)| l.clone()).collect();
        let values: Vec<f64> = ll.bins.iter().map(|(_, c)| *c as f64).collect();
        svg_category_bars(out, &labels, &values, PLOT_W, PLOT_H, "Length bin", "Count")?;
        table(out, "Data", |o| {
            writeln!(o, "<tr><th>Metric</th><th>Value</th></tr>")?;
            writeln!(o, "<tr><td>Min</td><td>{}</td></tr>", fmt_int(ll.min))?;
            writeln!(o, "<tr><td>Max</td><td>{}</td></tr>", fmt_int(ll.max))?;
            writeln!(o, "<tr><td>Mean</td><td>{:.1}</td></tr>", ll.mean)?;
            writeln!(o, "<tr><td>Median</td><td>{}</td></tr>", fmt_int(ll.median))?;
            writeln!(o, "<tr><td>N50</td><td>{}</td></tr>", fmt_int(ll.n50))?;
            writeln!(o, "<tr><td>N90</td><td>{}</td></tr>", fmt_int(ll.n90))?;
            writeln!(
                o,
                "<tr><td>Total bases</td><td>{}</td></tr>",
                fmt_int(ll.total_bases)
            )?;
            for (label, count) in &ll.bins {
                writeln!(
                    o,
                    "<tr><td>{}</td><td>{}</td></tr>",
                    esc(label),
                    fmt_int(*count)
                )?;
            }
            Ok(())
        })
    } else {
        module_desc(
            out,
            "Read length frequencies. Multiple peaks or long tails may indicate trimming or mixed read sources.",
        )?;
        let data: Vec<(f64, f64)> = metrics
            .length_dist
            .iter()
            .map(|r| (r.length as f64, r.count as f64))
            .collect();
        svg_bars(out, &data, PLOT_W, PLOT_H, "Length", "Count")?;
        table(out, "Data", |o| {
            writeln!(o, "<tr><th>Length</th><th>Count</th></tr>")?;
            for r in &metrics.length_dist {
                writeln!(
                    o,
                    "<tr><td>{}</td><td>{}</td></tr>",
                    r.length,
                    fmt_int(r.count)
                )?;
            }
            Ok(())
        })
    }
}

fn duplication(out: &mut String, metrics: &FinalMetrics, ctx: &FinalizeContext) -> Result<()> {
    module_desc(
        out,
        "Duplication levels counted exactly over a window of reads. High duplication usually means PCR over-amplification or low library complexity.",
    )?;
    if let Some(ref s) = metrics.duplication_summary {
        writeln!(
            out,
            "<p class=\"desc\"><b>Total deduplicated percentage:</b> {:.2}%</p>",
            s.total_deduplicated_percent
        )?;
        if !ctx.dup_window_complete {
            note(
                out,
                &format!(
                    "Estimated from the first {} reads of the file.",
                    fmt_int(ctx.dup_window_reads)
                ),
            )?;
        }
    }

    let labels: Vec<String> = metrics
        .duplication
        .iter()
        .map(|r| r.level.to_string())
        .collect();
    let dedup: Vec<f64> = metrics
        .duplication
        .iter()
        .map(|r| r.percent_of_deduplicated)
        .collect();
    let total: Vec<f64> = metrics
        .duplication
        .iter()
        .map(|r| r.percent_of_total)
        .collect();
    svg_category_lines(
        out,
        &labels,
        &[
            ("% Deduplicated", &dedup, "#1f77b4"),
            ("% Total", &total, "#d62728"),
        ],
        PLOT_W,
        PLOT_H,
        "Duplication level",
        "%",
    )?;

    table(out, "Data", |o| {
        writeln!(
            o,
            "<tr><th>Duplication Level</th><th>Percentage of deduplicated</th><th>Percentage of total</th></tr>"
        )?;
        for r in &metrics.duplication {
            writeln!(
                o,
                "<tr><td>{}</td><td>{:.2}</td><td>{:.2}</td></tr>",
                esc(r.level),
                r.percent_of_deduplicated,
                r.percent_of_total
            )?;
        }
        Ok(())
    })
}

fn overrep(out: &mut String, metrics: &FinalMetrics, ctx: &FinalizeContext) -> Result<()> {
    module_desc(
        out,
        "Sequences occurring far more often than expected. Common sources are adapters, primers and contamination.",
    )?;
    if !ctx.dup_window_complete {
        note(
            out,
            &format!(
                "Counted over the first {} reads of the file.",
                fmt_int(ctx.dup_window_reads)
            ),
        )?;
    }
    if metrics.overrepresented.is_empty() {
        writeln!(
            out,
            "<p class=\"desc\">No sequences accounted for more than 0.1% of reads.</p>"
        )?;
        return Ok(());
    }
    table(out, "Data", |o| {
        writeln!(
            o,
            "<tr><th>Sequence</th><th>Count</th><th>Percentage</th><th>Possible Source</th></tr>"
        )?;
        for r in &metrics.overrepresented {
            writeln!(
                o,
                "<tr><td class=\"seq\">{}</td><td>{}</td><td>{:.4}</td><td>{}</td></tr>",
                esc(&r.sequence),
                fmt_int(r.count),
                r.percent,
                esc(&r.source)
            )?;
        }
        Ok(())
    })
}

fn adapter_content(out: &mut String, metrics: &FinalMetrics, ctx: &FinalizeContext) -> Result<()> {
    if ctx.mode == Mode::Long {
        module_desc(
            out,
            "Fraction of reads containing a known adapter motif, searched on both strands. Elevated values suggest residual adapters or chimeric reads.",
        )?;
        table(out, "Data", |o| {
            writeln!(o, "<tr><th>Adapter</th><th>% of reads</th></tr>")?;
            if let Some(row) = metrics.adapter_content.first() {
                for (i, name) in metrics.adapter_names.iter().enumerate() {
                    writeln!(
                        o,
                        "<tr><td>{}</td><td>{:.2}</td></tr>",
                        esc(name),
                        row.values.get(i).copied().unwrap_or(0.0)
                    )?;
                }
            }
            Ok(())
        })?;
        return Ok(());
    }

    module_desc(
        out,
        "Cumulative adapter match rate by position, using 12 bp adapter prefixes as FastQC does. A rising curve toward read ends means adapter read-through.",
    )?;
    svg_adapter_lines(out, metrics, PLOT_W, PLOT_H)?;
    table(out, "Data", |o| {
        write!(o, "<tr><th>Position</th>")?;
        for name in &metrics.adapter_names {
            write!(o, "<th>{}</th>", esc(name))?;
        }
        writeln!(o, "</tr>")?;
        for r in &metrics.adapter_content {
            write!(o, "<tr><td>{}</td>", esc(&r.group.label()))?;
            for v in &r.values {
                write!(o, "<td>{:.2}</td>", v)?;
            }
            writeln!(o, "</tr>")?;
        }
        Ok(())
    })
}

fn kmer_content(out: &mut String, metrics: &FinalMetrics) -> Result<()> {
    module_desc(
        out,
        "K-mers enriched at particular positions. This is a heuristic: strong enrichment usually means adapters, primers or low-complexity sequence.",
    )?;
    if metrics.kmer_rows.is_empty() {
        writeln!(
            out,
            "<p class=\"desc\">No k-mer showed positional enrichment above 3x.</p>"
        )?;
        return Ok(());
    }
    table(out, "Data", |o| {
        writeln!(
            o,
            "<tr><th>Sequence</th><th>Count</th><th>P-value</th><th>Obs/Exp Max</th><th>Max Obs/Exp Position</th></tr>"
        )?;
        for r in &metrics.kmer_rows {
            writeln!(
                o,
                "<tr><td class=\"seq\">{}</td><td>{}</td><td>{:.3e}</td><td>{:.2}</td><td>{}</td></tr>",
                esc(&r.sequence),
                fmt_int(r.count),
                r.p_value,
                r.obs_exp,
                r.max_pos
            )?;
        }
        Ok(())
    })
}

// ---------------------------------------------------------------- plotting --

struct Frame {
    left: f64,
    top: f64,
    plot_w: f64,
    plot_h: f64,
}

fn open_plot(out: &mut String, w: f64, h: f64) -> Result<Frame> {
    writeln!(out, "<div class=\"plot\">")?;
    writeln!(
        out,
        "<svg width=\"{}\" height=\"{}\" viewBox=\"0 0 {} {}\">",
        w, h, w, h
    )?;
    let left = 50.0;
    let right = 20.0;
    let top = 12.0;
    let bottom = 34.0;
    let plot_w = w - left - right;
    let plot_h = h - top - bottom;
    writeln!(
        out,
        "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"#fff\" stroke=\"#ddd\"/>",
        left, top, plot_w, plot_h
    )?;
    Ok(Frame {
        left,
        top,
        plot_w,
        plot_h,
    })
}

fn close_plot(out: &mut String) -> Result<()> {
    writeln!(out, "</svg></div>")?;
    Ok(())
}

fn svg_boxplot(
    out: &mut String,
    rows: &[crate::core::metrics::PerBaseQualRow],
    w: f64,
    h: f64,
    max_q: f64,
) -> Result<()> {
    let f = open_plot(out, w, h)?;
    draw_y_bands(
        out,
        &f,
        0.0,
        max_q,
        &[
            (0.0, 20.0, "#f4c7c3"),
            (20.0, 28.0, "#ffe5b4"),
            (28.0, max_q.max(28.0), "#cdeccf"),
        ],
    )?;

    let n = rows.len().max(1) as f64;
    let x_step = f.plot_w / n;
    let y_scale = if max_q <= 0.0 { 1.0 } else { f.plot_h / max_q };
    draw_y_axis_ticks(out, &f, 0.0, max_q, 5, true)?;
    draw_category_axis(
        out,
        &f,
        &rows.iter().map(|r| r.group.label()).collect::<Vec<_>>(),
    )?;
    draw_axis_labels(out, &f, "Position", "Quality")?;

    for (i, r) in rows.iter().enumerate() {
        let x = f.left + (i as f64 + 0.5) * x_step;
        let y_m = f.top + f.plot_h - (r.median as f64 * y_scale);
        let y_lq = f.top + f.plot_h - (r.lower_quartile as f64 * y_scale);
        let y_uq = f.top + f.plot_h - (r.upper_quartile as f64 * y_scale);
        let y_p10 = f.top + f.plot_h - (r.p10 as f64 * y_scale);
        let y_p90 = f.top + f.plot_h - (r.p90 as f64 * y_scale);
        let box_w = (x_step * 0.8).max(1.0);
        let box_x = x - box_w / 2.0;
        let color = if r.median >= 28 {
            "#cdeccf"
        } else if r.median >= 20 {
            "#ffe5b4"
        } else {
            "#f4c7c3"
        };
        writeln!(
            out,
            "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"{}\" stroke=\"#666\"/>",
            box_x,
            y_uq,
            box_w,
            (y_lq - y_uq).max(0.0),
            color
        )?;
        writeln!(
            out,
            "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"#555\" stroke-width=\"1\"/>",
            x, y_p90, x, y_p10
        )?;
        let cap_w = (box_w * 0.6).max(1.0);
        for y in [y_p90, y_p10] {
            writeln!(
                out,
                "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"#555\" stroke-width=\"1\"/>",
                x - cap_w / 2.0,
                y,
                x + cap_w / 2.0,
                y
            )?;
        }
        writeln!(
            out,
            "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"#333\" stroke-width=\"1.2\"/>",
            box_x,
            y_m,
            box_x + box_w,
            y_m
        )?;
    }
    close_plot(out)
}

fn svg_tile_heatmap(out: &mut String, metrics: &FinalMetrics, w: f64, h: f64) -> Result<()> {
    let f = open_plot(out, w, h)?;
    if metrics.per_tile.is_empty() {
        close_plot(out)?;
        return Ok(());
    }

    let mut tiles: Vec<u32> = metrics.per_tile.iter().map(|r| r.tile).collect();
    tiles.sort_unstable();
    tiles.dedup();
    let mut groups: Vec<String> = Vec::new();
    for r in &metrics.per_tile {
        let label = r.group.label();
        if !groups.contains(&label) {
            groups.push(label);
        }
    }
    if tiles.is_empty() || groups.is_empty() {
        close_plot(out)?;
        return Ok(());
    }

    let cell_w = f.plot_w / groups.len() as f64;
    let cell_h = f.plot_h / tiles.len() as f64;
    let max_abs = metrics
        .per_tile
        .iter()
        .map(|r| r.mean_deviation.abs())
        .fold(1.0, f64::max);

    for r in &metrics.per_tile {
        let Some(row) = tiles.iter().position(|&t| t == r.tile) else {
            continue;
        };
        let label = r.group.label();
        let Some(col) = groups.iter().position(|g| *g == label) else {
            continue;
        };
        let x = f.left + col as f64 * cell_w;
        let y = f.top + row as f64 * cell_h;
        let t = (r.mean_deviation / max_abs).clamp(-1.0, 1.0);
        // Blue for above average, red for below.
        let color = if t >= 0.0 {
            let v = (255.0 - 120.0 * t) as u32;
            format!("#{:02x}{:02x}ff", v, v)
        } else {
            let v = (255.0 + 120.0 * t) as u32;
            format!("#ff{:02x}{:02x}", v, v)
        };
        writeln!(
            out,
            "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"{}\"/>",
            x,
            y,
            cell_w.max(0.5),
            cell_h.max(0.5),
            color
        )?;
    }

    draw_category_axis(out, &f, &groups)?;
    draw_axis_labels(out, &f, "Position", "Tile")?;
    close_plot(out)
}

fn svg_bars(
    out: &mut String,
    data: &[(f64, f64)],
    w: f64,
    h: f64,
    x_label: &str,
    y_label: &str,
) -> Result<()> {
    let f = open_plot(out, w, h)?;
    let max_y = data.iter().map(|(_, y)| *y).fold(0.0, f64::max);
    let x_min = data.first().map(|d| d.0).unwrap_or(0.0);
    let x_max = data.last().map(|d| d.0).unwrap_or(1.0);
    let bar_w = if data.is_empty() {
        1.0
    } else {
        (f.plot_w / data.len() as f64).max(1.0)
    };
    draw_y_axis_ticks(out, &f, 0.0, max_y, 4, false)?;
    draw_x_axis_ticks(out, &f, x_min, x_max, 5, false)?;
    draw_axis_labels(out, &f, x_label, y_label)?;
    for (i, (_xv, yv)) in data.iter().enumerate() {
        let x = f.left + (i as f64) * bar_w;
        let y = if max_y == 0.0 {
            0.0
        } else {
            yv / max_y * f.plot_h
        };
        writeln!(
            out,
            "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"#8ecae6\"/>",
            x,
            f.top + f.plot_h - y,
            bar_w,
            y
        )?;
    }
    close_plot(out)
}

fn svg_bars_xbands(
    out: &mut String,
    data: &[(f64, f64)],
    w: f64,
    h: f64,
    min_x: f64,
    max_x: f64,
    bands: &[(f64, f64, &str)],
    x_label: &str,
    y_label: &str,
) -> Result<()> {
    let f = open_plot(out, w, h)?;
    let x_range = (max_x - min_x).max(1.0);
    for (lo, hi, color) in bands {
        let start = ((*lo - min_x) / x_range).clamp(0.0, 1.0);
        let end = ((*hi - min_x) / x_range).clamp(0.0, 1.0);
        let x = f.left + start * f.plot_w;
        let w_band = (end - start) * f.plot_w;
        if w_band > 0.0 {
            writeln!(
                out,
                "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"{}\" opacity=\"0.18\"/>",
                x, f.top, w_band, f.plot_h, color
            )?;
        }
    }
    let max_y = data.iter().map(|(_, y)| *y).fold(0.0, f64::max);
    draw_y_axis_ticks(out, &f, 0.0, max_y, 4, true)?;
    draw_x_axis_ticks(out, &f, min_x, max_x, 5, true)?;
    draw_axis_labels(out, &f, x_label, y_label)?;
    for (xv, yv) in data.iter() {
        let x = f.left + ((*xv - min_x) / x_range).clamp(0.0, 1.0) * f.plot_w;
        let bar_w = (f.plot_w / x_range).max(1.0);
        let y = if max_y == 0.0 {
            0.0
        } else {
            yv / max_y * f.plot_h
        };
        writeln!(
            out,
            "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"#7db8da\"/>",
            x,
            f.top + f.plot_h - y,
            bar_w,
            y
        )?;
    }
    close_plot(out)
}

fn svg_category_bars(
    out: &mut String,
    labels: &[String],
    values: &[f64],
    w: f64,
    h: f64,
    x_label: &str,
    y_label: &str,
) -> Result<()> {
    let f = open_plot(out, w, h)?;
    let max_y = values.iter().copied().fold(0.0, f64::max);
    let bar_w = if values.is_empty() {
        1.0
    } else {
        f.plot_w / values.len() as f64
    };
    draw_y_axis_ticks(out, &f, 0.0, max_y, 4, true)?;
    draw_category_axis(out, &f, labels)?;
    draw_axis_labels(out, &f, x_label, y_label)?;
    for (i, v) in values.iter().enumerate() {
        let y = if max_y == 0.0 {
            0.0
        } else {
            v / max_y * f.plot_h
        };
        writeln!(
            out,
            "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"#8ecae6\"/>",
            f.left + i as f64 * bar_w + bar_w * 0.1,
            f.top + f.plot_h - y,
            (bar_w * 0.8).max(1.0),
            y
        )?;
    }
    close_plot(out)
}

fn svg_category_lines(
    out: &mut String,
    labels: &[String],
    series: &[(&str, &[f64], &str)],
    w: f64,
    h: f64,
    x_label: &str,
    y_label: &str,
) -> Result<()> {
    let f = open_plot(out, w, h)?;
    let max_y = series
        .iter()
        .flat_map(|(_, v, _)| v.iter().copied())
        .fold(0.0, f64::max)
        .max(1.0);
    draw_y_axis_ticks(out, &f, 0.0, max_y, 5, true)?;
    draw_category_axis(out, &f, labels)?;
    draw_axis_labels(out, &f, x_label, y_label)?;

    let n = labels.len().max(1);
    let step = f.plot_w / n as f64;
    for (name, values, color) in series {
        let mut path = String::new();
        for (i, v) in values.iter().enumerate() {
            let x = f.left + (i as f64 + 0.5) * step;
            let y = f.top + f.plot_h - (v / max_y * f.plot_h);
            if i == 0 {
                write!(path, "M {} {}", x, y)?;
            } else {
                write!(path, " L {} {}", x, y)?;
            }
        }
        writeln!(
            out,
            "<path d=\"{}\" fill=\"none\" stroke=\"{}\" stroke-width=\"1.5\"><title>{}</title></path>",
            path,
            color,
            esc(name)
        )?;
    }

    // Inline legend.
    let mut lx = f.left + 8.0;
    for (name, _, color) in series {
        writeln!(
            out,
            "<rect x=\"{}\" y=\"{}\" width=\"10\" height=\"3\" fill=\"{}\"/><text x=\"{}\" y=\"{}\" font-size=\"10\" fill=\"#444\">{}</text>",
            lx,
            f.top + 8.0,
            color,
            lx + 14.0,
            f.top + 12.0,
            esc(name)
        )?;
        lx += 14.0 + name.len() as f64 * 6.0 + 16.0;
    }
    close_plot(out)
}

fn svg_base_content(
    out: &mut String,
    rows: &[crate::core::metrics::PerBaseContentRow],
    w: f64,
    h: f64,
) -> Result<()> {
    let f = open_plot(out, w, h)?;
    let mut min_y: f64 = 100.0;
    let mut max_y: f64 = 0.0;
    for r in rows {
        min_y = min_y.min(r.a.min(r.c.min(r.g.min(r.t))));
        max_y = max_y.max(r.a.max(r.c.max(r.g.max(r.t))));
    }
    if rows.is_empty() {
        min_y = 0.0;
        max_y = 100.0;
    }
    let pad = ((max_y - min_y) * 0.2).max(2.0);
    let y_min = (min_y - pad).max(0.0);
    let y_max = (max_y + pad).min(100.0);
    draw_y_axis_ticks(out, &f, y_min, y_max, 5, true)?;
    draw_category_axis(
        out,
        &f,
        &rows.iter().map(|r| r.group.label()).collect::<Vec<_>>(),
    )?;
    draw_axis_labels(out, &f, "Position", "%")?;

    type ContentRow = crate::core::metrics::PerBaseContentRow;
    type ContentSeries = (&'static str, fn(&ContentRow) -> f64);
    let series: [ContentSeries; 4] = [
        ("#882255", |r: &ContentRow| r.g),
        ("#332288", |r: &ContentRow| r.a),
        ("#117733", |r: &ContentRow| r.t),
        ("#DDCC77", |r: &ContentRow| r.c),
    ];
    let n = rows.len().max(1);
    let step = f.plot_w / n as f64;
    for (color, get) in series.iter() {
        let mut path = String::new();
        for (i, r) in rows.iter().enumerate() {
            let x = f.left + (i as f64 + 0.5) * step;
            let y = f.top + f.plot_h - ((get(r) - y_min) / (y_max - y_min).max(1e-6) * f.plot_h);
            if i == 0 {
                write!(path, "M {} {}", x, y)?;
            } else {
                write!(path, " L {} {}", x, y)?;
            }
        }
        writeln!(
            out,
            "<path d=\"{}\" fill=\"none\" stroke=\"{}\" stroke-width=\"1.5\"/>",
            path, color
        )?;
    }
    close_plot(out)
}

fn svg_line_ybands(
    out: &mut String,
    data: &[(f64, f64)],
    w: f64,
    h: f64,
    min_y: f64,
    max_y: f64,
    color: &str,
    bands: &[(f64, f64, &str)],
    x_label: &str,
    y_label: &str,
) -> Result<()> {
    let f = open_plot(out, w, h)?;
    draw_y_bands(out, &f, min_y, max_y, bands)?;
    draw_y_axis_ticks(out, &f, min_y, max_y, 5, true)?;
    let x_min = data.first().map(|d| d.0).unwrap_or(0.0);
    let x_max = data.last().map(|d| d.0).unwrap_or(1.0);
    draw_x_axis_ticks(out, &f, x_min, x_max, 5, true)?;
    draw_axis_labels(out, &f, x_label, y_label)?;
    svg_polyline(out, data, &f, min_y, max_y, color)?;
    close_plot(out)
}

fn svg_adapter_lines(out: &mut String, metrics: &FinalMetrics, w: f64, h: f64) -> Result<()> {
    let f = open_plot(out, w, h)?;
    let n_adapters = metrics.adapter_names.len();
    let mut series: Vec<Vec<f64>> = vec![Vec::new(); n_adapters];
    for r in &metrics.adapter_content {
        for (i, s) in series.iter_mut().enumerate() {
            s.push(r.values.get(i).copied().unwrap_or(0.0));
        }
    }
    let max_y = series
        .iter()
        .flat_map(|s| s.iter().copied())
        .fold(0.0, f64::max)
        .max(1.0);
    draw_y_axis_ticks(out, &f, 0.0, max_y, 5, true)?;
    draw_category_axis(
        out,
        &f,
        &metrics
            .adapter_content
            .iter()
            .map(|r| r.group.label())
            .collect::<Vec<_>>(),
    )?;
    draw_axis_labels(out, &f, "Position", "%")?;

    let colors = [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
    ];
    let n = metrics.adapter_content.len().max(1);
    let step = f.plot_w / n as f64;
    let mut lx = f.left + 8.0;
    for (i, values) in series.iter().enumerate() {
        let color = colors[i % colors.len()];
        let mut path = String::new();
        for (j, v) in values.iter().enumerate() {
            let x = f.left + (j as f64 + 0.5) * step;
            let y = f.top + f.plot_h - (v / max_y * f.plot_h);
            if j == 0 {
                write!(path, "M {} {}", x, y)?;
            } else {
                write!(path, " L {} {}", x, y)?;
            }
        }
        writeln!(
            out,
            "<path d=\"{}\" fill=\"none\" stroke=\"{}\" stroke-width=\"1.5\"/>",
            path, color
        )?;
        let name = &metrics.adapter_names[i];
        writeln!(
            out,
            "<rect x=\"{}\" y=\"{}\" width=\"10\" height=\"3\" fill=\"{}\"/><text x=\"{}\" y=\"{}\" font-size=\"9\" fill=\"#444\">{}</text>",
            lx,
            f.top + 8.0,
            color,
            lx + 13.0,
            f.top + 11.5,
            esc(name)
        )?;
        lx += 13.0 + name.len() as f64 * 5.0 + 12.0;
    }
    close_plot(out)
}

fn svg_polyline(
    out: &mut String,
    data: &[(f64, f64)],
    f: &Frame,
    min_y: f64,
    max_y: f64,
    color: &str,
) -> Result<()> {
    if data.is_empty() {
        return Ok(());
    }
    let x_min = data.first().map(|d| d.0).unwrap_or(0.0);
    let x_max = data.last().map(|d| d.0).unwrap_or(1.0);
    let x_range = (x_max - x_min).max(1.0);
    let y_range = (max_y - min_y).max(1e-6);

    let mut path = String::new();
    for (i, (xv, yv)) in data.iter().enumerate() {
        let x = f.left + (*xv - x_min) / x_range * f.plot_w;
        let y = f.top + f.plot_h - ((*yv - min_y) / y_range * f.plot_h);
        if i == 0 {
            write!(path, "M {} {}", x, y)?;
        } else {
            write!(path, " L {} {}", x, y)?;
        }
    }
    writeln!(
        out,
        "<path d=\"{}\" fill=\"none\" stroke=\"{}\" stroke-width=\"1.5\"/>",
        path, color
    )?;
    Ok(())
}

fn draw_y_bands(
    out: &mut String,
    f: &Frame,
    min_y: f64,
    max_y: f64,
    bands: &[(f64, f64, &str)],
) -> Result<()> {
    let y_range = (max_y - min_y).max(1.0);
    for (lo, hi, color) in bands {
        let start = ((*lo - min_y) / y_range).clamp(0.0, 1.0);
        let end = ((*hi - min_y) / y_range).clamp(0.0, 1.0);
        let y1 = f.top + f.plot_h - end * f.plot_h;
        let y2 = f.top + f.plot_h - start * f.plot_h;
        let h = (y2 - y1).max(0.0);
        if h > 0.0 {
            writeln!(
                out,
                "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"{}\" opacity=\"0.18\"/>",
                f.left, y1, f.plot_w, h, color
            )?;
        }
    }
    Ok(())
}

fn draw_y_axis_ticks(
    out: &mut String,
    f: &Frame,
    min_y: f64,
    max_y: f64,
    ticks: usize,
    gridlines: bool,
) -> Result<()> {
    if ticks < 2 || (max_y - min_y).abs() < 1e-9 {
        return Ok(());
    }
    let (start, step, count) = nice_ticks(min_y, max_y, ticks);
    for i in 0..count {
        let v = start + step * i as f64;
        if v < min_y - 1e-9 || v > max_y + 1e-9 {
            continue;
        }
        let y = f.top + f.plot_h - ((v - min_y) / (max_y - min_y).max(1e-6)) * f.plot_h;
        if gridlines {
            writeln!(
                out,
                "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"#eee\"/>",
                f.left,
                y,
                f.left + f.plot_w,
                y
            )?;
        }
        writeln!(
            out,
            "<text x=\"{}\" y=\"{}\" font-size=\"10\" fill=\"#666\" text-anchor=\"end\" dominant-baseline=\"middle\">{}</text>",
            f.left - 4.0,
            y,
            fmt_tick(v)
        )?;
    }
    Ok(())
}

fn draw_x_axis_ticks(
    out: &mut String,
    f: &Frame,
    min_x: f64,
    max_x: f64,
    ticks: usize,
    gridlines: bool,
) -> Result<()> {
    if ticks < 2 || (max_x - min_x).abs() < 1e-9 {
        return Ok(());
    }
    let (start, step, count) = nice_ticks(min_x, max_x, ticks);
    for i in 0..count {
        let v = start + step * i as f64;
        if v < min_x - 1e-9 || v > max_x + 1e-9 {
            continue;
        }
        let x = f.left + ((v - min_x) / (max_x - min_x).max(1e-6)) * f.plot_w;
        if gridlines {
            writeln!(
                out,
                "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"#eee\"/>",
                x,
                f.top,
                x,
                f.top + f.plot_h
            )?;
        }
        writeln!(
            out,
            "<text x=\"{}\" y=\"{}\" font-size=\"10\" fill=\"#666\" text-anchor=\"middle\" dominant-baseline=\"hanging\">{}</text>",
            x,
            f.top + f.plot_h + 4.0,
            fmt_tick(v)
        )?;
    }
    Ok(())
}

/// Label a categorical x axis, thinning labels so they stay readable.
fn draw_category_axis(out: &mut String, f: &Frame, labels: &[String]) -> Result<()> {
    if labels.is_empty() {
        return Ok(());
    }
    let step = f.plot_w / labels.len() as f64;
    // Roughly 60 px per label keeps text from overlapping.
    let stride = ((60.0 / step).ceil() as usize).max(1);
    for (i, label) in labels.iter().enumerate() {
        if i % stride != 0 && i + 1 != labels.len() {
            continue;
        }
        let x = f.left + (i as f64 + 0.5) * step;
        writeln!(
            out,
            "<text x=\"{}\" y=\"{}\" font-size=\"10\" fill=\"#666\" text-anchor=\"middle\" dominant-baseline=\"hanging\">{}</text>",
            x,
            f.top + f.plot_h + 4.0,
            esc(label)
        )?;
    }
    Ok(())
}

fn draw_axis_labels(out: &mut String, f: &Frame, x_label: &str, y_label: &str) -> Result<()> {
    let x = f.left + f.plot_w / 2.0;
    let y = f.top + f.plot_h + 22.0;
    writeln!(
        out,
        "<text x=\"{}\" y=\"{}\" font-size=\"11\" fill=\"#444\" text-anchor=\"middle\">{}</text>",
        x,
        y,
        esc(x_label)
    )?;
    let yx = f.left - 34.0;
    let yy = f.top + f.plot_h / 2.0;
    writeln!(
        out,
        "<text x=\"{}\" y=\"{}\" font-size=\"11\" fill=\"#444\" text-anchor=\"middle\" transform=\"rotate(-90 {} {})\">{}</text>",
        yx,
        yy,
        yx,
        yy,
        esc(y_label)
    )?;
    Ok(())
}

fn legend_base_content(out: &mut String) -> Result<()> {
    writeln!(
        out,
        "<div class=\"desc\"><b>Legend:</b> <span style=\"display:inline-block;width:18px;height:4px;background:#882255;margin:0 6px 2px 6px;vertical-align:middle;\"></span><b>G</b> <span style=\"display:inline-block;width:18px;height:4px;background:#332288;margin:0 6px 2px 10px;vertical-align:middle;\"></span><b>A</b> <span style=\"display:inline-block;width:18px;height:4px;background:#117733;margin:0 6px 2px 10px;vertical-align:middle;\"></span><b>T</b> <span style=\"display:inline-block;width:18px;height:4px;background:#DDCC77;margin:0 6px 2px 10px;vertical-align:middle;\"></span><b>C</b></div>"
    )?;
    Ok(())
}

fn fmt_tick(v: f64) -> String {
    if (v - v.round()).abs() < 0.001 {
        format!("{}", v.round() as i64)
    } else if v.abs() < 10.0 {
        format!("{:.2}", v)
    } else {
        format!("{:.1}", v)
    }
}

fn nice_ticks(min: f64, max: f64, ticks: usize) -> (f64, f64, usize) {
    let range = (max - min).abs().max(1e-9);
    let rough = range / (ticks as f64 - 1.0);
    let mag = 10f64.powf(rough.abs().log10().floor());
    let norm = rough / mag;
    let step = if norm <= 1.0 {
        1.0
    } else if norm <= 2.0 {
        2.0
    } else if norm <= 5.0 {
        5.0
    } else {
        10.0
    } * mag;
    let start = (min / step).floor() * step;
    let end = (max / step).ceil() * step;
    let count = (((end - start) / step).round() as usize + 1).min(64);
    (start, step, count)
}

fn auto_range<I: Iterator<Item = f64>>(values: I, min_bound: f64, max_bound: f64) -> (f64, f64) {
    let mut min_v = f64::INFINITY;
    let mut max_v = f64::NEG_INFINITY;
    for v in values {
        if v < min_v {
            min_v = v;
        }
        if v > max_v {
            max_v = v;
        }
    }
    if !min_v.is_finite() || !max_v.is_finite() {
        return (min_bound, max_bound);
    }
    let span = (max_v - min_v).max(1e-6);
    let pad = (span * 0.2).max(1.0);
    let y_min = (min_v - pad).max(min_bound);
    let y_max = (max_v + pad).min(max_bound);
    if (y_max - y_min) < 1e-6 {
        (min_bound, max_bound)
    } else {
        (y_min, y_max)
    }
}

// ------------------------------------------------------- LaTeX figure hooks --

#[cfg(feature = "latex")]
fn extract_svg(s: &str) -> String {
    if let (Some(start), Some(end)) = (s.find("<svg"), s.rfind("</svg>")) {
        s[start..end + 6].to_string()
    } else {
        s.to_string()
    }
}

#[cfg(feature = "latex")]
fn simple_text_svg(title: &str, lines: &[String]) -> String {
    let mut out = String::new();
    out.push_str("<svg width=\"800\" height=\"260\" viewBox=\"0 0 800 260\">");
    out.push_str(
        "<rect x=\"0\" y=\"0\" width=\"800\" height=\"260\" fill=\"#fff\" stroke=\"#ddd\"/>",
    );
    let _ = write!(
        out,
        "<text x=\"16\" y=\"28\" font-size=\"14\" fill=\"#333\">{}</text>",
        esc(title)
    );
    let mut y = 54;
    for l in lines {
        let _ = write!(
            out,
            "<text x=\"16\" y=\"{}\" font-size=\"12\" fill=\"#333\">{}</text>",
            y,
            esc(l)
        );
        y += 18;
    }
    out.push_str("</svg>");
    out
}

#[cfg(feature = "latex")]
pub(crate) fn latex_svg_per_base_quality(metrics: &FinalMetrics) -> Result<String> {
    let max_q = metrics
        .per_base_qual
        .iter()
        .map(|r| r.p90 as f64)
        .fold(40.0, f64::max);
    let mut s = String::new();
    svg_boxplot(&mut s, &metrics.per_base_qual, PLOT_W, PLOT_H, max_q)?;
    Ok(extract_svg(&s))
}

#[cfg(feature = "latex")]
pub(crate) fn latex_svg_per_tile(metrics: &FinalMetrics) -> Result<String> {
    let mut s = String::new();
    svg_tile_heatmap(&mut s, metrics, PLOT_W, PLOT_H)?;
    Ok(extract_svg(&s))
}

#[cfg(feature = "latex")]
pub(crate) fn latex_svg_per_seq_quality(metrics: &FinalMetrics) -> Result<String> {
    let data: Vec<(f64, f64)> = metrics
        .per_seq_qual
        .iter()
        .map(|r| (r.mean_q as f64, r.count as f64))
        .collect();
    let mut s = String::new();
    svg_bars(&mut s, &data, PLOT_W, PLOT_H, "Mean Q", "Count")?;
    Ok(extract_svg(&s))
}

#[cfg(feature = "latex")]
pub(crate) fn latex_svg_duplication(metrics: &FinalMetrics) -> Result<String> {
    let labels: Vec<String> = metrics
        .duplication
        .iter()
        .map(|r| r.level.to_string())
        .collect();
    let dedup: Vec<f64> = metrics
        .duplication
        .iter()
        .map(|r| r.percent_of_deduplicated)
        .collect();
    let total: Vec<f64> = metrics
        .duplication
        .iter()
        .map(|r| r.percent_of_total)
        .collect();
    let mut s = String::new();
    svg_category_lines(
        &mut s,
        &labels,
        &[
            ("% Deduplicated", &dedup, "#1f77b4"),
            ("% Total", &total, "#d62728"),
        ],
        PLOT_W,
        PLOT_H,
        "Duplication level",
        "%",
    )?;
    Ok(extract_svg(&s))
}

#[cfg(feature = "latex")]
pub(crate) fn latex_svg_adapter_content(metrics: &FinalMetrics) -> Result<String> {
    let mut s = String::new();
    svg_adapter_lines(&mut s, metrics, PLOT_W, PLOT_H)?;
    Ok(extract_svg(&s))
}

#[cfg(feature = "latex")]
pub(crate) fn latex_svg_per_base_content(metrics: &FinalMetrics) -> Result<String> {
    let mut s = String::new();
    svg_base_content(&mut s, &metrics.per_base_content, PLOT_W, PLOT_H)?;
    Ok(extract_svg(&s))
}

#[cfg(feature = "latex")]
pub(crate) fn latex_svg_per_seq_gc(metrics: &FinalMetrics) -> Result<String> {
    let data: Vec<(f64, f64)> = metrics
        .per_seq_gc
        .iter()
        .map(|r| (r.gc as f64, r.count as f64))
        .collect();
    let mut s = String::new();
    svg_bars_xbands(
        &mut s,
        &data,
        PLOT_W,
        PLOT_H,
        0.0,
        100.0,
        &[(40.0, 60.0, "#cdeccf")],
        "GC%",
        "Count",
    )?;
    Ok(extract_svg(&s))
}

#[cfg(feature = "latex")]
pub(crate) fn latex_svg_per_base_n(metrics: &FinalMetrics) -> Result<String> {
    let data: Vec<(f64, f64)> = metrics
        .per_base_n
        .iter()
        .map(|r| (r.group.mid(), r.n_percent))
        .collect();
    let (y_min, y_max) = auto_range(data.iter().map(|(_, y)| *y), 0.0, 100.0);
    let mut s = String::new();
    svg_line_ybands(
        &mut s,
        &data,
        PLOT_W,
        PLOT_H,
        y_min,
        y_max,
        "#555",
        &[
            (0.0, 5.0, "#cdeccf"),
            (5.0, 20.0, "#ffe5b4"),
            (20.0, 100.0, "#f4c7c3"),
        ],
        "Position",
        "% N",
    )?;
    Ok(extract_svg(&s))
}

#[cfg(feature = "latex")]
pub(crate) fn latex_svg_per_seq_n(metrics: &FinalMetrics) -> Result<String> {
    let data: Vec<(f64, f64)> = metrics
        .per_seq_n
        .iter()
        .map(|r| (r.n_percent as f64, r.count as f64))
        .collect();
    let mut s = String::new();
    svg_bars(&mut s, &data, PLOT_W, PLOT_H, "N%", "Count")?;
    Ok(extract_svg(&s))
}

#[cfg(feature = "latex")]
pub(crate) fn latex_svg_length_dist(metrics: &FinalMetrics) -> Result<String> {
    let mut s = String::new();
    if let Some(ref ll) = metrics.long_length {
        let labels: Vec<String> = ll.bins.iter().map(|(l, _)| l.clone()).collect();
        let values: Vec<f64> = ll.bins.iter().map(|(_, c)| *c as f64).collect();
        svg_category_bars(
            &mut s,
            &labels,
            &values,
            PLOT_W,
            PLOT_H,
            "Length bin",
            "Count",
        )?;
    } else {
        let data: Vec<(f64, f64)> = metrics
            .length_dist
            .iter()
            .map(|r| (r.length as f64, r.count as f64))
            .collect();
        svg_bars(&mut s, &data, PLOT_W, PLOT_H, "Length", "Count")?;
    }
    Ok(extract_svg(&s))
}

#[cfg(feature = "latex")]
pub(crate) fn latex_svg_overrep(metrics: &FinalMetrics) -> Result<String> {
    let mut lines: Vec<String> = metrics
        .overrepresented
        .iter()
        .take(6)
        .map(|r| format!("{} ({:.2}%)", r.sequence, r.percent))
        .collect();
    if lines.is_empty() {
        lines.push("No overrepresented sequences".to_string());
    }
    Ok(simple_text_svg("Overrepresented sequences", &lines))
}

#[cfg(feature = "latex")]
pub(crate) fn latex_svg_kmer_content(metrics: &FinalMetrics) -> Result<String> {
    let mut lines: Vec<String> = metrics
        .kmer_rows
        .iter()
        .take(6)
        .map(|r| format!("{} (obs/exp {:.2})", r.sequence, r.obs_exp))
        .collect();
    if lines.is_empty() {
        lines.push("No enriched k-mers".to_string());
    }
    Ok(simple_text_svg("Kmer content", &lines))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn markup_is_escaped() {
        assert_eq!(esc("<script>"), "&lt;script&gt;");
        assert_eq!(esc("a&b"), "a&amp;b");
        assert_eq!(esc("\"q\""), "&quot;q&quot;");
    }

    #[test]
    fn nice_ticks_are_bounded() {
        let (_, _, count) = nice_ticks(0.0, 1e12, 5);
        assert!(count <= 64);
    }
}
