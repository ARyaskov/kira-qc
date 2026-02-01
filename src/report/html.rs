use crate::core::engine::RunOutput;
use crate::core::model::{Mode, Status};
use anyhow::{Context, Result};
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::{SystemTime, UNIX_EPOCH};

fn write_modern(path: &Path, output: &RunOutput) -> Result<()> {
    let metrics = output.agg.finalize(&output.ctx);
    let mut html = String::with_capacity(256 * 1024);
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);

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
        output.ctx.sample_name
    )?;
    writeln!(html, "<style>")?;
    writeln!(
        html,
        "body{{font-family:Arial,Helvetica,sans-serif;margin:20px;color:#222;background:#fff;}}"
    )?;
    writeln!(html, "h1{{margin:0 0 8px 0;font-size:24px;}}")?;
    writeln!(html, "h2{{margin:24px 0 8px 0;font-size:20px;}}")?;
    writeln!(
        html,
        ".meta{{color:#555;font-size:13px;margin-bottom:16px;}}"
    )?;
    writeln!(
        html,
        ".summary{{border-collapse:collapse;margin:12px 0 20px 0;width:100%;max-width:900px;}}"
    )?;
    writeln!(
        html,
        ".summary th,.summary td{{border:1px solid #ddd;padding:6px 10px;text-align:left;}}"
    )?;
    writeln!(html, ".pass{{color:#0a7a0a;font-weight:bold;}}")?;
    writeln!(html, ".warn{{color:#d98200;font-weight:bold;}}")?;
    writeln!(html, ".fail{{color:#c00000;font-weight:bold;}}")?;
    writeln!(
        html,
        ".module{{border-top:1px solid #eee;padding-top:8px;}}"
    )?;
    writeln!(html, ".plot{{margin:8px 0 6px 0;}}")?;
    writeln!(
        html,
        ".desc{{color:#444;font-size:13px;max-width:1000px;margin:4px 0 10px 0;}}"
    )?;
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
    writeln!(html, "details{{margin:6px 0 18px 0;}}")?;
    writeln!(html, "svg{{background:#fafafa;border:1px solid #e5e5e5;}}")?;
    writeln!(html, "</style>")?;
    writeln!(html, "</head>")?;
    writeln!(html, "<body>")?;

    let mode_label = match output.ctx.mode {
        Mode::Short => "Short-read (Illumina)",
        Mode::Long => "Long-read (ONT / PacBio)",
    };

    writeln!(html, "<h1>kira-qc report</h1>")?;
    writeln!(
        html,
        "<div class=\"meta\">Sample: <b>{}</b><br/>File: {}<br/>Mode: {}<br/>Timestamp: {} (unix: {})</div>",
        output.ctx.sample_name,
        output.ctx.file_name,
        mode_label,
        fmt_timestamp(ts),
        ts
    )?;

    writeln!(html, "<h2>Summary</h2>")?;
    writeln!(html, "<table class=\"summary\">")?;
    writeln!(html, "<tr><th>Status</th><th>Module</th></tr>")?;
    summary_row(&mut html, metrics.statuses.basic, "Basic Statistics")?;
    match output.ctx.mode {
        Mode::Short => {
            summary_row(
                &mut html,
                metrics.statuses.per_base_qual,
                "Per base sequence quality",
            )?;
            summary_row(
                &mut html,
                metrics.statuses.per_seq_qual,
                "Per sequence quality scores",
            )?;
            summary_row(
                &mut html,
                metrics.statuses.per_base_content,
                "Per base sequence content",
            )?;
            summary_row(
                &mut html,
                metrics.statuses.per_seq_gc,
                "Per sequence GC content",
            )?;
            summary_row(&mut html, metrics.statuses.per_base_n, "Per base N content")?;
            summary_row(
                &mut html,
                metrics.statuses.length_dist,
                "Sequence Length Distribution",
            )?;
            summary_row(
                &mut html,
                metrics.statuses.duplication,
                "Sequence Duplication Levels",
            )?;
            summary_row(
                &mut html,
                metrics.statuses.overrepresented,
                "Overrepresented sequences",
            )?;
            summary_row(
                &mut html,
                metrics.statuses.adapter_content,
                "Adapter Content",
            )?;
            #[cfg(not(feature = "no-kmer"))]
            summary_row(&mut html, metrics.statuses.kmer_content, "Kmer Content")?;
        }
        Mode::Long => {
            summary_row(
                &mut html,
                metrics.statuses.length_dist,
                "Sequence Length Distribution",
            )?;
            summary_row(
                &mut html,
                metrics.statuses.per_seq_qual,
                "Per sequence quality scores",
            )?;
            summary_row(
                &mut html,
                metrics.statuses.per_seq_gc,
                "Per sequence GC content",
            )?;
            summary_row(
                &mut html,
                metrics.statuses.per_seq_n,
                "Per sequence N content",
            )?;
            summary_row(
                &mut html,
                metrics.statuses.adapter_content,
                "Adapter Content",
            )?;
        }
    }
    writeln!(html, "</table>")?;

    module_basic_stats(&mut html, &metrics, &output.ctx.file_name)?;
    match output.ctx.mode {
        Mode::Short => {
            module_per_base_quality(&mut html, &metrics)?;
            module_per_seq_quality(&mut html, &metrics)?;
            module_per_base_content(&mut html, &metrics)?;
            module_per_seq_gc(&mut html, &metrics)?;
            module_per_base_n(&mut html, &metrics)?;
            module_length_dist_short(&mut html, &metrics)?;
            module_duplication(&mut html, &metrics)?;
            module_overrep(&mut html, &metrics)?;
            module_adapter_content_short(&mut html, &metrics)?;
            #[cfg(not(feature = "no-kmer"))]
            module_kmer_content(&mut html, &metrics)?;
        }
        Mode::Long => {
            module_length_dist_long(&mut html, &metrics)?;
            module_per_seq_quality(&mut html, &metrics)?;
            module_per_seq_gc(&mut html, &metrics)?;
            module_per_seq_n(&mut html, &metrics)?;
            module_adapter_content_long(&mut html, &metrics)?;
        }
    }

    html.push_str("<script>");
    html.push_str(r#"document.querySelectorAll('table.sortable').forEach(t=>{const h=t.querySelectorAll('th');h.forEach((th,i)=>{th.style.cursor='pointer';th.addEventListener('click',()=>{const rows=[...t.querySelectorAll('tr')].slice(1);const asc=th.getAttribute('data-asc')!=='true';rows.sort((a,b)=>{const av=a.children[i].innerText;const bv=b.children[i].innerText;const an=parseFloat(av);const bn=parseFloat(bv);if(!isNaN(an)&&!isNaN(bn)){return asc?an-bn:bn-an;}return asc?av.localeCompare(bv):bv.localeCompare(av);});th.setAttribute('data-asc',asc);rows.forEach(r=>t.appendChild(r));});});});"#);
    html.push_str("</script>");
    writeln!(html, "</body></html>")?;

    let mut w =
        BufWriter::new(File::create(path).with_context(|| "create fastqc_report.html failed")?);
    w.write_all(html.as_bytes())?;
    Ok(())
}

pub fn write(path: &Path, output: &RunOutput) -> Result<()> {
    let metrics = output.agg.finalize(&output.ctx);
    let mut html = String::with_capacity(256 * 1024);
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);

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
        "<title>FastQC-compatible report: {}</title>",
        output.ctx.sample_name
    )?;
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
        ".main{{flex:1;background:#fff;border:1px solid #ddd;border-radius:4px;box-shadow:0 1px 3px rgba(0,0,0,0.08);padding:16px 20px;}}"
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
    writeln!(html, ".plot{{margin:8px 0 6px 0;}}")?;
    writeln!(
        html,
        ".desc{{color:#444;font-size:13px;max-width:1000px;margin:4px 0 10px 0;}}"
    )?;
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
    writeln!(html, "</head>")?;
    writeln!(html, "<body>")?;

    let mode_label = match output.ctx.mode {
        Mode::Short => "Short-read (Illumina)",
        Mode::Long => "Long-read (ONT / PacBio)",
    };

    writeln!(html, "<div class=\"page\">")?;
    writeln!(html, "<aside class=\"sidebar\">")?;
    writeln!(html, "<h2 id=\"summary\">Summary</h2>")?;
    writeln!(html, "<ul>")?;
    sidebar_item(
        &mut html,
        metrics.statuses.basic,
        "Basic Statistics",
        module_id_basic(),
    )?;
    match output.ctx.mode {
        Mode::Short => {
            sidebar_item(
                &mut html,
                metrics.statuses.per_base_qual,
                "Per base sequence quality",
                module_id_per_base_qual(),
            )?;
            sidebar_item(
                &mut html,
                metrics.statuses.per_seq_qual,
                "Per sequence quality scores",
                module_id_per_seq_qual(),
            )?;
            sidebar_item(
                &mut html,
                metrics.statuses.per_base_content,
                "Per base sequence content",
                module_id_per_base_content(),
            )?;
            sidebar_item(
                &mut html,
                metrics.statuses.per_seq_gc,
                "Per sequence GC content",
                module_id_per_seq_gc(),
            )?;
            sidebar_item(
                &mut html,
                metrics.statuses.per_base_n,
                "Per base N content",
                module_id_per_base_n(),
            )?;
            sidebar_item(
                &mut html,
                metrics.statuses.length_dist,
                "Sequence Length Distribution",
                module_id_length_dist(),
            )?;
            sidebar_item(
                &mut html,
                metrics.statuses.duplication,
                "Sequence Duplication Levels",
                module_id_duplication(),
            )?;
            sidebar_item(
                &mut html,
                metrics.statuses.overrepresented,
                "Overrepresented sequences",
                module_id_overrep(),
            )?;
            sidebar_item(
                &mut html,
                metrics.statuses.adapter_content,
                "Adapter Content",
                module_id_adapter_content(),
            )?;
            #[cfg(not(feature = "no-kmer"))]
            sidebar_item(
                &mut html,
                metrics.statuses.kmer_content,
                "Kmer Content",
                module_id_kmer(),
            )?;
        }
        Mode::Long => {
            sidebar_item(
                &mut html,
                metrics.statuses.length_dist,
                "Sequence Length Distribution",
                module_id_length_dist(),
            )?;
            sidebar_item(
                &mut html,
                metrics.statuses.per_seq_qual,
                "Per sequence quality scores",
                module_id_per_seq_qual(),
            )?;
            sidebar_item(
                &mut html,
                metrics.statuses.per_seq_gc,
                "Per sequence GC content",
                module_id_per_seq_gc(),
            )?;
            sidebar_item(
                &mut html,
                metrics.statuses.per_seq_n,
                "Per sequence N content",
                module_id_per_seq_n(),
            )?;
            sidebar_item(
                &mut html,
                metrics.statuses.adapter_content,
                "Adapter Content",
                module_id_adapter_content(),
            )?;
        }
    }
    writeln!(html, "</ul>")?;
    writeln!(html, "</aside>")?;

    writeln!(html, "<main class=\"main\">")?;
    writeln!(html, "<h1>kira-qc FastQC-compatible Report</h1>")?;
    writeln!(
        html,
        "<div class=\"meta\">File: {}<br/>Mode: {}<br/>Timestamp: {} (unix: {})</div>",
        output.ctx.file_name,
        mode_label,
        fmt_timestamp(ts),
        ts
    )?;

    compat_basic_stats(&mut html, &metrics, &output.ctx.file_name)?;
    match output.ctx.mode {
        Mode::Short => {
            compat_per_base_quality(&mut html, &metrics)?;
            compat_per_seq_quality(&mut html, &metrics)?;
            compat_per_base_content(&mut html, &metrics)?;
            compat_per_seq_gc(&mut html, &metrics)?;
            compat_per_base_n(&mut html, &metrics)?;
            compat_length_dist_short(&mut html, &metrics)?;
            compat_duplication(&mut html, &metrics)?;
            compat_overrep(&mut html, &metrics)?;
            compat_adapter_content_short(&mut html, &metrics)?;
            #[cfg(not(feature = "no-kmer"))]
            compat_kmer_content(&mut html, &metrics)?;
        }
        Mode::Long => {
            compat_length_dist_long(&mut html, &metrics)?;
            compat_per_seq_quality(&mut html, &metrics)?;
            compat_per_seq_gc(&mut html, &metrics)?;
            compat_per_seq_n(&mut html, &metrics)?;
            compat_adapter_content_long(&mut html, &metrics)?;
        }
    }

    writeln!(html, "<div class=\"meta\">Produced by kira-qc</div>")?;
    writeln!(html, "</main>")?;
    writeln!(html, "</div>")?;
    writeln!(html, "</body></html>")?;

    let mut w =
        BufWriter::new(File::create(path).with_context(|| "create fastqc_compat.html failed")?);
    w.write_all(html.as_bytes())?;
    Ok(())
}

fn summary_row(out: &mut String, status: Status, name: &str) -> Result<()> {
    let class = status_class(status);
    writeln!(
        out,
        "<tr><td class=\"{}\">{}</td><td>{}</td></tr>",
        class,
        status.as_str_upper(),
        name
    )?;
    Ok(())
}

fn module_header(out: &mut String, status: Status, title: &str) -> Result<()> {
    let class = status_class(status);
    writeln!(out, "<div class=\"module\">")?;
    writeln!(out, "<h2 class=\"{}\">{}</h2>", class, title)?;
    Ok(())
}

fn module_desc(out: &mut String, text: &str) -> Result<()> {
    writeln!(out, "<p class=\"desc\">{}</p>", text)?;
    Ok(())
}

fn module_footer(out: &mut String) -> Result<()> {
    writeln!(out, "</div>")?;
    Ok(())
}

fn status_class(status: Status) -> &'static str {
    match status {
        Status::Pass => "pass",
        Status::Warn => "warn",
        Status::Fail => "fail",
    }
}

fn status_icon_svg(status: Status, size: u32) -> String {
    let (fill, mark) = match status {
        Status::Pass => ("#2e8b57", "M6 10 L10 14 L18 6"),
        Status::Warn => ("#e6a400", "M11 5 L11 13 M11 16 L11 18"),
        Status::Fail => ("#c00000", "M6 6 L18 18 M18 6 L6 18"),
    };
    format!(
        "<svg width=\"{s}\" height=\"{s}\" viewBox=\"0 0 24 24\" aria-hidden=\"true\"><circle cx=\"12\" cy=\"12\" r=\"11\" fill=\"{f}\"/><path d=\"{p}\" stroke=\"#fff\" stroke-width=\"2\" fill=\"none\" stroke-linecap=\"round\" stroke-linejoin=\"round\"/></svg>",
        s = size,
        f = fill,
        p = mark
    )
}

fn sidebar_item(out: &mut String, status: Status, name: &str, id: &str) -> Result<()> {
    writeln!(
        out,
        "<li>{} <a href=\"#{}\">{}</a></li>",
        status_icon_svg(status, 14),
        id,
        name
    )?;
    Ok(())
}

fn compat_section_header(out: &mut String, status: Status, title: &str, id: &str) -> Result<()> {
    writeln!(out, "<section id=\"{}\" class=\"module\">", id)?;
    writeln!(out, "<h2>{} {}</h2>", status_icon_svg(status, 16), title)?;
    Ok(())
}

fn compat_section_footer(out: &mut String) -> Result<()> {
    writeln!(
        out,
        "<a class=\"back\" href=\"#summary\">Back to Summary</a>"
    )?;
    writeln!(out, "</section>")?;
    Ok(())
}

fn module_id_basic() -> &'static str {
    "basic_statistics"
}
fn module_id_per_base_qual() -> &'static str {
    "per_base_sequence_quality"
}
fn module_id_per_seq_qual() -> &'static str {
    "per_sequence_quality_scores"
}
fn module_id_per_base_content() -> &'static str {
    "per_base_sequence_content"
}
fn module_id_per_seq_gc() -> &'static str {
    "per_sequence_gc_content"
}
fn module_id_per_base_n() -> &'static str {
    "per_base_n_content"
}
fn module_id_per_seq_n() -> &'static str {
    "per_sequence_n_content"
}
fn module_id_length_dist() -> &'static str {
    "sequence_length_distribution"
}
fn module_id_duplication() -> &'static str {
    "sequence_duplication_levels"
}
fn module_id_overrep() -> &'static str {
    "overrepresented_sequences"
}
fn module_id_adapter_content() -> &'static str {
    "adapter_content"
}
fn module_id_kmer() -> &'static str {
    "kmer_content"
}

fn table_with_summary<F>(out: &mut String, summary: &str, f: F) -> Result<()>
where
    F: FnOnce(&mut String) -> Result<()>,
{
    let mut tmp = String::new();
    f(&mut tmp)?;
    let open = "<details><summary>Table</summary>";
    let open_alt = "<details open><summary>Table</summary>";
    let repl = format!("<details><summary>{}</summary>", summary);
    let repl_open = format!("<details open><summary>{}</summary>", summary);
    let tmp = tmp.replace(open, &repl).replace(open_alt, &repl_open);
    out.push_str(&tmp);
    Ok(())
}

fn compat_basic_stats(
    out: &mut String,
    metrics: &crate::core::metrics::FinalMetrics,
    file: &str,
) -> Result<()> {
    compat_section_header(
        out,
        metrics.statuses.basic,
        "Basic Statistics",
        module_id_basic(),
    )?;
    module_desc(
        out,
        "Summary of input size, read counts, length range, and GC%. Large length ranges or unusual GC can indicate mixed libraries or contamination.",
    )?;
    writeln!(out, "<table class=\"bs-table\">")?;
    writeln!(out, "<tr><th>Measure</th><th>Value</th></tr>")?;
    writeln!(out, "<tr><td>Filename</td><td>{}</td></tr>", file)?;
    writeln!(
        out,
        "<tr><td>File type</td><td>{}</td></tr>",
        metrics.basic.file_type
    )?;
    writeln!(
        out,
        "<tr><td>Encoding</td><td>{}</td></tr>",
        metrics.basic.encoding
    )?;
    writeln!(
        out,
        "<tr><td>Total Sequences</td><td>{}</td></tr>",
        fmt_int(metrics.basic.total_sequences)
    )?;
    writeln!(
        out,
        "<tr><td>Filtered Sequences</td><td>{}</td></tr>",
        fmt_int(metrics.basic.filtered_sequences)
    )?;
    if metrics.basic.min_len == metrics.basic.max_len {
        writeln!(
            out,
            "<tr><td>Sequence length</td><td>{}</td></tr>",
            metrics.basic.min_len
        )?;
    } else {
        writeln!(
            out,
            "<tr><td>Sequence length</td><td>{}-{}</td></tr>",
            metrics.basic.min_len, metrics.basic.max_len
        )?;
    }
    writeln!(
        out,
        "<tr><td>%GC</td><td>{}</td></tr>",
        metrics.basic.gc_percent
    )?;
    writeln!(out, "</table>")?;
    compat_section_footer(out)
}

fn compat_per_base_quality(
    out: &mut String,
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<()> {
    compat_section_header(
        out,
        metrics.statuses.per_base_qual,
        "Per base sequence quality",
        module_id_per_base_qual(),
    )?;
    module_desc(
        out,
        "Shows quality score distributions at each base position. Systematic drops toward read ends often reflect sequencing degradation or adapter read-through.",
    )?;
    let (w, h) = (800.0, 260.0);
    let max_q = metrics
        .per_base_qual
        .iter()
        .map(|r| r.p90 as f64)
        .fold(40.0, f64::max);
    svg_boxplot(
        out,
        &metrics.per_base_qual,
        w,
        h,
        max_q,
        "Position",
        "Quality",
    )?;
    table_with_summary(out, "Data", |o| {
        table_per_base_quality(o, &metrics.per_base_qual)
    })?;
    compat_section_footer(out)
}

fn compat_per_seq_quality(
    out: &mut String,
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<()> {
    compat_section_header(
        out,
        metrics.statuses.per_seq_qual,
        "Per sequence quality scores",
        module_id_per_seq_qual(),
    )?;
    module_desc(
        out,
        "Shows the distribution of mean quality per read. A left-shifted distribution indicates overall low-quality reads or mixed data.",
    )?;
    let (w, h) = (800.0, 260.0);
    let data = metrics
        .per_seq_qual
        .iter()
        .map(|r| (r.mean_q as f64, r.count as f64))
        .collect::<Vec<_>>();
    svg_histogram_compat_bars(out, data.as_slice(), w, h, 0.0, 0.0, "Mean Q", "Count")?;
    table_with_summary(out, "Data", |o| {
        table_per_seq_quality(o, &metrics.per_seq_qual)
    })?;
    compat_section_footer(out)
}

fn compat_per_base_content(
    out: &mut String,
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<()> {
    compat_section_header(
        out,
        metrics.statuses.per_base_content,
        "Per base sequence content",
        module_id_per_base_content(),
    )?;
    module_desc(
        out,
        "Shows the percentage of each base at each position. Strong positional biases can indicate priming artifacts or residual adapters.",
    )?;
    let (w, h) = (800.0, 260.0);
    legend_base_content(out)?;
    svg_multi_line(out, &metrics.per_base_content, w, h, "Position", "%")?;
    table_with_summary(out, "Data", |o| {
        table_per_base_content(o, &metrics.per_base_content)
    })?;
    compat_section_footer(out)
}

fn compat_per_seq_gc(out: &mut String, metrics: &crate::core::metrics::FinalMetrics) -> Result<()> {
    compat_section_header(
        out,
        metrics.statuses.per_seq_gc,
        "Per sequence GC content",
        module_id_per_seq_gc(),
    )?;
    module_desc(
        out,
        "Shows the distribution of GC% across reads. Broad or multi-modal shapes can indicate contamination or mixed libraries.",
    )?;
    let (w, h) = (800.0, 260.0);
    let data = metrics
        .per_seq_gc
        .iter()
        .map(|r| (r.gc as f64, r.count as f64))
        .collect::<Vec<_>>();
    svg_histogram_xbands(
        out,
        data.as_slice(),
        w,
        h,
        0.0,
        100.0,
        &[(40.0, 60.0, "#cdeccf")],
        "GC%",
        "Count",
    )?;
    table_with_summary(out, "Data", |o| table_per_seq_gc(o, &metrics.per_seq_gc))?;
    compat_section_footer(out)
}

fn compat_per_base_n(out: &mut String, metrics: &crate::core::metrics::FinalMetrics) -> Result<()> {
    compat_section_header(
        out,
        metrics.statuses.per_base_n,
        "Per base N content",
        module_id_per_base_n(),
    )?;
    module_desc(
        out,
        "Shows the proportion of Ns at each position. Spikes or elevated Ns suggest base-calling issues or low-complexity regions.",
    )?;
    let (w, h) = (800.0, 260.0);
    let data = metrics
        .per_base_n
        .iter()
        .map(|r| (r.base as f64, r.n_percent))
        .collect::<Vec<_>>();
    let (y_min, y_max) = auto_range(data.iter().map(|(_, y)| *y), 0.0, 100.0);
    svg_single_line_ybands(
        out,
        data.as_slice(),
        w,
        h,
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
    table_with_summary(out, "Data", |o| table_per_base_n(o, &metrics.per_base_n))?;
    compat_section_footer(out)
}

fn compat_per_seq_n(out: &mut String, metrics: &crate::core::metrics::FinalMetrics) -> Result<()> {
    compat_section_header(
        out,
        metrics.statuses.per_seq_n,
        "Per sequence N content",
        module_id_per_seq_n(),
    )?;
    module_desc(
        out,
        "Shows the distribution of N% per read. Excess high-N reads indicate poor base-calling or low-quality segments.",
    )?;
    let data = metrics
        .per_seq_n
        .iter()
        .map(|r| (r.n_percent as f64, r.count as f64))
        .collect::<Vec<_>>();
    svg_histogram_xbands(
        out,
        data.as_slice(),
        800.0,
        260.0,
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
    table_with_summary(out, "Data", |o| table_per_seq_n(o, &metrics.per_seq_n))?;
    compat_section_footer(out)
}

fn compat_length_dist_short(
    out: &mut String,
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<()> {
    compat_section_header(
        out,
        metrics.statuses.length_dist,
        "Sequence Length Distribution",
        module_id_length_dist(),
    )?;
    module_desc(
        out,
        "Shows read length frequencies. Multiple peaks or long tails may indicate trimming or mixed read sources.",
    )?;
    let (w, h) = (800.0, 260.0);
    let data = metrics
        .length_dist
        .iter()
        .map(|r| (r.length as f64, r.count as f64))
        .collect::<Vec<_>>();
    svg_histogram_compat_bars(out, data.as_slice(), w, h, 0.0, 0.0, "Length", "Count")?;
    table_with_summary(out, "Data", |o| table_length_dist(o, &metrics.length_dist))?;
    compat_section_footer(out)
}

fn compat_length_dist_long(
    out: &mut String,
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<()> {
    compat_section_header(
        out,
        metrics.statuses.length_dist,
        "Sequence Length Distribution",
        module_id_length_dist(),
    )?;
    module_desc(
        out,
        "Shows read length frequencies using log-scaled bins. Very long tails or multiple modes may indicate mixed input or variable trimming.",
    )?;
    if let Some(ref ll) = metrics.long_length {
        let data = ll
            .bins
            .iter()
            .enumerate()
            .map(|(i, &c)| (i as f64 + 1.0, c as f64))
            .collect::<Vec<_>>();
        svg_histogram_compat_bars(
            out,
            data.as_slice(),
            800.0,
            260.0,
            0.0,
            0.0,
            "Length bin",
            "Count",
        )?;
        table_with_summary(out, "Data", |o| table_long_length(o, ll))?;
    }
    compat_section_footer(out)
}

fn compat_duplication(
    out: &mut String,
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<()> {
    compat_section_header(
        out,
        metrics.statuses.duplication,
        "Sequence Duplication Levels",
        module_id_duplication(),
    )?;
    module_desc(
        out,
        "Estimates duplication using a streaming heavy-hitter model. High duplication often indicates PCR over-amplification or low library complexity.",
    )?;
    let data = metrics
        .duplication
        .iter()
        .enumerate()
        .map(|(i, r)| (i as f64 + 1.0, r.relative))
        .collect::<Vec<_>>();
    svg_histogram_compat_bars(
        out,
        data.as_slice(),
        800.0,
        260.0,
        0.0,
        0.0,
        "Level",
        "Relative count",
    )?;
    table_with_summary(out, "Data", |o| table_duplication(o, &metrics.duplication))?;
    compat_section_footer(out)
}

fn compat_overrep(out: &mut String, metrics: &crate::core::metrics::FinalMetrics) -> Result<()> {
    compat_section_header(
        out,
        metrics.statuses.overrepresented,
        "Overrepresented sequences",
        module_id_overrep(),
    )?;
    module_desc(
        out,
        "Lists sequences occurring more often than expected. Common sources are adapters, primers, or contamination.",
    )?;
    table_with_summary(out, "Data", |o| table_overrep(o, &metrics.overrepresented))?;
    compat_section_footer(out)
}

fn compat_adapter_content_short(
    out: &mut String,
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<()> {
    compat_section_header(
        out,
        metrics.statuses.adapter_content,
        "Adapter Content",
        module_id_adapter_content(),
    )?;
    module_desc(
        out,
        "Shows adapter match percentages by position. Increasing signal toward read ends suggests adapter read-through.",
    )?;
    let (w, h) = (800.0, 260.0);
    svg_adapter_lines(out, &metrics.adapter_content, w, h, "Position", "%")?;
    table_with_summary(out, "Data", |o| {
        table_adapter_content(o, &metrics.adapter_content)
    })?;
    compat_section_footer(out)
}

fn compat_adapter_content_long(
    out: &mut String,
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<()> {
    compat_section_header(
        out,
        metrics.statuses.adapter_content,
        "Adapter Content",
        module_id_adapter_content(),
    )?;
    module_desc(
        out,
        "Reports the fraction of reads containing common adapter motifs. Elevated percentages suggest residual adapters or chimeric reads.",
    )?;
    table_with_summary(out, "Data", |o| {
        table_adapter_summary(o, &metrics.adapter_content)
    })?;
    compat_section_footer(out)
}

#[cfg(not(feature = "no-kmer"))]
fn compat_kmer_content(
    out: &mut String,
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<()> {
    compat_section_header(
        out,
        metrics.statuses.kmer_content,
        "Kmer Content",
        module_id_kmer(),
    )?;
    module_desc(
        out,
        "Reports k-mers enriched at specific positions. Strong enrichment can indicate adapters or sequence bias.",
    )?;
    table_with_summary(out, "Data", |o| table_kmer(o, &metrics.kmer_rows))?;
    compat_section_footer(out)
}

fn module_basic_stats(
    out: &mut String,
    metrics: &crate::core::metrics::FinalMetrics,
    file: &str,
) -> Result<()> {
    module_header(out, metrics.statuses.basic, "Basic Statistics")?;
    module_desc(
        out,
        "Summary of input size, read counts, length range, and GC%. Large length ranges or unusual GC can indicate mixed libraries or contamination.",
    )?;
    writeln!(out, "<table class=\"table\">")?;
    writeln!(out, "<tr><th>Measure</th><th>Value</th></tr>")?;
    writeln!(out, "<tr><td>Filename</td><td>{}</td></tr>", file)?;
    writeln!(
        out,
        "<tr><td>File type</td><td>{}</td></tr>",
        metrics.basic.file_type
    )?;
    writeln!(
        out,
        "<tr><td>Encoding</td><td>{}</td></tr>",
        metrics.basic.encoding
    )?;
    writeln!(
        out,
        "<tr><td>Total Sequences</td><td>{}</td></tr>",
        fmt_int(metrics.basic.total_sequences)
    )?;
    writeln!(
        out,
        "<tr><td>Filtered Sequences</td><td>{}</td></tr>",
        fmt_int(metrics.basic.filtered_sequences)
    )?;
    if metrics.basic.min_len == metrics.basic.max_len {
        writeln!(
            out,
            "<tr><td>Sequence length</td><td>{}</td></tr>",
            metrics.basic.min_len
        )?;
    } else {
        writeln!(
            out,
            "<tr><td>Sequence length</td><td>{}-{}</td></tr>",
            metrics.basic.min_len, metrics.basic.max_len
        )?;
    }
    writeln!(
        out,
        "<tr><td>%GC</td><td>{}</td></tr>",
        metrics.basic.gc_percent
    )?;
    writeln!(out, "</table>")?;
    module_footer(out)
}

fn module_per_base_quality(
    out: &mut String,
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<()> {
    module_header(
        out,
        metrics.statuses.per_base_qual,
        "Per base sequence quality",
    )?;
    module_desc(
        out,
        "Shows quality score distributions at each base position. Systematic drops toward read ends often reflect sequencing degradation or adapter read-through.",
    )?;
    let (w, h) = (800.0, 260.0);
    let max_q = metrics
        .per_base_qual
        .iter()
        .map(|r| r.p90 as f64)
        .fold(40.0, f64::max);
    svg_boxplot(
        out,
        &metrics.per_base_qual,
        w,
        h,
        max_q,
        "Position",
        "Quality",
    )?;
    table_per_base_quality(out, &metrics.per_base_qual)?;
    module_footer(out)
}

fn module_per_seq_quality(
    out: &mut String,
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<()> {
    module_header(
        out,
        metrics.statuses.per_seq_qual,
        "Per sequence quality scores",
    )?;
    module_desc(
        out,
        "Shows the distribution of mean quality per read. A left-shifted distribution indicates overall low-quality reads or mixed data.",
    )?;
    let (w, h) = (800.0, 260.0);
    let data = metrics
        .per_seq_qual
        .iter()
        .map(|r| (r.mean_q as f64, r.count as f64))
        .collect::<Vec<_>>();
    svg_histogram_xbands(
        out,
        data.as_slice(),
        w,
        h,
        0.0,
        0.0,
        &[
            (0.0, 20.0, "#f4c7c3"),
            (20.0, 28.0, "#ffe5b4"),
            (28.0, 60.0, "#cdeccf"),
        ],
        "Mean Q",
        "Count",
    )?;
    table_per_seq_quality(out, &metrics.per_seq_qual)?;
    module_footer(out)
}

fn module_per_base_content(
    out: &mut String,
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<()> {
    module_header(
        out,
        metrics.statuses.per_base_content,
        "Per base sequence content",
    )?;
    module_desc(
        out,
        "Shows the percentage of each base at each position. Strong positional biases can indicate priming artifacts or residual adapters.",
    )?;
    let (w, h) = (800.0, 260.0);
    legend_base_content(out)?;
    svg_multi_line(out, &metrics.per_base_content, w, h, "Position", "%")?;
    table_per_base_content(out, &metrics.per_base_content)?;
    module_footer(out)
}

fn module_per_seq_gc(out: &mut String, metrics: &crate::core::metrics::FinalMetrics) -> Result<()> {
    module_header(out, metrics.statuses.per_seq_gc, "Per sequence GC content")?;
    module_desc(
        out,
        "Shows the distribution of GC% across reads. Broad or multi-modal shapes can indicate contamination or mixed libraries.",
    )?;
    let (w, h) = (800.0, 260.0);
    let data = metrics
        .per_seq_gc
        .iter()
        .map(|r| (r.gc as f64, r.count as f64))
        .collect::<Vec<_>>();
    svg_histogram_xbands(
        out,
        data.as_slice(),
        w,
        h,
        0.0,
        100.0,
        &[(40.0, 60.0, "#cdeccf")],
        "GC%",
        "Count",
    )?;
    table_per_seq_gc(out, &metrics.per_seq_gc)?;
    module_footer(out)
}

fn module_per_base_n(out: &mut String, metrics: &crate::core::metrics::FinalMetrics) -> Result<()> {
    module_header(out, metrics.statuses.per_base_n, "Per base N content")?;
    module_desc(
        out,
        "Shows the proportion of Ns at each position. Spikes or elevated Ns suggest base-calling issues or low-complexity regions.",
    )?;
    let (w, h) = (800.0, 260.0);
    let data = metrics
        .per_base_n
        .iter()
        .map(|r| (r.base as f64, r.n_percent))
        .collect::<Vec<_>>();
    let (y_min, y_max) = auto_range(data.iter().map(|(_, y)| *y), 0.0, 100.0);
    svg_single_line_ybands(
        out,
        data.as_slice(),
        w,
        h,
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
    table_per_base_n(out, &metrics.per_base_n)?;
    module_footer(out)
}

fn module_per_seq_n(out: &mut String, metrics: &crate::core::metrics::FinalMetrics) -> Result<()> {
    module_header(out, metrics.statuses.per_seq_n, "Per sequence N content")?;
    module_desc(
        out,
        "Shows the distribution of N% per read. Excess high-N reads indicate poor base-calling or low-quality segments.",
    )?;
    let data = metrics
        .per_seq_n
        .iter()
        .map(|r| (r.n_percent as f64, r.count as f64))
        .collect::<Vec<_>>();
    svg_histogram_xbands(
        out,
        data.as_slice(),
        800.0,
        260.0,
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
    table_per_seq_n(out, &metrics.per_seq_n)?;
    module_footer(out)
}

fn module_length_dist_short(
    out: &mut String,
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<()> {
    module_header(
        out,
        metrics.statuses.length_dist,
        "Sequence Length Distribution",
    )?;
    module_desc(
        out,
        "Shows read length frequencies. Multiple peaks or long tails may indicate trimming or mixed read sources.",
    )?;
    let (w, h) = (800.0, 260.0);
    let data = metrics
        .length_dist
        .iter()
        .map(|r| (r.length as f64, r.count as f64))
        .collect::<Vec<_>>();
    svg_histogram(out, data.as_slice(), w, h, 0.0, 0.0, "Length", "Count")?;
    table_length_dist(out, &metrics.length_dist)?;
    module_footer(out)
}

fn module_length_dist_long(
    out: &mut String,
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<()> {
    module_header(
        out,
        metrics.statuses.length_dist,
        "Sequence Length Distribution",
    )?;
    module_desc(
        out,
        "Shows read length frequencies using log-scaled bins. Very long tails or multiple modes may indicate mixed input or variable trimming.",
    )?;
    if let Some(ref ll) = metrics.long_length {
        let data = ll
            .bins
            .iter()
            .enumerate()
            .map(|(i, &c)| (i as f64 + 1.0, c as f64))
            .collect::<Vec<_>>();
        svg_histogram(
            out,
            data.as_slice(),
            800.0,
            260.0,
            0.0,
            0.0,
            "Length bin",
            "Count",
        )?;
        table_long_length(out, ll)?;
    }
    module_footer(out)
}

fn module_duplication(
    out: &mut String,
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<()> {
    module_header(
        out,
        metrics.statuses.duplication,
        "Sequence Duplication Levels",
    )?;
    module_desc(
        out,
        "Estimates duplication using a streaming heavy-hitter model. High duplication often indicates PCR over-amplification or low library complexity.",
    )?;
    let data = metrics
        .duplication
        .iter()
        .enumerate()
        .map(|(i, r)| (i as f64 + 1.0, r.relative))
        .collect::<Vec<_>>();
    svg_histogram(
        out,
        data.as_slice(),
        800.0,
        260.0,
        0.0,
        0.0,
        "Level",
        "Relative count",
    )?;
    table_duplication(out, &metrics.duplication)?;
    module_footer(out)
}

fn module_overrep(out: &mut String, metrics: &crate::core::metrics::FinalMetrics) -> Result<()> {
    module_header(
        out,
        metrics.statuses.overrepresented,
        "Overrepresented sequences",
    )?;
    module_desc(
        out,
        "Lists sequences occurring more often than expected. Common sources are adapters, primers, or contamination.",
    )?;
    table_overrep(out, &metrics.overrepresented)?;
    module_footer(out)
}

fn module_adapter_content_short(
    out: &mut String,
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<()> {
    module_header(out, metrics.statuses.adapter_content, "Adapter Content")?;
    module_desc(
        out,
        "Shows adapter match percentages by position. Increasing signal toward read ends suggests adapter read-through.",
    )?;
    let (w, h) = (800.0, 260.0);
    svg_adapter_lines(out, &metrics.adapter_content, w, h, "Position", "%")?;
    table_adapter_content(out, &metrics.adapter_content)?;
    module_footer(out)
}

fn module_adapter_content_long(
    out: &mut String,
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<()> {
    module_header(out, metrics.statuses.adapter_content, "Adapter Content")?;
    module_desc(
        out,
        "Reports the fraction of reads containing common adapter motifs. Elevated percentages suggest residual adapters or chimeric reads.",
    )?;
    table_adapter_summary(out, &metrics.adapter_content)?;
    module_footer(out)
}

#[cfg(not(feature = "no-kmer"))]
fn module_kmer_content(
    out: &mut String,
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<()> {
    module_header(out, metrics.statuses.kmer_content, "Kmer Content")?;
    module_desc(
        out,
        "Reports k-mers enriched at specific positions. Strong enrichment can indicate adapters or sequence bias.",
    )?;
    table_kmer(out, &metrics.kmer_rows)?;
    module_footer(out)
}

fn svg_boxplot(
    out: &mut String,
    rows: &[crate::core::metrics::PerBaseQualRow],
    w: f64,
    h: f64,
    max_q: f64,
    x_label: &str,
    y_label: &str,
) -> Result<()> {
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
    // Background quality bands (FastQC-like).
    draw_y_bands(
        out,
        left,
        top,
        plot_w,
        plot_h,
        0.0,
        max_q,
        &[
            (0.0, 20.0, "#f4c7c3"),
            (20.0, 28.0, "#ffe5b4"),
            (28.0, max_q.max(28.0), "#cdeccf"),
        ],
    )?;

    let n = rows.len().max(1) as f64;
    let x_step = plot_w / n;
    let y_scale = if max_q <= 0.0 { 1.0 } else { plot_h / max_q };
    draw_y_axis_ticks(out, left, top, plot_w, plot_h, 0.0, max_q, 5)?;
    draw_y_axis_ticks_right(out, left, top, plot_w, plot_h, 0.0, max_q, 5)?;
    draw_x_axis_ticks(out, left, top, plot_w, plot_h, 1.0, n, 5)?;
    draw_axis_labels(out, left, top, plot_w, plot_h, x_label, y_label)?;

    for (i, r) in rows.iter().enumerate() {
        let x = left + (i as f64 + 0.5) * x_step;
        let y_m = top + plot_h - (r.median as f64 * y_scale);
        let y_lq = top + plot_h - (r.lower_quartile as f64 * y_scale);
        let y_uq = top + plot_h - (r.upper_quartile as f64 * y_scale);
        let y_p10 = top + plot_h - (r.p10 as f64 * y_scale);
        let y_p90 = top + plot_h - (r.p90 as f64 * y_scale);
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
        // Whiskers (10th-90th) and median line.
        writeln!(
            out,
            "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"#555\" stroke-width=\"1\"/>",
            x, y_p90, x, y_p10
        )?;
        let cap_w = (box_w * 0.6).max(1.0);
        writeln!(
            out,
            "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"#555\" stroke-width=\"1\"/>",
            x - cap_w / 2.0,
            y_p90,
            x + cap_w / 2.0,
            y_p90
        )?;
        writeln!(
            out,
            "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"#555\" stroke-width=\"1\"/>",
            x - cap_w / 2.0,
            y_p10,
            x + cap_w / 2.0,
            y_p10
        )?;
        writeln!(
            out,
            "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"#333\" stroke-width=\"1.2\"/>",
            box_x,
            y_m,
            box_x + box_w,
            y_m
        )?;
    }
    writeln!(out, "</svg></div>")?;
    Ok(())
}

fn svg_histogram(
    out: &mut String,
    data: &[(f64, f64)],
    w: f64,
    h: f64,
    min_x: f64,
    max_x: f64,
    x_label: &str,
    y_label: &str,
) -> Result<()> {
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
    let max_y = data.iter().map(|(_, y)| *y).fold(0.0, f64::max);
    let (x_min, x_max) = if min_x == max_x {
        let min_b = data.first().map(|d| d.0).unwrap_or(0.0);
        let max_b = data.last().map(|d| d.0).unwrap_or(1.0);
        auto_range(data.iter().map(|(x, _)| *x), min_b, max_b)
    } else {
        (min_x, max_x)
    };
    let bar_w = if data.is_empty() {
        1.0
    } else {
        plot_w / data.len() as f64
    };
    draw_y_axis_ticks(out, left, top, plot_w, plot_h, 0.0, max_y, 4)?;
    draw_y_axis_ticks_right(out, left, top, plot_w, plot_h, 0.0, max_y, 4)?;
    draw_x_axis_ticks(out, left, top, plot_w, plot_h, x_min, x_max, 5)?;
    draw_axis_labels(out, left, top, plot_w, plot_h, x_label, y_label)?;
    for (i, (_xv, yv)) in data.iter().enumerate() {
        let x = left + (i as f64) * bar_w;
        let y = if max_y == 0.0 {
            0.0
        } else {
            yv / max_y * plot_h
        };
        let y0 = top + plot_h - y;
        writeln!(
            out,
            "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"#7db8da\"/>",
            x,
            y0,
            bar_w.max(1.0),
            y
        )?;
    }
    writeln!(out, "</svg></div>")?;
    Ok(())
}

fn svg_histogram_compat_bars(
    out: &mut String,
    data: &[(f64, f64)],
    w: f64,
    h: f64,
    min_x: f64,
    max_x: f64,
    x_label: &str,
    y_label: &str,
) -> Result<()> {
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
    let max_y = data.iter().map(|(_, y)| *y).fold(0.0, f64::max);
    let (x_min, x_max) = if min_x == max_x {
        let min_b = data.first().map(|d| d.0).unwrap_or(0.0);
        let max_b = data.last().map(|d| d.0).unwrap_or(1.0);
        auto_range(data.iter().map(|(x, _)| *x), min_b, max_b)
    } else {
        (min_x, max_x)
    };
    let bar_w = if data.is_empty() {
        1.0
    } else {
        (plot_w / data.len() as f64).max(1.0)
    };
    draw_y_axis_labels_only(out, left, top, plot_w, plot_h, 0.0, max_y, 4)?;
    draw_x_axis_labels_only(out, left, top, plot_w, plot_h, x_min, x_max, 5)?;
    draw_axis_labels(out, left, top, plot_w, plot_h, x_label, y_label)?;
    for (i, (_xv, yv)) in data.iter().enumerate() {
        let x = left + (i as f64) * bar_w;
        let y = if max_y == 0.0 {
            0.0
        } else {
            yv / max_y * plot_h
        };
        let y0 = top + plot_h - y;
        writeln!(
            out,
            "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"#8ecae6\"/>",
            x, y0, bar_w, y
        )?;
    }
    writeln!(out, "</svg></div>")?;
    Ok(())
}

fn svg_histogram_xbands(
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
    let (x_min, x_max) = if min_x == max_x {
        let min_b = data.first().map(|d| d.0).unwrap_or(0.0);
        let max_b = data.last().map(|d| d.0).unwrap_or(1.0);
        auto_range(data.iter().map(|(x, _)| *x), min_b, max_b)
    } else {
        (min_x, max_x)
    };
    let x_range = (x_max - x_min).max(1.0);
    for (lo, hi, color) in bands {
        let start = ((*lo - x_min) / x_range).clamp(0.0, 1.0);
        let end = ((*hi - x_min) / x_range).clamp(0.0, 1.0);
        let x = left + start * plot_w;
        let w_band = (end - start) * plot_w;
        if w_band > 0.0 {
            writeln!(
                out,
                "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"{}\" opacity=\"0.18\"/>",
                x, top, w_band, plot_h, color
            )?;
        }
    }
    let max_y = data.iter().map(|(_, y)| *y).fold(0.0, f64::max);
    let bar_w = if data.is_empty() {
        1.0
    } else {
        plot_w / data.len() as f64
    };
    draw_y_axis_ticks(out, left, top, plot_w, plot_h, 0.0, max_y, 4)?;
    draw_y_axis_ticks_right(out, left, top, plot_w, plot_h, 0.0, max_y, 4)?;
    draw_x_axis_ticks(out, left, top, plot_w, plot_h, x_min, x_max, 5)?;
    draw_axis_labels(out, left, top, plot_w, plot_h, x_label, y_label)?;
    for (i, (_xv, yv)) in data.iter().enumerate() {
        let x = left + (i as f64) * bar_w;
        let y = if max_y == 0.0 {
            0.0
        } else {
            yv / max_y * plot_h
        };
        let y0 = top + plot_h - y;
        writeln!(
            out,
            "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"#7db8da\"/>",
            x,
            y0,
            bar_w.max(1.0),
            y
        )?;
    }
    writeln!(out, "</svg></div>")?;
    Ok(())
}

fn svg_multi_line(
    out: &mut String,
    rows: &[crate::core::metrics::PerBaseContentRow],
    w: f64,
    h: f64,
    x_label: &str,
    y_label: &str,
) -> Result<()> {
    let data_g = rows
        .iter()
        .map(|r| (r.base as f64, r.g))
        .collect::<Vec<_>>();
    let data_a = rows
        .iter()
        .map(|r| (r.base as f64, r.a))
        .collect::<Vec<_>>();
    let data_t = rows
        .iter()
        .map(|r| (r.base as f64, r.t))
        .collect::<Vec<_>>();
    let data_c = rows
        .iter()
        .map(|r| (r.base as f64, r.c))
        .collect::<Vec<_>>();
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
    let mut min_y: f64 = 100.0;
    let mut max_y: f64 = 0.0;
    for r in rows {
        min_y = min_y.min(r.a.min(r.c.min(r.g.min(r.t))));
        max_y = max_y.max(r.a.max(r.c.max(r.g.max(r.t))));
    }
    let pad = ((max_y - min_y) * 0.2).max(2.0);
    let y_min = (min_y - pad).max(0.0);
    let y_max = (max_y + pad).min(100.0);
    draw_y_axis_ticks(out, left, top, plot_w, plot_h, y_min, y_max, 5)?;
    draw_y_axis_ticks_right(out, left, top, plot_w, plot_h, y_min, y_max, 5)?;
    draw_x_axis_ticks(out, left, top, plot_w, plot_h, 1.0, rows.len() as f64, 5)?;
    draw_axis_labels(out, left, top, plot_w, plot_h, x_label, y_label)?;
    // FastQC line colours (Tol scheme): #882255, #332288, #117733, #DDCC77
    svg_line(
        out, &data_g, left, top, plot_w, plot_h, y_min, y_max, "#882255",
    )?;
    svg_line(
        out, &data_a, left, top, plot_w, plot_h, y_min, y_max, "#332288",
    )?;
    svg_line(
        out, &data_t, left, top, plot_w, plot_h, y_min, y_max, "#117733",
    )?;
    svg_line(
        out, &data_c, left, top, plot_w, plot_h, y_min, y_max, "#DDCC77",
    )?;
    writeln!(out, "</svg></div>")?;
    Ok(())
}

fn svg_single_line(
    out: &mut String,
    data: &[(f64, f64)],
    w: f64,
    h: f64,
    min_y: f64,
    max_y: f64,
    color: &str,
    x_label: &str,
    y_label: &str,
) -> Result<()> {
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
    draw_y_axis_ticks(out, left, top, plot_w, plot_h, min_y, max_y, 5)?;
    draw_y_axis_ticks_right(out, left, top, plot_w, plot_h, min_y, max_y, 5)?;
    let x_min = data.first().map(|d| d.0).unwrap_or(0.0);
    let x_max = data.last().map(|d| d.0).unwrap_or(1.0);
    draw_x_axis_ticks(out, left, top, plot_w, plot_h, x_min, x_max, 5)?;
    draw_axis_labels(out, left, top, plot_w, plot_h, x_label, y_label)?;
    svg_line(out, data, left, top, plot_w, plot_h, min_y, max_y, color)?;
    writeln!(out, "</svg></div>")?;
    Ok(())
}

fn svg_single_line_ybands(
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
    draw_y_bands(out, left, top, plot_w, plot_h, min_y, max_y, bands)?;
    draw_y_axis_ticks(out, left, top, plot_w, plot_h, min_y, max_y, 5)?;
    draw_y_axis_ticks_right(out, left, top, plot_w, plot_h, min_y, max_y, 5)?;
    let x_min = data.first().map(|d| d.0).unwrap_or(0.0);
    let x_max = data.last().map(|d| d.0).unwrap_or(1.0);
    draw_x_axis_ticks(out, left, top, plot_w, plot_h, x_min, x_max, 5)?;
    draw_axis_labels(out, left, top, plot_w, plot_h, x_label, y_label)?;
    svg_line(out, data, left, top, plot_w, plot_h, min_y, max_y, color)?;
    writeln!(out, "</svg></div>")?;
    Ok(())
}

fn draw_y_axis_ticks(
    out: &mut String,
    left: f64,
    top: f64,
    plot_w: f64,
    plot_h: f64,
    min_y: f64,
    max_y: f64,
    ticks: usize,
) -> Result<()> {
    if ticks < 2 || (max_y - min_y).abs() < 1e-9 {
        return Ok(());
    }
    let (start, step, count) = nice_ticks(min_y, max_y, ticks);
    for i in 0..count {
        let v = start + step * i as f64;
        let y = top + plot_h - ((v - min_y) / (max_y - min_y).max(1e-6)) * plot_h;
        writeln!(
            out,
            "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"#eee\"/>",
            left,
            y,
            left + plot_w,
            y
        )?;
        writeln!(
            out,
            "<text x=\"{}\" y=\"{}\" font-size=\"10\" fill=\"#666\" text-anchor=\"end\" dominant-baseline=\"middle\">{}</text>",
            left - 4.0,
            y,
            fmt_tick(v)
        )?;
    }
    Ok(())
}

fn draw_y_axis_labels_only(
    out: &mut String,
    left: f64,
    top: f64,
    _plot_w: f64,
    plot_h: f64,
    min_y: f64,
    max_y: f64,
    ticks: usize,
) -> Result<()> {
    if ticks < 2 || (max_y - min_y).abs() < 1e-9 {
        return Ok(());
    }
    let (start, step, count) = nice_ticks(min_y, max_y, ticks);
    for i in 0..count {
        let v = start + step * i as f64;
        let y = top + plot_h - ((v - min_y) / (max_y - min_y).max(1e-6)) * plot_h;
        writeln!(
            out,
            "<text x=\"{}\" y=\"{}\" font-size=\"10\" fill=\"#666\" text-anchor=\"end\" dominant-baseline=\"middle\">{}</text>",
            left - 4.0,
            y,
            fmt_tick(v)
        )?;
    }
    Ok(())
}

fn draw_y_axis_ticks_right(
    out: &mut String,
    left: f64,
    top: f64,
    plot_w: f64,
    plot_h: f64,
    min_y: f64,
    max_y: f64,
    ticks: usize,
) -> Result<()> {
    if ticks < 2 || (max_y - min_y).abs() < 1e-9 {
        return Ok(());
    }
    let (start, step, count) = nice_ticks(min_y, max_y, ticks);
    for i in 0..count {
        let v = start + step * i as f64;
        let y = top + plot_h - ((v - min_y) / (max_y - min_y).max(1e-6)) * plot_h;
        writeln!(
            out,
            "<text x=\"{}\" y=\"{}\" font-size=\"10\" fill=\"#666\" text-anchor=\"start\" dominant-baseline=\"middle\">{}</text>",
            left + plot_w + 4.0,
            y,
            fmt_tick(v)
        )?;
    }
    Ok(())
}

fn draw_x_axis_ticks(
    out: &mut String,
    left: f64,
    top: f64,
    plot_w: f64,
    plot_h: f64,
    min_x: f64,
    max_x: f64,
    ticks: usize,
) -> Result<()> {
    if ticks < 2 || (max_x - min_x).abs() < 1e-9 {
        return Ok(());
    }
    let (start, step, count) = nice_ticks(min_x, max_x, ticks);
    for i in 0..count {
        let v = start + step * i as f64;
        let x = left + ((v - min_x) / (max_x - min_x).max(1e-6)) * plot_w;
        writeln!(
            out,
            "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"#eee\"/>",
            x,
            top,
            x,
            top + plot_h
        )?;
        writeln!(
            out,
            "<text x=\"{}\" y=\"{}\" font-size=\"10\" fill=\"#666\" text-anchor=\"middle\" dominant-baseline=\"hanging\">{}</text>",
            x,
            top + plot_h + 4.0,
            fmt_tick(v)
        )?;
    }
    Ok(())
}

fn draw_x_axis_labels_only(
    out: &mut String,
    left: f64,
    top: f64,
    plot_w: f64,
    plot_h: f64,
    min_x: f64,
    max_x: f64,
    ticks: usize,
) -> Result<()> {
    if ticks < 2 || (max_x - min_x).abs() < 1e-9 {
        return Ok(());
    }
    let (start, step, count) = nice_ticks(min_x, max_x, ticks);
    for i in 0..count {
        let v = start + step * i as f64;
        let x = left + ((v - min_x) / (max_x - min_x).max(1e-6)) * plot_w;
        writeln!(
            out,
            "<text x=\"{}\" y=\"{}\" font-size=\"10\" fill=\"#666\" text-anchor=\"middle\" dominant-baseline=\"hanging\">{}</text>",
            x,
            top + plot_h + 4.0,
            fmt_tick(v)
        )?;
    }
    Ok(())
}

fn draw_axis_labels(
    out: &mut String,
    left: f64,
    top: f64,
    plot_w: f64,
    plot_h: f64,
    x_label: &str,
    y_label: &str,
) -> Result<()> {
    let x = left + plot_w / 2.0;
    let y = top + plot_h + 22.0;
    writeln!(
        out,
        "<text x=\"{}\" y=\"{}\" font-size=\"11\" fill=\"#444\" text-anchor=\"middle\">{}</text>",
        x, y, x_label
    )?;
    let yx = left - 28.0;
    let yy = top + plot_h / 2.0;
    writeln!(
        out,
        "<text x=\"{}\" y=\"{}\" font-size=\"11\" fill=\"#444\" text-anchor=\"middle\" transform=\"rotate(-90 {} {})\">{}</text>",
        yx, yy, yx, yy, y_label
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

fn fmt_int(v: u64) -> String {
    let s = v.to_string();
    let mut out = String::with_capacity(s.len() + s.len() / 3);
    let len = s.len();
    for (i, ch) in s.chars().enumerate() {
        if i != 0 && (len - i) % 3 == 0 {
            out.push(',');
        }
        out.push(ch);
    }
    out
}

fn fmt_timestamp(ts: u64) -> String {
    let days = (ts / 86_400) as i64;
    let secs = (ts % 86_400) as u32;
    let hour = secs / 3_600;
    let min = (secs % 3_600) / 60;
    let sec = secs % 60;

    let z = days + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = mp + if mp < 10 { 3 } else { -9 };
    let year = y + if m <= 2 { 1 } else { 0 };

    format!(
        "{:04}-{:02}-{:02} {:02}:{:02}:{:02} UTC",
        year, m, d, hour, min, sec
    )
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
    let count = ((end - start) / step).round() as usize + 1;
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

fn legend_base_content(out: &mut String) -> Result<()> {
    writeln!(
        out,
        "<div class=\"desc\"><b>Legend:</b> <span style=\"display:inline-block;width:18px;height:4px;background:#882255;margin:0 6px 2px 6px;vertical-align:middle;\"></span><b>G</b> <span style=\"display:inline-block;width:18px;height:4px;background:#332288;margin:0 6px 2px 10px;vertical-align:middle;\"></span><b>A</b> <span style=\"display:inline-block;width:18px;height:4px;background:#117733;margin:0 6px 2px 10px;vertical-align:middle;\"></span><b>T</b> <span style=\"display:inline-block;width:18px;height:4px;background:#DDCC77;margin:0 6px 2px 10px;vertical-align:middle;\"></span><b>C</b></div>"
    )?;
    Ok(())
}

fn draw_y_bands(
    out: &mut String,
    left: f64,
    top: f64,
    plot_w: f64,
    plot_h: f64,
    min_y: f64,
    max_y: f64,
    bands: &[(f64, f64, &str)],
) -> Result<()> {
    let y_range = (max_y - min_y).max(1.0);
    for (lo, hi, color) in bands {
        let start = ((*lo - min_y) / y_range).clamp(0.0, 1.0);
        let end = ((*hi - min_y) / y_range).clamp(0.0, 1.0);
        let y1 = top + plot_h - end * plot_h;
        let y2 = top + plot_h - start * plot_h;
        let h = (y2 - y1).max(0.0);
        if h > 0.0 {
            writeln!(
                out,
                "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"{}\" opacity=\"0.18\"/>",
                left, y1, plot_w, h, color
            )?;
        }
    }
    Ok(())
}

fn svg_line(
    out: &mut String,
    data: &[(f64, f64)],
    left: f64,
    top: f64,
    plot_w: f64,
    plot_h: f64,
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
    let y_range = (max_y - min_y).max(1.0);

    let mut path = String::new();
    for (i, (xv, yv)) in data.iter().enumerate() {
        let x = left + (*xv - x_min) / x_range * plot_w;
        let y = top + plot_h - ((*yv - min_y) / y_range * plot_h);
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

fn table_per_base_quality(
    out: &mut String,
    rows: &[crate::core::metrics::PerBaseQualRow],
) -> Result<()> {
    writeln!(
        out,
        "<details><summary>Table</summary><table class=\"table\">"
    )?;
    writeln!(
        out,
        "<tr><th>Base</th><th>Mean</th><th>Median</th><th>Lower Quartile</th><th>Upper Quartile</th><th>10th Percentile</th><th>90th Percentile</th></tr>"
    )?;
    for r in rows {
        writeln!(
            out,
            "<tr><td>{}</td><td>{:.1}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>",
            r.base, r.mean, r.median, r.lower_quartile, r.upper_quartile, r.p10, r.p90
        )?;
    }
    writeln!(out, "</table></details>")?;
    Ok(())
}

fn table_per_seq_quality(
    out: &mut String,
    rows: &[crate::core::metrics::PerSeqQualRow],
) -> Result<()> {
    writeln!(
        out,
        "<details><summary>Table</summary><table class=\"table\">"
    )?;
    writeln!(out, "<tr><th>Quality</th><th>Count</th></tr>")?;
    for r in rows {
        writeln!(out, "<tr><td>{}</td><td>{}</td></tr>", r.mean_q, r.count)?;
    }
    writeln!(out, "</table></details>")?;
    Ok(())
}

fn table_per_seq_n(out: &mut String, rows: &[crate::core::metrics::PerSeqNRow]) -> Result<()> {
    writeln!(
        out,
        "<details><summary>Table</summary><table class=\"table\">"
    )?;
    writeln!(out, "<tr><th>N%</th><th>Count</th></tr>")?;
    for r in rows {
        writeln!(out, "<tr><td>{}</td><td>{}</td></tr>", r.n_percent, r.count)?;
    }
    writeln!(out, "</table></details>")?;
    Ok(())
}

fn table_per_base_content(
    out: &mut String,
    rows: &[crate::core::metrics::PerBaseContentRow],
) -> Result<()> {
    writeln!(
        out,
        "<details><summary>Table</summary><table class=\"table\">"
    )?;
    writeln!(
        out,
        "<tr><th>Base</th><th>G</th><th>A</th><th>T</th><th>C</th></tr>"
    )?;
    for r in rows {
        writeln!(
            out,
            "<tr><td>{}</td><td>{:.1}</td><td>{:.1}</td><td>{:.1}</td><td>{:.1}</td></tr>",
            r.base, r.g, r.a, r.t, r.c
        )?;
    }
    writeln!(out, "</table></details>")?;
    Ok(())
}

fn table_per_seq_gc(out: &mut String, rows: &[crate::core::metrics::PerSeqGcRow]) -> Result<()> {
    writeln!(
        out,
        "<details><summary>Table</summary><table class=\"table\">"
    )?;
    writeln!(out, "<tr><th>GC%</th><th>Count</th></tr>")?;
    for r in rows {
        writeln!(out, "<tr><td>{}</td><td>{}</td></tr>", r.gc, r.count)?;
    }
    writeln!(out, "</table></details>")?;
    Ok(())
}

fn table_per_base_n(out: &mut String, rows: &[crate::core::metrics::PerBaseNRow]) -> Result<()> {
    writeln!(
        out,
        "<details><summary>Table</summary><table class=\"table\">"
    )?;
    writeln!(out, "<tr><th>Base</th><th>N%</th></tr>")?;
    for r in rows {
        writeln!(
            out,
            "<tr><td>{}</td><td>{:.1}</td></tr>",
            r.base, r.n_percent
        )?;
    }
    writeln!(out, "</table></details>")?;
    Ok(())
}

fn table_length_dist(out: &mut String, rows: &[crate::core::metrics::LengthDistRow]) -> Result<()> {
    writeln!(
        out,
        "<details><summary>Table</summary><table class=\"table\">"
    )?;
    writeln!(out, "<tr><th>Length</th><th>Count</th></tr>")?;
    for r in rows {
        writeln!(out, "<tr><td>{}</td><td>{}</td></tr>", r.length, r.count)?;
    }
    writeln!(out, "</table></details>")?;
    Ok(())
}

fn table_long_length(out: &mut String, ll: &crate::core::metrics::LongLengthSummary) -> Result<()> {
    writeln!(
        out,
        "<details><summary>Table</summary><table class=\"table\">"
    )?;
    writeln!(out, "<tr><th>Metric</th><th>Value</th></tr>")?;
    writeln!(out, "<tr><td>Min</td><td>{}</td></tr>", ll.min)?;
    writeln!(out, "<tr><td>Max</td><td>{}</td></tr>", ll.max)?;
    writeln!(out, "<tr><td>Mean</td><td>{:.1}</td></tr>", ll.mean)?;
    writeln!(out, "<tr><td>N50</td><td>{}</td></tr>", ll.n50)?;
    writeln!(out, "<tr><td>N90</td><td>{}</td></tr>", ll.n90)?;
    writeln!(out, "</table>")?;
    writeln!(out, "<table class=\"table\">")?;
    writeln!(out, "<tr><th>Length Bin</th><th>Count</th></tr>")?;
    for i in 0..ll.bins.len() {
        writeln!(
            out,
            "<tr><td>{}</td><td>{}</td></tr>",
            ll.labels[i], ll.bins[i]
        )?;
    }
    writeln!(out, "</table></details>")?;
    Ok(())
}

fn table_duplication(
    out: &mut String,
    rows: &[crate::core::metrics::DuplicationRow],
) -> Result<()> {
    writeln!(
        out,
        "<details><summary>Table</summary><table class=\"table\">"
    )?;
    writeln!(
        out,
        "<tr><th>Duplication Level</th><th>Relative Count</th></tr>"
    )?;
    for r in rows {
        writeln!(
            out,
            "<tr><td>{}</td><td>{:.2}</td></tr>",
            r.level.as_str(),
            r.relative
        )?;
    }
    writeln!(out, "</table></details>")?;
    Ok(())
}

fn table_overrep(out: &mut String, rows: &[crate::core::metrics::OverrepRow]) -> Result<()> {
    writeln!(
        out,
        "<details><summary>Table</summary><table class=\"table\">"
    )?;
    writeln!(
        out,
        "<tr><th>Sequence</th><th>Count</th><th>Percentage</th><th>Possible Source</th></tr>"
    )?;
    for r in rows {
        writeln!(
            out,
            "<tr><td>{}</td><td>{}</td><td>{:.2}</td><td>{}</td></tr>",
            r.sequence, r.count, r.percent, r.source
        )?;
    }
    writeln!(out, "</table></details>")?;
    Ok(())
}

fn table_adapter_content(
    out: &mut String,
    rows: &[crate::core::metrics::AdapterRow],
) -> Result<()> {
    writeln!(
        out,
        "<details><summary>Table</summary><table class=\"table\">"
    )?;
    write!(out, "<tr><th>Position</th>")?;
    for name in crate::core::metrics::ADAPTERS {
        write!(out, "<th>{}</th>", name)?;
    }
    writeln!(out, "</tr>")?;
    for r in rows {
        write!(out, "<tr><td>{}</td>", r.position)?;
        for v in r.values.iter() {
            write!(out, "<td>{:.1}</td>", v)?;
        }
        writeln!(out, "</tr>")?;
    }
    writeln!(out, "</table></details>")?;
    Ok(())
}

fn table_adapter_summary(
    out: &mut String,
    rows: &[crate::core::metrics::AdapterRow],
) -> Result<()> {
    writeln!(
        out,
        "<details><summary>Table</summary><table class=\"table\">"
    )?;
    write!(out, "<tr><th>Adapter</th>")?;
    for name in crate::core::metrics::ADAPTERS {
        write!(out, "<th>{}</th>", name)?;
    }
    writeln!(out, "</tr>")?;
    if let Some(r) = rows.first() {
        write!(out, "<tr><td>Any</td>")?;
        for v in r.values.iter() {
            write!(out, "<td>{:.1}</td>", v)?;
        }
        writeln!(out, "</tr>")?;
    }
    writeln!(out, "</table></details>")?;
    Ok(())
}

#[cfg(not(feature = "no-kmer"))]
fn table_kmer(out: &mut String, rows: &[crate::core::metrics::KmerRow]) -> Result<()> {
    writeln!(
        out,
        "<details open><summary>Table</summary><table class=\"table sortable\" data-sortable=\"true\">"
    )?;
    writeln!(
        out,
        "<tr><th>Sequence</th><th>Count</th><th>PValue</th><th>Obs/Exp Max</th><th>Max Obs/Exp Position</th></tr>"
    )?;
    for r in rows {
        writeln!(
            out,
            "<tr><td>{}</td><td>{}</td><td>{:.2e}</td><td>{:.2}</td><td>{}</td></tr>",
            r.sequence, r.count, r.p_value, r.obs_exp, r.max_pos
        )?;
    }
    writeln!(out, "</table></details>")?;
    Ok(())
}

pub(crate) fn latex_svg_per_base_quality(
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<String> {
    let (w, h) = (800.0, 260.0);
    let max_q = metrics
        .per_base_qual
        .iter()
        .map(|r| r.p90 as f64)
        .fold(40.0, f64::max);
    let mut s = String::new();
    svg_boxplot(
        &mut s,
        &metrics.per_base_qual,
        w,
        h,
        max_q,
        "Position",
        "Quality",
    )?;
    Ok(extract_svg(&s))
}

pub(crate) fn latex_svg_per_seq_quality(
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<String> {
    let (w, h) = (800.0, 260.0);
    let data = metrics
        .per_seq_qual
        .iter()
        .map(|r| (r.mean_q as f64, r.count as f64))
        .collect::<Vec<_>>();
    let mut s = String::new();
    svg_histogram_compat_bars(&mut s, data.as_slice(), w, h, 0.0, 0.0, "Mean Q", "Count")?;
    Ok(extract_svg(&s))
}

pub(crate) fn latex_svg_duplication(
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<String> {
    let data = metrics
        .duplication
        .iter()
        .enumerate()
        .map(|(i, r)| (i as f64 + 1.0, r.relative))
        .collect::<Vec<_>>();
    let mut s = String::new();
    svg_histogram_compat_bars(
        &mut s,
        data.as_slice(),
        800.0,
        260.0,
        0.0,
        0.0,
        "Level",
        "Relative count",
    )?;
    Ok(extract_svg(&s))
}

pub(crate) fn latex_svg_adapter_content(
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<String> {
    let mut s = String::new();
    svg_adapter_lines(
        &mut s,
        &metrics.adapter_content,
        800.0,
        260.0,
        "Position",
        "%",
    )?;
    Ok(extract_svg(&s))
}

pub(crate) fn latex_svg_per_base_content(
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<String> {
    let mut s = String::new();
    svg_multi_line(
        &mut s,
        &metrics.per_base_content,
        800.0,
        260.0,
        "Position",
        "%",
    )?;
    Ok(extract_svg(&s))
}

pub(crate) fn latex_svg_per_seq_gc(metrics: &crate::core::metrics::FinalMetrics) -> Result<String> {
    let data = metrics
        .per_seq_gc
        .iter()
        .map(|r| (r.gc as f64, r.count as f64))
        .collect::<Vec<_>>();
    let mut s = String::new();
    svg_histogram_compat_bars(
        &mut s,
        data.as_slice(),
        800.0,
        260.0,
        0.0,
        100.0,
        "GC%",
        "Count",
    )?;
    Ok(extract_svg(&s))
}

pub(crate) fn latex_svg_per_base_n(metrics: &crate::core::metrics::FinalMetrics) -> Result<String> {
    let data = metrics
        .per_base_n
        .iter()
        .map(|r| (r.base as f64, r.n_percent))
        .collect::<Vec<_>>();
    let (y_min, y_max) = auto_range(data.iter().map(|(_, y)| *y), 0.0, 100.0);
    let mut s = String::new();
    svg_single_line_ybands(
        &mut s,
        data.as_slice(),
        800.0,
        260.0,
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

pub(crate) fn latex_svg_per_seq_n(metrics: &crate::core::metrics::FinalMetrics) -> Result<String> {
    let data = metrics
        .per_seq_n
        .iter()
        .map(|r| (r.n_percent as f64, r.count as f64))
        .collect::<Vec<_>>();
    let mut s = String::new();
    svg_histogram_compat_bars(
        &mut s,
        data.as_slice(),
        800.0,
        260.0,
        0.0,
        100.0,
        "N%",
        "Count",
    )?;
    Ok(extract_svg(&s))
}

pub(crate) fn latex_svg_length_dist(
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<String> {
    let mut s = String::new();
    if let Some(ref ll) = metrics.long_length {
        let data = ll
            .bins
            .iter()
            .enumerate()
            .map(|(i, &c)| (i as f64 + 1.0, c as f64))
            .collect::<Vec<_>>();
        svg_histogram_compat_bars(
            &mut s,
            data.as_slice(),
            800.0,
            260.0,
            0.0,
            0.0,
            "Length bin",
            "Count",
        )?;
    } else {
        let data = metrics
            .length_dist
            .iter()
            .map(|r| (r.length as f64, r.count as f64))
            .collect::<Vec<_>>();
        svg_histogram_compat_bars(
            &mut s,
            data.as_slice(),
            800.0,
            260.0,
            0.0,
            0.0,
            "Length",
            "Count",
        )?;
    }
    Ok(extract_svg(&s))
}

pub(crate) fn latex_svg_overrep(metrics: &crate::core::metrics::FinalMetrics) -> Result<String> {
    let mut lines = Vec::new();
    for r in metrics.overrepresented.iter().take(6) {
        lines.push(format!("{} ({:.2}%)", r.sequence, r.percent));
    }
    if lines.is_empty() {
        lines.push("No overrepresented sequences".to_string());
    }
    Ok(simple_text_svg("Overrepresented sequences", &lines))
}

#[cfg(not(feature = "no-kmer"))]
pub(crate) fn latex_svg_kmer_content(
    metrics: &crate::core::metrics::FinalMetrics,
) -> Result<String> {
    let mut lines = Vec::new();
    for r in metrics.kmer_rows.iter().take(6) {
        lines.push(format!("{} (obs/exp {:.2})", r.sequence, r.obs_exp));
    }
    if lines.is_empty() {
        lines.push("No enriched k-mers".to_string());
    }
    Ok(simple_text_svg("Kmer content", &lines))
}

fn extract_svg(s: &str) -> String {
    if let (Some(start), Some(end)) = (s.find("<svg"), s.rfind("</svg>")) {
        s[start..end + 6].to_string()
    } else {
        s.to_string()
    }
}

fn simple_text_svg(title: &str, lines: &[String]) -> String {
    let mut out = String::new();
    let w = 800;
    let h = 260;
    out.push_str(&format!(
        "<svg width=\"{}\" height=\"{}\" viewBox=\"0 0 {} {}\">",
        w, h, w, h
    ));
    out.push_str(
        "<rect x=\"0\" y=\"0\" width=\"800\" height=\"260\" fill=\"#fff\" stroke=\"#ddd\"/>",
    );
    out.push_str(&format!(
        "<text x=\"16\" y=\"28\" font-size=\"14\" fill=\"#333\">{}</text>",
        title
    ));
    let mut y = 54;
    for l in lines {
        out.push_str(&format!(
            "<text x=\"16\" y=\"{}\" font-size=\"12\" fill=\"#333\">{}</text>",
            y,
            escape_svg(l)
        ));
        y += 18;
    }
    out.push_str("</svg>");
    out
}

fn escape_svg(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

fn svg_adapter_lines(
    out: &mut String,
    rows: &[crate::core::metrics::AdapterRow],
    w: f64,
    h: f64,
    x_label: &str,
    y_label: &str,
) -> Result<()> {
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
    let mut series: Vec<Vec<(f64, f64)>> = vec![Vec::new(); crate::core::metrics::ADAPTERS.len()];
    for r in rows {
        for i in 0..series.len() {
            series[i].push((r.position as f64, r.values[i]));
        }
    }
    let (y_min, y_max) = auto_range(
        series.iter().flat_map(|s| s.iter().map(|(_, y)| *y)),
        0.0,
        100.0,
    );
    draw_y_axis_ticks(out, left, top, plot_w, plot_h, y_min, y_max, 5)?;
    draw_y_axis_ticks_right(out, left, top, plot_w, plot_h, y_min, y_max, 5)?;
    let x_min = series
        .first()
        .and_then(|s| s.first())
        .map(|d| d.0)
        .unwrap_or(0.0);
    let x_max = series
        .first()
        .and_then(|s| s.last())
        .map(|d| d.0)
        .unwrap_or(1.0);
    draw_x_axis_ticks(out, left, top, plot_w, plot_h, x_min, x_max, 5)?;
    draw_axis_labels(out, left, top, plot_w, plot_h, x_label, y_label)?;
    let colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"];
    for i in 0..series.len() {
        let color = colors[i % colors.len()];
        svg_line(
            out, &series[i], left, top, plot_w, plot_h, y_min, y_max, color,
        )?;
    }
    writeln!(out, "</svg></div>")?;
    Ok(())
}
