//! Machine-readable report.
//!
//! `fastqc_data.txt` is a legacy format that every consumer has to parse with
//! bespoke code. This emits the same numbers as JSON so pipelines can read them
//! directly.

use crate::core::metrics::FinalMetrics;
use crate::core::model::{FinalizeContext, Mode};
use crate::report::{ModuleId, active_modules, format_timestamp, report_timestamp};
use crate::{FASTQC_FORMAT_VERSION, KIRA_QC_VERSION};
use anyhow::{Context, Result};
use std::fmt::Write as FmtWrite;
use std::fs;
use std::path::Path;

pub fn write(path: &Path, metrics: &FinalMetrics, ctx: &FinalizeContext) -> Result<()> {
    let text = render(metrics, ctx);
    fs::write(path, text).with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}

pub fn render(metrics: &FinalMetrics, ctx: &FinalizeContext) -> String {
    let mut o = String::with_capacity(64 * 1024);
    let ts = report_timestamp();

    o.push_str("{\n");
    kv_str(&mut o, 1, "tool", "kira-qc", true);
    kv_str(&mut o, 1, "version", KIRA_QC_VERSION, true);
    kv_str(
        &mut o,
        1,
        "fastqc_format_version",
        FASTQC_FORMAT_VERSION,
        true,
    );
    kv_str(&mut o, 1, "generated_at", &format_timestamp(ts), true);
    kv_str(&mut o, 1, "filename", &ctx.file_name, true);
    kv_str(&mut o, 1, "sample_name", &ctx.sample_name, true);
    kv_str(
        &mut o,
        1,
        "mode",
        match ctx.mode {
            Mode::Short => "short",
            Mode::Long => "long",
        },
        true,
    );

    // Statuses, keyed by module anchor so they are stable identifiers.
    let modules = active_modules(metrics, ctx);
    indent(&mut o, 1);
    o.push_str("\"statuses\": {\n");
    for (i, m) in modules.iter().enumerate() {
        indent(&mut o, 2);
        let _ = write!(
            o,
            "\"{}\": \"{}\"",
            m.anchor(),
            m.status(metrics).as_str_lower()
        );
        o.push_str(if i + 1 < modules.len() { ",\n" } else { "\n" });
    }
    indent(&mut o, 1);
    o.push_str("},\n");

    // Basic statistics.
    let b = &metrics.basic;
    indent(&mut o, 1);
    o.push_str("\"basic_statistics\": {\n");
    kv_str(&mut o, 2, "file_type", b.file_type, true);
    kv_str(&mut o, 2, "encoding", b.encoding, true);
    kv_num(&mut o, 2, "total_sequences", b.total_sequences as f64, true);
    kv_num(&mut o, 2, "total_bases", b.total_bases as f64, true);
    kv_num(
        &mut o,
        2,
        "sequences_flagged_as_poor_quality",
        b.poor_quality_sequences as f64,
        true,
    );
    kv_num(&mut o, 2, "min_length", b.min_len as f64, true);
    kv_num(&mut o, 2, "max_length", b.max_len as f64, true);
    kv_num(&mut o, 2, "percent_gc", b.gc_percent as f64, false);
    o.push('\n');
    indent(&mut o, 1);
    o.push_str("},\n");

    if modules.contains(&ModuleId::PerBaseQual) {
        array_of_objects(
            &mut o,
            1,
            "per_base_quality",
            &metrics.per_base_qual,
            |o, r| {
                kv_str(o, 3, "base", &r.group.label(), true);
                kv_num(o, 3, "mean", r.mean, true);
                kv_num(o, 3, "median", r.median, true);
                kv_num(o, 3, "lower_quartile", r.lower_quartile, true);
                kv_num(o, 3, "upper_quartile", r.upper_quartile, true);
                kv_num(o, 3, "percentile_10", r.p10, true);
                kv_num(o, 3, "percentile_90", r.p90, false);
            },
        );
    }

    if modules.contains(&ModuleId::PerTile) {
        array_of_objects(&mut o, 1, "per_tile_quality", &metrics.per_tile, |o, r| {
            kv_num(o, 3, "tile", r.tile as f64, true);
            kv_str(o, 3, "base", &r.group.label(), true);
            kv_num(o, 3, "mean_deviation", r.mean_deviation, false);
        });
    }

    array_of_objects(
        &mut o,
        1,
        "per_sequence_quality",
        &metrics.per_seq_qual,
        |o, r| {
            kv_num(o, 3, "quality", r.mean_q as f64, true);
            kv_num(o, 3, "count", r.count as f64, false);
        },
    );

    if modules.contains(&ModuleId::PerBaseContent) {
        array_of_objects(
            &mut o,
            1,
            "per_base_content",
            &metrics.per_base_content,
            |o, r| {
                kv_str(o, 3, "base", &r.group.label(), true);
                kv_num(o, 3, "g", r.g, true);
                kv_num(o, 3, "a", r.a, true);
                kv_num(o, 3, "t", r.t, true);
                kv_num(o, 3, "c", r.c, false);
            },
        );
    }

    array_of_objects(&mut o, 1, "per_sequence_gc", &metrics.per_seq_gc, |o, r| {
        kv_num(o, 3, "gc_percent", r.gc as f64, true);
        kv_num(o, 3, "count", r.count as f64, false);
    });

    if modules.contains(&ModuleId::PerBaseN) {
        array_of_objects(&mut o, 1, "per_base_n", &metrics.per_base_n, |o, r| {
            kv_str(o, 3, "base", &r.group.label(), true);
            kv_num(o, 3, "n_percent", r.n_percent, false);
        });
    }

    if modules.contains(&ModuleId::PerSeqN) {
        array_of_objects(&mut o, 1, "per_sequence_n", &metrics.per_seq_n, |o, r| {
            kv_num(o, 3, "n_percent", r.n_percent as f64, true);
            kv_num(o, 3, "count", r.count as f64, false);
        });
    }

    if let Some(ref ll) = metrics.long_length {
        indent(&mut o, 1);
        o.push_str("\"length_summary\": {\n");
        kv_num(&mut o, 2, "min", ll.min as f64, true);
        kv_num(&mut o, 2, "max", ll.max as f64, true);
        kv_num(&mut o, 2, "mean", ll.mean, true);
        kv_num(&mut o, 2, "median", ll.median as f64, true);
        kv_num(&mut o, 2, "n50", ll.n50 as f64, true);
        kv_num(&mut o, 2, "n90", ll.n90 as f64, true);
        kv_num(&mut o, 2, "total_bases", ll.total_bases as f64, true);
        indent(&mut o, 2);
        o.push_str("\"bins\": [\n");
        for (i, (label, count)) in ll.bins.iter().enumerate() {
            indent(&mut o, 3);
            let _ = write!(o, "{{\"range\": {}, \"count\": {}}}", quote(label), count);
            o.push_str(if i + 1 < ll.bins.len() { ",\n" } else { "\n" });
        }
        indent(&mut o, 2);
        o.push_str("]\n");
        indent(&mut o, 1);
        o.push_str("},\n");
    } else {
        array_of_objects(
            &mut o,
            1,
            "length_distribution",
            &metrics.length_dist,
            |o, r| {
                kv_num(o, 3, "length", r.length as f64, true);
                kv_num(o, 3, "count", r.count as f64, false);
            },
        );
    }

    if modules.contains(&ModuleId::Duplication) {
        indent(&mut o, 1);
        o.push_str("\"duplication\": {\n");
        let dedup = metrics
            .duplication_summary
            .as_ref()
            .map(|s| s.total_deduplicated_percent)
            .unwrap_or(100.0);
        kv_num(&mut o, 2, "total_deduplicated_percent", dedup, true);
        kv_num(
            &mut o,
            2,
            "reads_analysed",
            ctx.dup_window_reads as f64,
            true,
        );
        kv_bool(&mut o, 2, "whole_file", ctx.dup_window_complete, true);
        indent(&mut o, 2);
        o.push_str("\"levels\": [\n");
        for (i, r) in metrics.duplication.iter().enumerate() {
            indent(&mut o, 3);
            let _ = write!(
                o,
                "{{\"level\": {}, \"percent_of_deduplicated\": {}, \"percent_of_total\": {}}}",
                quote(r.level),
                num(r.percent_of_deduplicated),
                num(r.percent_of_total)
            );
            o.push_str(if i + 1 < metrics.duplication.len() {
                ",\n"
            } else {
                "\n"
            });
        }
        indent(&mut o, 2);
        o.push_str("]\n");
        indent(&mut o, 1);
        o.push_str("},\n");

        array_of_objects(
            &mut o,
            1,
            "overrepresented_sequences",
            &metrics.overrepresented,
            |o, r| {
                kv_str(o, 3, "sequence", &r.sequence, true);
                kv_num(o, 3, "count", r.count as f64, true);
                kv_num(o, 3, "percent", r.percent, true);
                kv_str(o, 3, "possible_source", &r.source, false);
            },
        );
    }

    // Adapter content: one object per position group, adapters keyed by name.
    indent(&mut o, 1);
    o.push_str("\"adapter_content\": [\n");
    for (i, row) in metrics.adapter_content.iter().enumerate() {
        indent(&mut o, 2);
        o.push_str("{\n");
        kv_str(&mut o, 3, "position", &row.group.label(), true);
        indent(&mut o, 3);
        o.push_str("\"adapters\": {");
        for (j, name) in metrics.adapter_names.iter().enumerate() {
            if j > 0 {
                o.push_str(", ");
            }
            let _ = write!(
                o,
                "{}: {}",
                quote(name),
                num(row.values.get(j).copied().unwrap_or(0.0))
            );
        }
        o.push_str("}\n");
        indent(&mut o, 2);
        o.push('}');
        o.push_str(if i + 1 < metrics.adapter_content.len() {
            ",\n"
        } else {
            "\n"
        });
    }
    indent(&mut o, 1);
    o.push(']');

    if modules.contains(&ModuleId::KmerContent) {
        o.push_str(",\n");
        array_of_objects_last(&mut o, 1, "kmer_content", &metrics.kmer_rows, |o, r| {
            kv_str(o, 3, "sequence", &r.sequence, true);
            kv_num(o, 3, "count", r.count as f64, true);
            kv_num(o, 3, "p_value", r.p_value, true);
            kv_num(o, 3, "obs_exp_max", r.obs_exp, true);
            kv_num(o, 3, "max_obs_exp_position", r.max_pos as f64, false);
        });
    } else {
        o.push('\n');
    }

    o.push_str("}\n");
    o
}

fn indent(o: &mut String, level: usize) {
    for _ in 0..level {
        o.push_str("  ");
    }
}

fn kv_str(o: &mut String, level: usize, key: &str, value: &str, comma: bool) {
    indent(o, level);
    let _ = write!(o, "\"{}\": {}", key, quote(value));
    if comma {
        o.push(',');
    }
    o.push('\n');
}

fn kv_num(o: &mut String, level: usize, key: &str, value: f64, comma: bool) {
    indent(o, level);
    let _ = write!(o, "\"{}\": {}", key, num(value));
    if comma {
        o.push(',');
    }
    o.push('\n');
}

fn kv_bool(o: &mut String, level: usize, key: &str, value: bool, comma: bool) {
    indent(o, level);
    let _ = write!(o, "\"{}\": {}", key, value);
    if comma {
        o.push(',');
    }
    o.push('\n');
}

fn array_of_objects<T, F>(o: &mut String, level: usize, key: &str, rows: &[T], mut body: F)
where
    F: FnMut(&mut String, &T),
{
    array_of_objects_inner(o, level, key, rows, &mut body);
    o.push_str(",\n");
}

fn array_of_objects_last<T, F>(o: &mut String, level: usize, key: &str, rows: &[T], mut body: F)
where
    F: FnMut(&mut String, &T),
{
    array_of_objects_inner(o, level, key, rows, &mut body);
    o.push('\n');
}

fn array_of_objects_inner<T, F>(o: &mut String, level: usize, key: &str, rows: &[T], body: &mut F)
where
    F: FnMut(&mut String, &T),
{
    indent(o, level);
    let _ = write!(o, "\"{}\": [", key);
    if rows.is_empty() {
        o.push(']');
        return;
    }
    o.push('\n');
    for (i, r) in rows.iter().enumerate() {
        indent(o, level + 1);
        o.push_str("{\n");
        body(o, r);
        indent(o, level + 1);
        o.push('}');
        if i + 1 < rows.len() {
            o.push(',');
        }
        o.push('\n');
    }
    indent(o, level);
    o.push(']');
}

/// JSON string literal with the escapes the spec requires.
pub fn quote(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{08}' => out.push_str("\\b"),
            '\u{0c}' => out.push_str("\\f"),
            c if (c as u32) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

/// JSON numbers must be finite; NaN and infinities become `null`.
fn num(v: f64) -> String {
    if !v.is_finite() {
        return "null".to_string();
    }
    if v == v.trunc() && v.abs() < 1e15 {
        format!("{}", v as i64)
    } else {
        format!("{:.6}", v)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strings_are_escaped() {
        assert_eq!(quote("a\"b"), "\"a\\\"b\"");
        assert_eq!(quote("a\\b"), "\"a\\\\b\"");
        assert_eq!(quote("a\nb"), "\"a\\nb\"");
        assert_eq!(quote("\u{1}"), "\"\\u0001\"");
    }

    #[test]
    fn numbers_are_finite() {
        assert_eq!(num(f64::NAN), "null");
        assert_eq!(num(f64::INFINITY), "null");
        assert_eq!(num(3.0), "3");
        assert_eq!(num(2.5), "2.500000");
    }
}
