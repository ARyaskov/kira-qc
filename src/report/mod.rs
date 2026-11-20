pub mod fastqc_txt;
pub mod html;
pub mod json;
#[cfg(feature = "latex")]
pub mod latex;
pub mod summary_txt;
pub mod zip;

use crate::core::metrics::FinalMetrics;
use crate::core::model::{FinalizeContext, Mode, Status};

/// Modules the report writers know about, in FastQC's order.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ModuleId {
    Basic,
    PerBaseQual,
    PerTile,
    PerSeqQual,
    PerBaseContent,
    PerSeqGc,
    PerBaseN,
    PerSeqN,
    LengthDist,
    Duplication,
    Overrepresented,
    AdapterContent,
    KmerContent,
}

impl ModuleId {
    /// FastQC's module title, used verbatim in `fastqc_data.txt` and `summary.txt`.
    pub fn title(self) -> &'static str {
        match self {
            ModuleId::Basic => "Basic Statistics",
            ModuleId::PerBaseQual => "Per base sequence quality",
            ModuleId::PerTile => "Per tile sequence quality",
            ModuleId::PerSeqQual => "Per sequence quality scores",
            ModuleId::PerBaseContent => "Per base sequence content",
            ModuleId::PerSeqGc => "Per sequence GC content",
            ModuleId::PerBaseN => "Per base N content",
            ModuleId::PerSeqN => "Per sequence N content",
            ModuleId::LengthDist => "Sequence Length Distribution",
            ModuleId::Duplication => "Sequence Duplication Levels",
            ModuleId::Overrepresented => "Overrepresented sequences",
            ModuleId::AdapterContent => "Adapter Content",
            ModuleId::KmerContent => "Kmer Content",
        }
    }

    /// Anchor used by the HTML report.
    pub fn anchor(self) -> &'static str {
        match self {
            ModuleId::Basic => "basic_statistics",
            ModuleId::PerBaseQual => "per_base_sequence_quality",
            ModuleId::PerTile => "per_tile_sequence_quality",
            ModuleId::PerSeqQual => "per_sequence_quality_scores",
            ModuleId::PerBaseContent => "per_base_sequence_content",
            ModuleId::PerSeqGc => "per_sequence_gc_content",
            ModuleId::PerBaseN => "per_base_n_content",
            ModuleId::PerSeqN => "per_sequence_n_content",
            ModuleId::LengthDist => "sequence_length_distribution",
            ModuleId::Duplication => "sequence_duplication_levels",
            ModuleId::Overrepresented => "overrepresented_sequences",
            ModuleId::AdapterContent => "adapter_content",
            ModuleId::KmerContent => "kmer_content",
        }
    }

    pub fn status(self, metrics: &FinalMetrics) -> Status {
        let s = &metrics.statuses;
        match self {
            ModuleId::Basic => s.basic,
            ModuleId::PerBaseQual => s.per_base_qual,
            ModuleId::PerTile => s.per_tile,
            ModuleId::PerSeqQual => s.per_seq_qual,
            ModuleId::PerBaseContent => s.per_base_content,
            ModuleId::PerSeqGc => s.per_seq_gc,
            ModuleId::PerBaseN => s.per_base_n,
            ModuleId::PerSeqN => s.per_seq_n,
            ModuleId::LengthDist => s.length_dist,
            ModuleId::Duplication => s.duplication,
            ModuleId::Overrepresented => s.overrepresented,
            ModuleId::AdapterContent => s.adapter_content,
            ModuleId::KmerContent => s.kmer_content,
        }
    }
}

/// The modules this run actually reports, honouring mode, per-tile
/// availability, the K-mer switch and any `ignore` directives from a limits
/// file. Every writer walks this list so the text, summary, HTML and JSON
/// outputs can never disagree about which modules exist.
pub fn active_modules(metrics: &FinalMetrics, ctx: &FinalizeContext) -> Vec<ModuleId> {
    let candidates: &[ModuleId] = match ctx.mode {
        Mode::Short => &[
            ModuleId::Basic,
            ModuleId::PerBaseQual,
            ModuleId::PerTile,
            ModuleId::PerSeqQual,
            ModuleId::PerBaseContent,
            ModuleId::PerSeqGc,
            ModuleId::PerBaseN,
            ModuleId::LengthDist,
            ModuleId::Duplication,
            ModuleId::Overrepresented,
            ModuleId::AdapterContent,
            ModuleId::KmerContent,
        ],
        Mode::Long => &[
            ModuleId::Basic,
            ModuleId::LengthDist,
            ModuleId::PerSeqQual,
            ModuleId::PerSeqGc,
            ModuleId::PerSeqN,
            ModuleId::AdapterContent,
        ],
    };

    candidates
        .iter()
        .copied()
        .filter(|m| match m {
            // Only Illumina headers carry tile numbers.
            ModuleId::PerTile => metrics.per_tile_available,
            ModuleId::KmerContent => ctx.kmer_enabled,
            _ => true,
        })
        .filter(|m| !ctx.is_ignored(m.anchor()))
        .collect()
}

/// Every module anchor, for validating limits files.
pub const ALL_MODULE_ANCHORS: [&str; 13] = [
    "basic_statistics",
    "per_base_sequence_quality",
    "per_tile_sequence_quality",
    "per_sequence_quality_scores",
    "per_base_sequence_content",
    "per_sequence_gc_content",
    "per_base_n_content",
    "per_sequence_n_content",
    "sequence_length_distribution",
    "sequence_duplication_levels",
    "overrepresented_sequences",
    "adapter_content",
    "kmer_content",
];

/// Translate a FastQC limits-file module key into a module anchor.
///
/// FastQC's `limits.txt` uses short keys (`duplication`, `kmer`, `tile`), not
/// module titles, so matching on the title alone silently ignored every
/// directive. Module titles and anchors are accepted too.
pub fn resolve_limits_key(key: &str) -> Option<&'static str> {
    let normalized = crate::core::model::normalize_module_name(key);
    let anchor = match normalized.as_str() {
        "quality_base" | "quality_base_lower" | "quality_base_median" => {
            "per_base_sequence_quality"
        }
        "tile" => "per_tile_sequence_quality",
        "quality_sequence" => "per_sequence_quality_scores",
        "sequence" => "per_base_sequence_content",
        "gc_sequence" => "per_sequence_gc_content",
        "n_content" => "per_base_n_content",
        "sequence_length" => "sequence_length_distribution",
        "duplication" => "sequence_duplication_levels",
        "overrepresented" => "overrepresented_sequences",
        "adapter" => "adapter_content",
        "kmer" => "kmer_content",
        "basic" => "basic_statistics",
        other => {
            // Also accept the module title or anchor verbatim.
            return ALL_MODULE_ANCHORS.iter().copied().find(|a| *a == other);
        }
    };
    Some(anchor)
}

/// Report timestamp, overridable through `SOURCE_DATE_EPOCH` so that a run can
/// be reproduced byte for byte.
pub fn report_timestamp() -> u64 {
    if let Ok(v) = std::env::var("SOURCE_DATE_EPOCH")
        && let Ok(secs) = v.trim().parse::<u64>()
    {
        return secs;
    }
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0)
}

/// Civil date/time from a Unix timestamp (Howard Hinnant's algorithm).
pub fn format_timestamp(ts: u64) -> String {
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

/// Thousands separators for report tables.
pub fn fmt_int(v: u64) -> String {
    let s = v.to_string();
    let mut out = String::with_capacity(s.len() + s.len() / 3);
    let len = s.len();
    for (i, ch) in s.chars().enumerate() {
        if i != 0 && (len - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(ch);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn timestamp_formatting_is_stable() {
        assert_eq!(format_timestamp(0), "1970-01-01 00:00:00 UTC");
        assert_eq!(format_timestamp(1_700_000_000), "2023-11-14 22:13:20 UTC");
    }

    #[test]
    fn integers_get_separators() {
        assert_eq!(fmt_int(0), "0");
        assert_eq!(fmt_int(999), "999");
        assert_eq!(fmt_int(1_000), "1,000");
        assert_eq!(fmt_int(1_234_567), "1,234,567");
    }
}
