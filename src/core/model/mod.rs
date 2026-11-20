use rustc_hash::FxHashMap;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Mode {
    Short,
    Long,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Status {
    Pass,
    Warn,
    Fail,
}

impl Status {
    pub fn as_str_lower(self) -> &'static str {
        match self {
            Status::Pass => "pass",
            Status::Warn => "warn",
            Status::Fail => "fail",
        }
    }

    pub fn as_str_upper(self) -> &'static str {
        match self {
            Status::Pass => "PASS",
            Status::Warn => "WARN",
            Status::Fail => "FAIL",
        }
    }

    /// Keep the most severe of two statuses.
    pub fn escalate(&mut self, other: Status) {
        let rank = |s: Status| match s {
            Status::Pass => 0u8,
            Status::Warn => 1,
            Status::Fail => 2,
        };
        if rank(other) > rank(*self) {
            *self = other;
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Encoding {
    Sanger,
    Illumina15,
}

impl Encoding {
    pub fn as_str(self) -> &'static str {
        match self {
            Encoding::Sanger => "Sanger / Illumina 1.9",
            Encoding::Illumina15 => "Illumina 1.5",
        }
    }

    pub fn offset(self) -> u8 {
        match self {
            Encoding::Sanger => 33,
            Encoding::Illumina15 => 64,
        }
    }
}

pub struct FinalizeContext {
    pub encoding: Encoding,
    pub file_name: String,
    pub sample_name: String,
    pub mode: Mode,
    /// Emit one row per base instead of FastQC base groups.
    pub nogroup: bool,
    /// K-mer module is opt-in (FastQC disables it by default too).
    pub kmer_enabled: bool,
    /// Reads fed to the duplication / overrepresentation estimators.
    pub dup_window_reads: u64,
    /// `true` when the duplication window covered the whole file.
    pub dup_window_complete: bool,
    /// Modules suppressed via a FastQC limits file (`ignore` directives).
    pub ignored_modules: Vec<String>,
}

impl FinalizeContext {
    /// `module` is a module anchor such as `sequence_duplication_levels`.
    pub fn is_ignored(&self, module: &str) -> bool {
        if self.ignored_modules.is_empty() {
            return false;
        }
        self.ignored_modules.iter().any(|m| m == module)
    }
}

/// `Per base sequence quality` -> `per_base_sequence_quality`.
pub fn normalize_module_name(name: &str) -> String {
    name.trim()
        .to_ascii_lowercase()
        .chars()
        .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
        .collect::<String>()
        .split('_')
        .filter(|s| !s.is_empty())
        .collect::<Vec<_>>()
        .join("_")
}

pub const MAX_Q: usize = 93;

pub type QualHist = [u64; MAX_Q + 1];

/// `10^(-q/10)` for every representable Phred score.
#[inline]
pub fn error_probability(q: usize) -> f64 {
    error_probability_table()[q.min(MAX_Q)]
}

/// Lookup table for `10^(-q/10)`, built once per process.
#[inline]
pub fn error_probability_table() -> &'static [f64; MAX_Q + 1] {
    static TABLE: std::sync::OnceLock<[f64; MAX_Q + 1]> = std::sync::OnceLock::new();
    TABLE.get_or_init(|| {
        let mut table = [0.0f64; MAX_Q + 1];
        for (q, slot) in table.iter_mut().enumerate() {
            *slot = 10f64.powf(-(q as f64) / 10.0);
        }
        table
    })
}

/// Mean quality of a read expressed as a Phred score derived from the mean
/// error probability. This is the convention used by long-read QC tools
/// (NanoPlot, cramino, chopper); a plain arithmetic mean of Phred scores
/// systematically overestimates read accuracy.
pub fn mean_q_from_error_sum(error_sum: f64, bases: u64) -> f64 {
    if bases == 0 || error_sum <= 0.0 {
        return 0.0;
    }
    let mean_err = error_sum / bases as f64;
    if mean_err <= 0.0 {
        return MAX_Q as f64;
    }
    (-10.0 * mean_err.log10()).clamp(0.0, MAX_Q as f64)
}

pub fn quantile_from_hist(hist: &[u64], q: f64) -> u8 {
    let mut total: u64 = 0;
    for &v in hist {
        total += v;
    }
    if total == 0 {
        return 0;
    }
    let mut rank = (q * total as f64).ceil() as u64;
    if rank < 1 {
        rank = 1;
    }
    let mut cum: u64 = 0;
    for (i, &v) in hist.iter().enumerate() {
        cum += v;
        if cum >= rank {
            return i as u8;
        }
    }
    (hist.len() - 1) as u8
}

/// One row of a per-base module. FastQC collapses positions into widening
/// groups once reads get long, which keeps report size bounded and is what
/// downstream parsers (MultiQC) expect to see.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BaseGroup {
    /// 1-based inclusive.
    pub start: u32,
    /// 1-based inclusive.
    pub end: u32,
}

impl BaseGroup {
    pub fn label(&self) -> String {
        if self.start == self.end {
            self.start.to_string()
        } else {
            format!("{}-{}", self.start, self.end)
        }
    }

    /// Plotting x coordinate.
    pub fn mid(&self) -> f64 {
        (self.start as f64 + self.end as f64) / 2.0
    }

    /// 0-based half-open range into per-position accumulators.
    pub fn range(&self) -> std::ops::Range<usize> {
        (self.start as usize - 1)..(self.end as usize)
    }
}

/// FastQC's default linear grouping, with positions 1-9 shown individually.
pub fn make_base_groups(max_len: usize, nogroup: bool) -> Vec<BaseGroup> {
    if max_len == 0 {
        return Vec::new();
    }
    if nogroup || max_len <= 75 {
        return (1..=max_len)
            .map(|p| BaseGroup {
                start: p as u32,
                end: p as u32,
            })
            .collect();
    }

    let mut interval = 1;
    for magnitude in std::iter::successors(Some(1usize), |scale| scale.checked_mul(10)) {
        if let Some(width) = [2, 5, 10]
            .into_iter()
            .map(|base| base * magnitude)
            .find(|&width| (max_len - 9).div_ceil(width) < 66)
        {
            interval = width;
            break;
        }
    }

    let mut groups = Vec::with_capacity(75);
    let mut start = 1usize;
    while start <= max_len {
        let end = match start {
            1..=9 => start,
            10 if interval > 10 => interval - 1,
            _ => start.saturating_add(interval - 1),
        }
        .min(max_len);
        groups.push(BaseGroup {
            start: start as u32,
            end: end as u32,
        });
        if end == max_len {
            break;
        }
        start = end + 1;
    }
    groups
}

/// Exact read-length histogram with bounded memory: a dense array for the
/// lengths short-read data actually uses and a sparse map for the long tail.
/// Exactness matters because N50/N90 are computed from it.
#[derive(Clone, Debug, Default)]
pub struct LengthHistogram {
    dense: Vec<u64>,
    sparse: FxHashMap<u32, u64>,
}

const DENSE_MAX: usize = 65_536;

impl LengthHistogram {
    pub fn new() -> Self {
        Self::default()
    }

    #[inline]
    pub fn add(&mut self, len: usize) {
        if len < DENSE_MAX {
            if self.dense.len() <= len {
                self.dense.resize(len + 1, 0);
            }
            self.dense[len] += 1;
        } else {
            *self.sparse.entry(len as u32).or_insert(0) += 1;
        }
    }

    pub fn merge(&mut self, other: &LengthHistogram) {
        if self.dense.len() < other.dense.len() {
            self.dense.resize(other.dense.len(), 0);
        }
        for (i, &c) in other.dense.iter().enumerate() {
            self.dense[i] += c;
        }
        for (&len, &c) in &other.sparse {
            *self.sparse.entry(len).or_insert(0) += c;
        }
    }

    pub fn is_empty(&self) -> bool {
        self.dense.iter().all(|&c| c == 0) && self.sparse.is_empty()
    }

    /// Ascending `(length, count)` pairs, skipping empty bins.
    pub fn iter_sorted(&self) -> Vec<(u64, u64)> {
        let mut out: Vec<(u64, u64)> = self
            .dense
            .iter()
            .enumerate()
            .filter(|(_, c)| **c > 0)
            .map(|(len, &c)| (len as u64, c))
            .collect();
        let mut tail: Vec<(u64, u64)> = self
            .sparse
            .iter()
            .map(|(&len, &c)| (len as u64, c))
            .collect();
        tail.sort_unstable_by_key(|&(len, _)| len);
        out.extend(tail);
        out
    }

    pub fn total_reads(&self) -> u64 {
        self.dense.iter().sum::<u64>() + self.sparse.values().sum::<u64>()
    }

    pub fn total_bases(&self) -> u64 {
        self.dense
            .iter()
            .enumerate()
            .map(|(len, &c)| len as u64 * c)
            .sum::<u64>()
            + self
                .sparse
                .iter()
                .map(|(&len, &c)| len as u64 * c)
                .sum::<u64>()
    }

    pub fn min_len(&self) -> u64 {
        self.iter_sorted().first().map(|&(l, _)| l).unwrap_or(0)
    }

    pub fn max_len(&self) -> u64 {
        self.iter_sorted().last().map(|&(l, _)| l).unwrap_or(0)
    }

    pub fn distinct_lengths(&self) -> usize {
        self.iter_sorted().len()
    }

    /// Exact N50 / N90: the length at which the cumulative sum of bases from
    /// the longest read downwards reaches `frac` of all bases.
    pub fn nxx(&self, frac: f64) -> u64 {
        let total_bases = self.total_bases();
        if total_bases == 0 {
            return 0;
        }
        let target = (total_bases as f64 * frac).ceil() as u64;
        let sorted = self.iter_sorted();
        let mut acc: u64 = 0;
        for &(len, count) in sorted.iter().rev() {
            acc += len * count;
            if acc >= target {
                return len;
            }
        }
        0
    }

    /// Read-count quantile (median read length and friends).
    pub fn quantile(&self, frac: f64) -> u64 {
        let total = self.total_reads();
        if total == 0 {
            return 0;
        }
        let target = (total as f64 * frac).ceil().max(1.0) as u64;
        let mut acc: u64 = 0;
        for (len, count) in self.iter_sorted() {
            acc += count;
            if acc >= target {
                return len;
            }
        }
        0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn error_table_matches_powf() {
        for q in [0usize, 1, 10, 20, 30, 40, 60, 93] {
            let expect = 10f64.powf(-(q as f64) / 10.0);
            let got = error_probability(q);
            assert!(
                (got - expect).abs() <= expect * 1e-12,
                "q={q}: {got} vs {expect}"
            );
        }
    }

    #[test]
    fn base_groups_are_ungrouped_for_short_reads() {
        let groups = make_base_groups(75, false);
        assert_eq!(groups.len(), 75);
        assert_eq!(groups[0], BaseGroup { start: 1, end: 1 });
        assert_eq!(groups[74], BaseGroup { start: 75, end: 75 });
    }

    #[test]
    fn base_groups_widen_for_long_reads() {
        let groups = make_base_groups(101, false);
        // first nine positions stay individual
        for (i, g) in groups.iter().take(9).enumerate() {
            assert_eq!(g.start, i as u32 + 1);
            assert_eq!(g.end, i as u32 + 1);
        }
        assert_eq!(groups[9], BaseGroup { start: 10, end: 11 });
        // groups must tile the whole range without gaps or overlap
        let mut expect_next = 1u32;
        for g in &groups {
            assert_eq!(g.start, expect_next);
            expect_next = g.end + 1;
        }
        assert_eq!(expect_next, 102);
    }

    #[test]
    fn nogroup_forces_one_row_per_base() {
        let groups = make_base_groups(300, true);
        assert_eq!(groups.len(), 300);
    }

    #[test]
    fn n50_is_exact() {
        let mut h = LengthHistogram::new();
        // 1 x 100, 2 x 50 => 200 bases, half is 100 -> N50 = 100
        h.add(100);
        h.add(50);
        h.add(50);
        assert_eq!(h.total_bases(), 200);
        assert_eq!(h.nxx(0.5), 100);
        assert_eq!(h.nxx(0.9), 50);
    }

    #[test]
    fn sparse_tail_is_counted() {
        let mut h = LengthHistogram::new();
        h.add(10);
        h.add(1_000_000);
        assert_eq!(h.total_reads(), 2);
        assert_eq!(h.total_bases(), 1_000_010);
        assert_eq!(h.max_len(), 1_000_000);
        assert_eq!(h.min_len(), 10);
    }

    #[test]
    fn merge_is_order_independent() {
        let mut a = LengthHistogram::new();
        a.add(10);
        a.add(200_000);
        let mut b = LengthHistogram::new();
        b.add(10);
        b.add(300_000);

        let mut ab = a.clone();
        ab.merge(&b);
        let mut ba = b.clone();
        ba.merge(&a);
        assert_eq!(ab.iter_sorted(), ba.iter_sorted());
    }
}
