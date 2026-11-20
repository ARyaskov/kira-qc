//! Sequence duplication and overrepresentation.
//!
//! FastQC keeps an exact map of observed sequences (truncated to 50 bp) and
//! derives both the duplication histogram and the overrepresented-sequence
//! table from it. The previous implementation here used a Space-Saving
//! heavy-hitter sketch, which reports every sequence as a heavy hitter once the
//! number of distinct sequences exceeds the sketch capacity: on any real
//! library with more than `DUP_TRACK_CAP` distinct reads it produced a
//! duplication level of 100 % and a spurious FAIL.
//!
//! This module counts exactly, over a bounded prefix of the file (the
//! "duplication window"). Because the window is a prefix defined by global read
//! index, the result does not depend on how chunks are distributed across
//! worker threads.

use rustc_hash::FxHashMap;

/// Reads fed to the estimator by default. Counting is exact within the window.
pub const DEFAULT_DUP_WINDOW: u64 = 1_000_000;

/// Hard ceiling on distinct sequences in a duplication window.
pub const DUP_DISTINCT_CAP: usize = 2_000_000;

/// FastQC truncates sequences longer than this before hashing them.
pub const DUP_TRUNCATE_OVER: usize = 75;
/// ... down to this many bases.
pub const DUP_TRUNCATE_TO: usize = 50;

/// FastQC's duplication levels: exact counts 1-9, then widening buckets.
pub const DUP_LEVEL_LABELS: [&str; 16] = [
    "1", "2", "3", "4", "5", "6", "7", "8", "9", ">10", ">50", ">100", ">500", ">1k", ">5k",
    ">10k+",
];

pub fn dup_level_index(count: u64) -> usize {
    match count {
        0 => 0,
        1..=9 => count as usize - 1,
        10..=49 => 9,
        50..=99 => 10,
        100..=499 => 11,
        500..=999 => 12,
        1_000..=4_999 => 13,
        5_000..=9_999 => 14,
        _ => 15,
    }
}

#[derive(Clone, Debug)]
pub struct DuplicationRow {
    pub level: &'static str,
    /// Share of the deduplicated (distinct) sequences at this level, in percent.
    pub percent_of_deduplicated: f64,
    /// Share of all observed sequences at this level, in percent.
    pub percent_of_total: f64,
}

#[derive(Clone, Debug)]
pub struct OverrepRow {
    pub sequence: String,
    pub count: u64,
    pub percent: f64,
    pub source: String,
}

/// One tracked sequence.
#[derive(Clone, Copy, Debug)]
struct Slot {
    hash: u64,
    count: u64,
    /// Global index of the read that first produced this sequence. Used to make
    /// truncation deterministic when the distinct cap is reached.
    first_index: u64,
    seq_off: u32,
    seq_len: u16,
}

/// Exact sequence counter over a bounded window of reads.
#[derive(Clone, Debug)]
pub struct SeqTracker {
    distinct_cap: usize,
    index: FxHashMap<u64, u32>,
    slots: Vec<Slot>,
    /// Packed storage for the truncated sequences; avoids one allocation per slot.
    arena: Vec<u8>,
    /// Reads offered to the tracker (i.e. inside the window).
    observed_reads: u64,
    /// Reads dropped because the distinct cap was hit.
    dropped_reads: u64,
}

impl Default for SeqTracker {
    fn default() -> Self {
        Self::new()
    }
}

impl SeqTracker {
    pub fn new() -> Self {
        Self {
            distinct_cap: DUP_DISTINCT_CAP,
            index: FxHashMap::default(),
            slots: Vec::new(),
            arena: Vec::new(),
            observed_reads: 0,
            dropped_reads: 0,
        }
    }

    /// Reserve for the share of the window this worker is likely to see.
    pub fn with_capacity(entries: usize) -> Self {
        let entries = entries.min(DUP_DISTINCT_CAP);
        Self {
            distinct_cap: DUP_DISTINCT_CAP,
            index: FxHashMap::with_capacity_and_hasher(entries, Default::default()),
            slots: Vec::with_capacity(entries),
            arena: Vec::with_capacity(entries * DUP_TRUNCATE_TO),
            observed_reads: 0,
            dropped_reads: 0,
        }
    }

    pub fn observed_reads(&self) -> u64 {
        self.observed_reads
    }

    pub fn dropped_reads(&self) -> u64 {
        self.dropped_reads
    }

    pub(crate) fn ensure_complete(&self) -> anyhow::Result<()> {
        anyhow::ensure!(
            self.dropped_reads == 0,
            "duplication window exceeded the {} distinct sequence limit; use a smaller --dup-window",
            self.distinct_cap
        );
        Ok(())
    }

    pub fn distinct(&self) -> usize {
        self.slots.len()
    }

    pub fn is_empty(&self) -> bool {
        self.slots.is_empty()
    }

    /// Record one read. `read_index` is the read's global position in the file.
    pub fn add(&mut self, seq: &[u8], read_index: u64) {
        let truncated = truncate(seq);
        let hash = hash_seq(truncated);
        self.observed_reads += 1;

        if let Some(&slot) = self.index.get(&hash) {
            let s = &mut self.slots[slot as usize];
            s.count += 1;
            if read_index < s.first_index {
                s.first_index = read_index;
            }
            return;
        }

        if self.slots.len() >= self.distinct_cap {
            self.dropped_reads += 1;
            return;
        }

        let seq_off = self.arena.len() as u32;
        self.arena.extend_from_slice(truncated);
        let idx = self.slots.len() as u32;
        self.slots.push(Slot {
            hash,
            count: 1,
            first_index: read_index,
            seq_off,
            seq_len: truncated.len() as u16,
        });
        self.index.insert(hash, idx);
    }

    pub fn merge(&mut self, other: &SeqTracker) {
        self.observed_reads += other.observed_reads;
        self.dropped_reads += other.dropped_reads;

        // Use first-seen order within this partial result. A result that
        // exceeds the cap is rejected before it can be reported.
        let mut order: Vec<u32> = (0..other.slots.len() as u32).collect();
        order.sort_unstable_by_key(|&i| {
            let s = &other.slots[i as usize];
            (s.first_index, s.hash)
        });

        for i in order {
            let s = other.slots[i as usize];
            let seq = &other.arena[s.seq_off as usize..s.seq_off as usize + s.seq_len as usize];
            match self.index.get(&s.hash) {
                Some(&slot) => {
                    let target = &mut self.slots[slot as usize];
                    target.count += s.count;
                    if s.first_index < target.first_index {
                        target.first_index = s.first_index;
                    }
                }
                None => {
                    if self.slots.len() >= self.distinct_cap {
                        self.dropped_reads += s.count;
                        continue;
                    }
                    let seq_off = self.arena.len() as u32;
                    self.arena.extend_from_slice(seq);
                    let idx = self.slots.len() as u32;
                    self.slots.push(Slot {
                        hash: s.hash,
                        count: s.count,
                        first_index: s.first_index,
                        seq_off,
                        seq_len: s.seq_len,
                    });
                    self.index.insert(s.hash, idx);
                }
            }
        }
    }

    /// `(sequence, count, first_index)` in deterministic order.
    pub fn entries_sorted(&self) -> Vec<(&[u8], u64, u64)> {
        let mut out: Vec<(&[u8], u64, u64)> = self
            .slots
            .iter()
            .map(|s| {
                let seq = &self.arena[s.seq_off as usize..s.seq_off as usize + s.seq_len as usize];
                (seq, s.count, s.first_index)
            })
            .collect();
        out.sort_unstable_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(b.0)));
        out
    }

    /// Duplication histogram over the tracked window.
    pub fn duplication_levels(&self) -> DuplicationSummary {
        let mut distinct_at_level = [0u64; 16];
        let mut reads_at_level = [0u64; 16];
        let mut total_distinct: u64 = 0;
        let mut total_reads: u64 = 0;

        for s in &self.slots {
            let idx = dup_level_index(s.count);
            distinct_at_level[idx] += 1;
            reads_at_level[idx] += s.count;
            total_distinct += 1;
            total_reads += s.count;
        }

        let mut rows = Vec::with_capacity(16);
        for i in 0..16 {
            let pct_dedup = if total_distinct == 0 {
                0.0
            } else {
                distinct_at_level[i] as f64 * 100.0 / total_distinct as f64
            };
            let pct_total = if total_reads == 0 {
                0.0
            } else {
                reads_at_level[i] as f64 * 100.0 / total_reads as f64
            };
            rows.push(DuplicationRow {
                level: DUP_LEVEL_LABELS[i],
                percent_of_deduplicated: pct_dedup,
                percent_of_total: pct_total,
            });
        }

        let total_deduplicated_percent = if total_reads == 0 {
            100.0
        } else {
            total_distinct as f64 * 100.0 / total_reads as f64
        };

        DuplicationSummary {
            rows,
            total_deduplicated_percent,
            observed_reads: total_reads,
            distinct_reads: total_distinct,
        }
    }
}

#[derive(Clone, Debug)]
pub struct DuplicationSummary {
    pub rows: Vec<DuplicationRow>,
    /// FastQC's headline number: what fraction of the library would remain
    /// after deduplication.
    pub total_deduplicated_percent: f64,
    pub observed_reads: u64,
    pub distinct_reads: u64,
}

/// FastQC status thresholds: warn when more than 20 % of the library is
/// duplicated, fail beyond 50 %.
pub fn duplication_status(total_deduplicated_percent: f64) -> crate::core::model::Status {
    use crate::core::model::Status;
    let duplicated = 100.0 - total_deduplicated_percent;
    if duplicated > 50.0 {
        Status::Fail
    } else if duplicated > 20.0 {
        Status::Warn
    } else {
        Status::Pass
    }
}

#[inline]
pub fn truncate(seq: &[u8]) -> &[u8] {
    if seq.len() > DUP_TRUNCATE_OVER {
        &seq[..DUP_TRUNCATE_TO]
    } else {
        seq
    }
}

/// Case-folded hash. Folding with `& 0xDF` maps lowercase ACGT onto uppercase.
#[inline]
pub fn hash_seq(seq: &[u8]) -> u64 {
    let mut buf: [u8; 128] = [0; 128];
    if seq.len() <= buf.len() {
        for (i, &b) in seq.iter().enumerate() {
            buf[i] = b & 0xDF;
        }
        return xxhash_rust::xxh3::xxh3_64(&buf[..seq.len()]);
    }
    let folded: Vec<u8> = seq.iter().map(|&b| b & 0xDF).collect();
    xxhash_rust::xxh3::xxh3_64(&folded)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_distinct_limit_overflow_invalidates_direct_and_merged_counts() {
        let tracker = || SeqTracker {
            distinct_cap: 1,
            ..SeqTracker::new()
        };
        let mut first = tracker();
        first.add(b"AAAA", 0);
        first.add(b"AAAA", 1);
        assert!(first.ensure_complete().is_ok());
        let mut second = tracker();
        second.add(b"CCCC", 2);
        let mut merged = first.clone();
        merged.merge(&second);
        second.merge(&first);
        first.add(b"CCCC", 2);
        for result in [first, second, merged] {
            assert!(
                result
                    .ensure_complete()
                    .unwrap_err()
                    .to_string()
                    .contains("--dup-window")
            );
        }
    }

    #[test]
    fn counts_are_exact() {
        let mut t = SeqTracker::new();
        t.add(b"ACGTACGTAC", 0);
        t.add(b"ACGTACGTAC", 1);
        t.add(b"TTTTTTTTTT", 2);
        assert_eq!(t.distinct(), 2);
        assert_eq!(t.observed_reads(), 3);

        let summary = t.duplication_levels();
        // one sequence seen twice, one seen once
        assert_eq!(summary.distinct_reads, 2);
        assert_eq!(summary.observed_reads, 3);
        assert!((summary.total_deduplicated_percent - 66.666).abs() < 0.01);
    }

    #[test]
    fn high_complexity_library_is_not_flagged() {
        // The old sketch reported 100 % duplication here; every read is unique.
        let mut t = SeqTracker::new();
        for i in 0..50_000u64 {
            let seq = format!("{:A>50}", i);
            t.add(seq.as_bytes(), i);
        }
        let summary = t.duplication_levels();
        assert_eq!(summary.distinct_reads, summary.observed_reads);
        assert!((summary.total_deduplicated_percent - 100.0).abs() < 1e-9);
        assert_eq!(
            duplication_status(summary.total_deduplicated_percent),
            crate::core::model::Status::Pass
        );
    }

    #[test]
    fn merge_is_order_independent() {
        let mut a = SeqTracker::new();
        a.add(b"AAAA", 0);
        a.add(b"CCCC", 1);
        let mut b = SeqTracker::new();
        b.add(b"CCCC", 2);
        b.add(b"GGGG", 3);

        let mut ab = a.clone();
        ab.merge(&b);
        let mut ba = b.clone();
        ba.merge(&a);

        assert_eq!(ab.entries_sorted(), ba.entries_sorted());
        assert_eq!(ab.observed_reads(), 4);
    }

    #[test]
    fn long_sequences_are_truncated_like_fastqc() {
        let long = vec![b'A'; 200];
        assert_eq!(truncate(&long).len(), DUP_TRUNCATE_TO);
        let short = vec![b'A'; 60];
        assert_eq!(truncate(&short).len(), 60);
    }

    #[test]
    fn case_folding_matches() {
        assert_eq!(hash_seq(b"acgt"), hash_seq(b"ACGT"));
        assert_ne!(hash_seq(b"ACGT"), hash_seq(b"ACGA"));
    }

    #[test]
    fn duplication_levels_bucket_correctly() {
        assert_eq!(dup_level_index(1), 0);
        assert_eq!(dup_level_index(9), 8);
        assert_eq!(dup_level_index(10), 9);
        assert_eq!(dup_level_index(49), 9);
        assert_eq!(dup_level_index(50), 10);
        assert_eq!(dup_level_index(20_000), 15);
    }
}
