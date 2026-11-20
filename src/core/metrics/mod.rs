use crate::core::fastq::ReadView;
use crate::core::model::{
    BaseGroup, Encoding, FinalizeContext, LengthHistogram, MAX_Q, Mode, QualHist, Status,
    error_probability_table, make_base_groups, mean_q_from_error_sum, quantile_from_hist,
};
use crate::simd;
use std::sync::Arc;
use std::time::Duration;

pub mod adapter_content;
pub mod basic;
pub mod duplication;
pub mod kmer_content;
pub mod length_dist;
pub mod overrepresented;
pub mod per_base_content;
pub mod per_base_n;
pub mod per_base_qual;
pub mod per_seq_gc;
pub mod per_seq_n;
pub mod per_seq_qual;
pub mod per_tile_qual;

pub use adapter_content::{Adapter, AdapterRow, AdapterSet, MAX_ADAPTERS};
pub use basic::BasicStats;
pub use duplication::{
    DEFAULT_DUP_WINDOW, DuplicationRow, DuplicationSummary, OverrepRow, SeqTracker,
};
pub use kmer_content::KmerRow;
pub use length_dist::{LengthDistRow, LongLengthSummary};
pub use overrepresented::ContaminantList;
pub use per_base_content::PerBaseContentRow;
pub use per_base_n::PerBaseNRow;
pub use per_base_qual::PerBaseQualRow;
pub use per_seq_gc::PerSeqGcRow;
pub use per_seq_n::PerSeqNRow;
pub use per_seq_qual::PerSeqQualRow;
pub use per_tile_qual::{PerTileQual, PerTileRow};

#[derive(Clone, Debug, Default)]
pub struct BaseCounts {
    pub a: u64,
    pub c: u64,
    pub g: u64,
    pub t: u64,
    pub n: u64,
}

impl BaseCounts {
    fn zero() -> Self {
        Self::default()
    }

    fn add_assign(&mut self, other: &BaseCounts) {
        self.a += other.a;
        self.c += other.c;
        self.g += other.g;
        self.t += other.t;
        self.n += other.n;
    }

    fn called(&self) -> u64 {
        self.a + self.c + self.g + self.t
    }

    fn total(&self) -> u64 {
        self.called() + self.n
    }
}

/// Optional per-worker timing breakdown, gated on `KIRA_STATS=1`.
#[derive(Clone, Debug, Default)]
pub struct UpdateTimings {
    pub kmer: Duration,
    pub kmer_updates: u64,
}

#[derive(Clone, Debug)]
pub struct AggConfig {
    pub mode: Mode,
    /// K-mer content is opt-in.
    pub kmer_enabled: bool,
    /// Reads fed to the duplication estimator; 0 means "every read".
    pub dup_window: u64,
    /// Collect per-tile quality when read headers carry tile numbers.
    pub per_tile_enabled: bool,
    pub adapters: Arc<AdapterSet>,
    pub contaminants: Arc<ContaminantList>,
}

impl AggConfig {
    pub fn defaults_for(mode: Mode) -> Self {
        Self {
            mode,
            kmer_enabled: false,
            dup_window: DEFAULT_DUP_WINDOW,
            per_tile_enabled: mode == Mode::Short,
            adapters: AdapterSet::default_for(mode),
            contaminants: ContaminantList::builtin(),
        }
    }
}

impl Default for AggConfig {
    fn default() -> Self {
        Self::defaults_for(Mode::Short)
    }
}

#[derive(Clone, Debug)]
pub struct Agg {
    pub cfg: AggConfig,
    pub total_reads: u64,
    pub total_bases: u64,
    pub gc_bases: u64,
    pub n_bases: u64,
    pub empty_reads: u64,
    pub length_hist: LengthHistogram,
    pub per_pos_qual: Vec<QualHist>,
    pub per_pos_base: Vec<BaseCounts>,
    pub per_seq_mean_q_hist: Vec<u64>,
    pub per_seq_gc_hist: [u64; 101],
    pub per_seq_n_hist: [u64; 101],
    pub reads_n_gt10: u64,
    pub reads_n_gt20: u64,
    /// Flat `positions * adapter_stride` matrix.
    pub adapter_counts: Vec<u64>,
    pub adapter_stride: usize,
    pub adapter_reads_any: Vec<u64>,
    pub seq_tracker: SeqTracker,
    pub per_tile: PerTileQual,
    pub kmer_bins: Vec<kmer_content::BinCounts>,
    pub kmer_bin_counts: [u64; kmer_content::BINS],
    pub kmer_total: u64,
}

/// Base slot lookup: 0=A 1=C 2=G 3=T 4=N 5=other.
const fn build_base_lut() -> [u8; 256] {
    let mut lut = [5u8; 256];
    lut[b'A' as usize] = 0;
    lut[b'C' as usize] = 1;
    lut[b'G' as usize] = 2;
    lut[b'T' as usize] = 3;
    lut[b'N' as usize] = 4;
    lut
}
static BASE_LUT: [u8; 256] = build_base_lut();

impl Agg {
    pub fn new(cfg: AggConfig) -> Self {
        let kmer_bins = if cfg.mode == Mode::Short && cfg.kmer_enabled {
            (0..kmer_content::BINS)
                .map(|_| kmer_content::BinCounts::new())
                .collect()
        } else {
            Vec::new()
        };
        let stride = cfg.adapters.len();
        let adapter_reads_any = vec![0u64; stride];

        Self {
            cfg,
            total_reads: 0,
            total_bases: 0,
            gc_bases: 0,
            n_bases: 0,
            empty_reads: 0,
            length_hist: LengthHistogram::new(),
            per_pos_qual: Vec::new(),
            per_pos_base: Vec::new(),
            per_seq_mean_q_hist: vec![0u64; MAX_Q + 1],
            per_seq_gc_hist: [0u64; 101],
            per_seq_n_hist: [0u64; 101],
            reads_n_gt10: 0,
            reads_n_gt20: 0,
            adapter_counts: Vec::new(),
            adapter_stride: stride,
            adapter_reads_any,
            seq_tracker: SeqTracker::new(),
            per_tile: PerTileQual::new(),
            kmer_bins,
            kmer_bin_counts: [0u64; kmer_content::BINS],
            kmer_total: 0,
        }
    }

    pub fn mode(&self) -> Mode {
        self.cfg.mode
    }

    pub fn update_read(&mut self, read: &ReadView<'_>, phred_offset: u8, read_index: u64) {
        self.update_read_with_timings(read, phred_offset, read_index, None);
    }

    pub fn update_read_timed(
        &mut self,
        read: &ReadView<'_>,
        phred_offset: u8,
        read_index: u64,
        timing: &mut UpdateTimings,
    ) {
        self.update_read_with_timings(read, phred_offset, read_index, Some(timing));
    }

    fn update_read_with_timings(
        &mut self,
        read: &ReadView<'_>,
        phred_offset: u8,
        read_index: u64,
        timing: Option<&mut UpdateTimings>,
    ) {
        let len = read.seq.len();

        // Zero-length records still count towards the library, and a length of
        // zero is what makes the length-distribution module fail.
        self.total_reads += 1;
        self.length_hist.add(len);
        if len == 0 {
            self.empty_reads += 1;
            return;
        }
        self.total_bases += len as u64;

        match self.cfg.mode {
            Mode::Short => self.update_short(read, len, phred_offset, read_index, timing),
            Mode::Long => self.update_long(read, len, phred_offset),
        }
    }

    /// Short-read path: one fused pass over sequence and quality feeds the
    /// per-position accumulators, the base composition and the mean quality.
    fn update_short(
        &mut self,
        read: &ReadView<'_>,
        len: usize,
        phred_offset: u8,
        read_index: u64,
        timing: Option<&mut UpdateTimings>,
    ) {
        if self.per_pos_qual.len() < len {
            self.per_pos_qual.resize(len, [0u64; MAX_Q + 1]);
        }
        if self.per_pos_base.len() < len {
            self.per_pos_base.resize(len, BaseCounts::zero());
        }

        let seq = &read.seq[..len];
        let qual = &read.qual[..len];
        let per_pos_base = &mut self.per_pos_base[..len];
        let per_pos_qual = &mut self.per_pos_qual[..len];

        let mut gc: u64 = 0;
        let mut n_count: u64 = 0;
        let mut sum_q: u64 = 0;

        for i in 0..len {
            let bc = &mut per_pos_base[i];
            match BASE_LUT[(seq[i] & 0xDF) as usize] {
                0 => bc.a += 1,
                1 => {
                    bc.c += 1;
                    gc += 1;
                }
                2 => {
                    bc.g += 1;
                    gc += 1;
                }
                3 => bc.t += 1,
                4 => {
                    bc.n += 1;
                    n_count += 1;
                }
                _ => {}
            }

            let q = (qual[i].saturating_sub(phred_offset) as usize).min(MAX_Q);
            per_pos_qual[i][q] += 1;
            sum_q += q as u64;
        }

        self.gc_bases += gc;
        self.n_bases += n_count;

        // FastQC averages Phred scores arithmetically for short reads.
        let mean_q = sum_q / len as u64;
        self.per_seq_mean_q_hist[(mean_q as usize).min(MAX_Q)] += 1;

        let gc_percent = ((gc * 100) + (len as u64 / 2)) / len as u64;
        self.per_seq_gc_hist[(gc_percent as usize).min(100)] += 1;

        let stride = self.adapter_stride;
        if self.adapter_counts.len() < len * stride {
            self.adapter_counts.resize(len * stride, 0);
        }
        self.cfg
            .adapters
            .scan(read.seq, &mut self.adapter_counts, stride);

        if self.in_dup_window(read_index) {
            self.seq_tracker.add(read.seq, read_index);
        }

        if self.cfg.per_tile_enabled
            && let Some(tile) = read.tile
        {
            self.per_tile.add(tile, qual, phred_offset);
        }

        if self.cfg.kmer_enabled && len >= kmer_content::K {
            kmer_content::update_kmers(
                read.seq,
                len,
                &mut self.kmer_bins,
                &mut self.kmer_bin_counts,
                &mut self.kmer_total,
                timing,
            );
        } else {
            let _ = timing;
        }
    }

    /// Long-read path: no per-position modules, and read quality is the Phred
    /// score of the mean error probability rather than the mean Phred score.
    fn update_long(&mut self, read: &ReadView<'_>, len: usize, phred_offset: u8) {
        let (_a, c, g, _t, n) = simd::count_bases(read.seq);
        let gc = c as u64 + g as u64;
        let n_count = n as u64;
        self.gc_bases += gc;
        self.n_bases += n_count;

        let mut qhist = [0u32; MAX_Q + 1];
        for &q in &read.qual[..len] {
            qhist[(q.saturating_sub(phred_offset) as usize).min(MAX_Q)] += 1;
        }
        let table = error_probability_table();
        let mut error_sum = 0.0f64;
        for (q, &count) in qhist.iter().enumerate() {
            if count > 0 {
                error_sum += table[q] * count as f64;
            }
        }
        let mean_q = mean_q_from_error_sum(error_sum, len as u64);
        self.per_seq_mean_q_hist[(mean_q.round() as usize).min(MAX_Q)] += 1;

        let gc_percent = ((gc * 100) + (len as u64 / 2)) / len as u64;
        self.per_seq_gc_hist[(gc_percent as usize).min(100)] += 1;

        let n_percent = ((n_count * 100) + (len as u64 / 2)) / len as u64;
        self.per_seq_n_hist[n_percent.min(100) as usize] += 1;
        if n_count * 100 > len as u64 * 20 {
            self.reads_n_gt20 += 1;
        }
        if n_count * 100 > len as u64 * 10 {
            self.reads_n_gt10 += 1;
        }

        let mut hits = [false; MAX_ADAPTERS];
        let n_adapters = self.adapter_stride;
        self.cfg
            .adapters
            .scan_any(read.seq, &mut hits[..n_adapters]);
        for (i, hit) in hits[..n_adapters].iter().enumerate() {
            if *hit {
                self.adapter_reads_any[i] += 1;
            }
        }
    }

    #[inline]
    fn in_dup_window(&self, read_index: u64) -> bool {
        self.cfg.dup_window == 0 || read_index < self.cfg.dup_window
    }

    pub fn merge(&mut self, other: &Agg) {
        self.total_reads += other.total_reads;
        self.total_bases += other.total_bases;
        self.gc_bases += other.gc_bases;
        self.n_bases += other.n_bases;
        self.empty_reads += other.empty_reads;
        self.length_hist.merge(&other.length_hist);

        if self.per_pos_qual.len() < other.per_pos_qual.len() {
            self.per_pos_qual
                .resize(other.per_pos_qual.len(), [0u64; MAX_Q + 1]);
        }
        if self.per_pos_base.len() < other.per_pos_base.len() {
            self.per_pos_base
                .resize(other.per_pos_base.len(), BaseCounts::zero());
        }
        for (i, hist) in other.per_pos_qual.iter().enumerate() {
            let target = &mut self.per_pos_qual[i];
            for q in 0..=MAX_Q {
                target[q] += hist[q];
            }
        }
        for (i, bc) in other.per_pos_base.iter().enumerate() {
            self.per_pos_base[i].add_assign(bc);
        }
        for i in 0..self.per_seq_mean_q_hist.len() {
            self.per_seq_mean_q_hist[i] += other.per_seq_mean_q_hist[i];
        }
        for i in 0..self.per_seq_gc_hist.len() {
            self.per_seq_gc_hist[i] += other.per_seq_gc_hist[i];
        }
        for i in 0..self.per_seq_n_hist.len() {
            self.per_seq_n_hist[i] += other.per_seq_n_hist[i];
        }
        self.reads_n_gt10 += other.reads_n_gt10;
        self.reads_n_gt20 += other.reads_n_gt20;

        match self.cfg.mode {
            Mode::Short => {
                self.seq_tracker.merge(&other.seq_tracker);
                self.per_tile.merge(&other.per_tile);
                if self.adapter_counts.len() < other.adapter_counts.len() {
                    self.adapter_counts.resize(other.adapter_counts.len(), 0);
                }
                for (t, &s) in self.adapter_counts.iter_mut().zip(&other.adapter_counts) {
                    *t += s;
                }
                if self.cfg.kmer_enabled && !other.kmer_bins.is_empty() {
                    for b in 0..kmer_content::BINS {
                        self.kmer_bins[b].merge(&other.kmer_bins[b]);
                        self.kmer_bin_counts[b] += other.kmer_bin_counts[b];
                    }
                    self.kmer_total += other.kmer_total;
                }
            }
            Mode::Long => {
                for (t, &s) in self
                    .adapter_reads_any
                    .iter_mut()
                    .zip(&other.adapter_reads_any)
                {
                    *t += s;
                }
            }
        }
    }

    /// Number of reads that reach each 0-based position.
    fn reads_reaching(&self, positions: usize) -> Vec<u64> {
        let mut lengths = self.length_hist.iter_sorted().into_iter().rev().peekable();
        let mut reaching = vec![0; positions];
        let mut count = 0;
        for (position, value) in reaching.iter_mut().enumerate().rev() {
            while let Some(&(length, reads)) = lengths.peek() {
                if length <= position as u64 {
                    break;
                }
                count += reads;
                lengths.next();
            }
            *value = count;
        }
        reaching
    }

    pub fn finalize(&self, ctx: &FinalizeContext) -> FinalMetrics {
        let min_len = if self.total_reads == 0 {
            0
        } else {
            self.length_hist.min_len()
        };
        let max_len = if self.total_reads == 0 {
            0
        } else {
            self.length_hist.max_len()
        };
        let gc_percent = if ctx.mode == Mode::Short {
            let called = self.per_pos_base.iter().map(BaseCounts::called).sum();
            (self.gc_bases * 100).checked_div(called).unwrap_or(0)
        } else {
            ((self.gc_bases * 100) + self.total_bases / 2)
                .checked_div(self.total_bases)
                .unwrap_or(0)
        } as u32;

        let basic = BasicStats {
            file_type: "Conventional base calls",
            encoding: match ctx.encoding {
                Encoding::Sanger => Encoding::Sanger.as_str(),
                Encoding::Illumina15 => Encoding::Illumina15.as_str(),
            },
            total_sequences: self.total_reads,
            total_bases: self.total_bases,
            poor_quality_sequences: 0,
            min_len,
            max_len,
            gc_percent,
        };

        let groups = if ctx.mode == Mode::Short {
            make_base_groups(self.per_pos_base.len(), ctx.nogroup)
        } else {
            Vec::new()
        };

        let (per_base_qual, per_base_qual_status) = self.build_per_base_quality(&groups);
        let (per_base_content, per_base_content_status) = self.build_per_base_content(&groups);
        let (per_base_n, per_base_n_status) = self.build_per_base_n(&groups);
        let (per_tile, per_tile_status) = self.per_tile.rows(&groups);

        let per_seq_qual: Vec<PerSeqQualRow> = self
            .per_seq_mean_q_hist
            .iter()
            .enumerate()
            .filter(|(_, c)| **c > 0)
            .map(|(q, &count)| PerSeqQualRow {
                mean_q: q as u8,
                count,
            })
            .collect();
        let per_seq_qual_status = self.per_seq_quality_status(ctx.mode);

        let per_seq_gc: Vec<PerSeqGcRow> = self
            .per_seq_gc_hist
            .iter()
            .enumerate()
            .filter(|(_, c)| **c > 0)
            .map(|(gc, &count)| PerSeqGcRow {
                gc: gc as u8,
                count,
            })
            .collect();
        let per_seq_gc_status = self.per_seq_gc_status();

        let mut per_seq_n = Vec::new();
        let mut per_seq_n_status = Status::Pass;
        if ctx.mode == Mode::Long {
            for (i, &count) in self.per_seq_n_hist.iter().enumerate() {
                if count > 0 {
                    per_seq_n.push(PerSeqNRow {
                        n_percent: i as u8,
                        count,
                    });
                }
            }
            if self.total_reads > 0 {
                let gt20 = self.reads_n_gt20 as f64 / self.total_reads as f64 * 100.0;
                let gt10 = self.reads_n_gt10 as f64 / self.total_reads as f64 * 100.0;
                if gt20 > 5.0 {
                    per_seq_n_status = Status::Fail;
                } else if gt10 > 5.0 {
                    per_seq_n_status = Status::Warn;
                }
            }
        }

        let mut length_dist = Vec::new();
        let mut long_length = None;
        if ctx.mode == Mode::Short {
            for (length, count) in self.length_hist.iter_sorted() {
                if count > 0 {
                    length_dist.push(LengthDistRow { length, count });
                }
            }
        } else {
            long_length = Some(self.build_long_length());
        }
        let length_dist_status = self.length_status();

        let (duplication, overrepresented, duplication_status, overrep_status, dup_summary) =
            self.build_duplication(ctx);

        let (adapter_content, adapter_status) = self.build_adapter_content(ctx, &groups);

        let (kmer_rows, kmer_status) = if ctx.mode == Mode::Short && ctx.kmer_enabled {
            kmer_content::build_rows(
                &self.kmer_bins,
                &self.kmer_bin_counts,
                self.kmer_total,
                max_len,
            )
        } else {
            (Vec::new(), Status::Pass)
        };

        let statuses = Statuses {
            basic: Status::Pass,
            per_base_qual: per_base_qual_status,
            per_tile: per_tile_status,
            per_seq_qual: per_seq_qual_status,
            per_base_content: per_base_content_status,
            per_seq_gc: per_seq_gc_status,
            per_base_n: per_base_n_status,
            length_dist: length_dist_status,
            duplication: duplication_status,
            overrepresented: overrep_status,
            adapter_content: adapter_status,
            per_seq_n: per_seq_n_status,
            kmer_content: kmer_status,
        };

        FinalMetrics {
            basic,
            groups,
            per_base_qual,
            per_tile,
            per_tile_available: self.per_tile.is_available(),
            per_seq_qual,
            per_base_content,
            per_seq_gc,
            per_base_n,
            length_dist,
            duplication,
            duplication_summary: dup_summary,
            overrepresented,
            adapter_content,
            adapter_names: self.cfg.adapters.names(),
            per_seq_n,
            long_length,
            kmer_rows,
            statuses,
        }
    }

    fn build_per_base_quality(&self, groups: &[BaseGroup]) -> (Vec<PerBaseQualRow>, Status) {
        let mut rows = Vec::with_capacity(groups.len());
        let mut status = Status::Pass;
        for g in groups {
            let mut mean_sum = 0.0;
            let mut positions = 0;
            let mut percentiles = [0.0; 5];
            let mut covered = 0;
            for i in g.range() {
                let Some(h) = self.per_pos_qual.get(i) else {
                    continue;
                };
                let (total, sum) = h.iter().enumerate().fold((0, 0), |(total, sum), (q, &n)| {
                    (total + n, sum + n * q as u64)
                });
                if total != 0 {
                    mean_sum += sum as f64 / total as f64;
                    positions += 1;
                }
                if total > 100 {
                    for (sum, value) in percentiles
                        .iter_mut()
                        .zip(per_base_qual::position_percentiles(h, total))
                    {
                        *sum += value;
                    }
                    covered += 1;
                }
            }
            if positions == 0 {
                continue;
            }
            let mean = mean_sum / positions as f64;
            let [p10, lq, median, uq, p90] = percentiles.map(|v| v / covered as f64);

            // FastQC: warn on lower quartile < 10 or median < 25,
            // fail on lower quartile < 5 or median < 20.
            if lq < 5.0 || median < 20.0 {
                status = Status::Fail;
            } else if (lq < 10.0 || median < 25.0) && status != Status::Fail {
                status = Status::Warn;
            }

            rows.push(PerBaseQualRow {
                group: *g,
                mean,
                median,
                lower_quartile: lq,
                upper_quartile: uq,
                p10,
                p90,
            });
        }
        (rows, status)
    }

    fn build_per_base_content(&self, groups: &[BaseGroup]) -> (Vec<PerBaseContentRow>, Status) {
        let mut rows = Vec::with_capacity(groups.len());
        let mut status = Status::Pass;
        for g in groups {
            let mut acc = BaseCounts::zero();
            for i in g.range() {
                if let Some(bc) = self.per_pos_base.get(i) {
                    acc.add_assign(bc);
                }
            }
            let denom = acc.called();
            if denom == 0 {
                continue;
            }
            let d = denom as f64;
            let a = acc.a as f64 * 100.0 / d;
            let c = acc.c as f64 * 100.0 / d;
            let g_pct = acc.g as f64 * 100.0 / d;
            let t = acc.t as f64 * 100.0 / d;

            // FastQC compares the two complementary pairs, not the distance
            // from 25 %: a GC-rich genome is not a per-base bias.
            let bias = (a - t).abs().max((g_pct - c).abs());
            if bias > 20.0 {
                status = Status::Fail;
            } else if bias > 10.0 && status != Status::Fail {
                status = Status::Warn;
            }

            rows.push(PerBaseContentRow {
                group: *g,
                g: g_pct,
                a,
                t,
                c,
            });
        }
        (rows, status)
    }

    fn build_per_base_n(&self, groups: &[BaseGroup]) -> (Vec<PerBaseNRow>, Status) {
        let mut rows = Vec::with_capacity(groups.len());
        let mut status = Status::Pass;
        for g in groups {
            let mut acc = BaseCounts::zero();
            for i in g.range() {
                if let Some(bc) = self.per_pos_base.get(i) {
                    acc.add_assign(bc);
                }
            }
            let total = acc.total();
            if total == 0 {
                continue;
            }
            let n_percent = acc.n as f64 * 100.0 / total as f64;
            if n_percent > 20.0 {
                status = Status::Fail;
            } else if n_percent > 5.0 && status != Status::Fail {
                status = Status::Warn;
            }
            rows.push(PerBaseNRow {
                group: *g,
                n_percent,
            });
        }
        (rows, status)
    }

    /// FastQC uses the mode of the per-read mean quality distribution.
    fn per_seq_quality_status(&self, mode: Mode) -> Status {
        if self.total_reads == 0 {
            return Status::Pass;
        }
        match mode {
            Mode::Short => {
                let mode_q = self
                    .per_seq_mean_q_hist
                    .iter()
                    .enumerate()
                    .rev()
                    .max_by_key(|(_, c)| **c)
                    .map(|(q, _)| q)
                    .unwrap_or(0);
                if mode_q <= 20 {
                    Status::Fail
                } else if mode_q <= 27 {
                    Status::Warn
                } else {
                    Status::Pass
                }
            }
            Mode::Long => {
                // Long-read basecallers top out far below Q27; judge on the
                // median read quality instead.
                let median = quantile_from_hist(&self.per_seq_mean_q_hist, 0.5);
                if median < 7 {
                    Status::Fail
                } else if median < 10 {
                    Status::Warn
                } else {
                    Status::Pass
                }
            }
        }
    }

    /// FastQC compares the observed GC distribution with a normal fitted to it
    /// and sums the absolute deviation.
    fn per_seq_gc_status(&self) -> Status {
        let total: u64 = self.per_seq_gc_hist.iter().sum();
        if total == 0 {
            return Status::Pass;
        }
        let total_f = total as f64;
        let mean: f64 = self
            .per_seq_gc_hist
            .iter()
            .enumerate()
            .map(|(gc, &c)| gc as f64 * c as f64)
            .sum::<f64>()
            / total_f;
        let variance: f64 = self
            .per_seq_gc_hist
            .iter()
            .enumerate()
            .map(|(gc, &c)| {
                let d = gc as f64 - mean;
                d * d * c as f64
            })
            .sum::<f64>()
            / total_f;
        let stdev = variance.sqrt();
        if stdev <= 0.0 {
            return Status::Pass;
        }

        let norm = 1.0 / (stdev * (2.0 * std::f64::consts::PI).sqrt());
        let mut deviation = 0.0f64;
        for (gc, &count) in self.per_seq_gc_hist.iter().enumerate() {
            let d = gc as f64 - mean;
            let theoretical = total_f * norm * (-(d * d) / (2.0 * variance)).exp();
            deviation += (count as f64 - theoretical).abs();
        }
        let deviation_percent = deviation * 100.0 / total_f;

        if deviation_percent > 30.0 {
            Status::Fail
        } else if deviation_percent > 15.0 {
            Status::Warn
        } else {
            Status::Pass
        }
    }

    /// FastQC: warn when reads have differing lengths, fail when any read is empty.
    fn length_status(&self) -> Status {
        if self.total_reads == 0 {
            return Status::Pass;
        }
        if self.empty_reads > 0 {
            return Status::Fail;
        }
        if self.length_hist.distinct_lengths() > 1 {
            Status::Warn
        } else {
            Status::Pass
        }
    }

    fn build_duplication(
        &self,
        ctx: &FinalizeContext,
    ) -> (
        Vec<DuplicationRow>,
        Vec<OverrepRow>,
        Status,
        Status,
        Option<DuplicationSummary>,
    ) {
        if ctx.mode != Mode::Short || self.seq_tracker.is_empty() {
            return (Vec::new(), Vec::new(), Status::Pass, Status::Pass, None);
        }

        let summary = self.seq_tracker.duplication_levels();
        let dup_status = duplication::duplication_status(summary.total_deduplicated_percent);

        let window_reads = summary.observed_reads.max(1);
        let mut overrep = Vec::new();
        for (seq, count, _first) in self.seq_tracker.entries() {
            if u128::from(count) * 1_000 >= u128::from(window_reads) {
                overrep.push(OverrepRow {
                    sequence: String::from_utf8_lossy(seq).into_owned(),
                    count,
                    percent: count as f64 * 100.0 / window_reads as f64,
                    source: self.cfg.contaminants.classify(seq),
                });
            }
        }
        overrep.sort_unstable_by(|a, b| {
            b.count
                .cmp(&a.count)
                .then_with(|| a.sequence.cmp(&b.sequence))
        });
        let overrep_status =
            overrepresented::overrep_status(overrep.first().map_or(0.0, |row| row.percent));

        (
            summary.rows.clone(),
            overrep,
            dup_status,
            overrep_status,
            Some(summary),
        )
    }

    fn build_adapter_content(
        &self,
        ctx: &FinalizeContext,
        groups: &[BaseGroup],
    ) -> (Vec<AdapterRow>, Status) {
        let stride = self.adapter_stride;
        let mut rows = Vec::new();
        let mut max_percent = 0.0f64;

        if ctx.mode == Mode::Long {
            let total = self.total_reads.max(1);
            let values: Vec<f64> = self
                .adapter_reads_any
                .iter()
                .map(|&c| c as f64 * 100.0 / total as f64)
                .collect();
            for &v in &values {
                if v > max_percent {
                    max_percent = v;
                }
            }
            rows.push(AdapterRow {
                group: BaseGroup { start: 1, end: 1 },
                values,
            });
            return (rows, adapter_content::adapter_status(max_percent));
        }

        if stride == 0 {
            return (rows, Status::Pass);
        }
        let positions_tracked = self.adapter_counts.len() / stride;
        let reaching = self.reads_reaching(positions_tracked);

        for g in groups {
            let mut values = vec![0.0f64; stride];
            let mut positions = 0u32;
            for i in g.range() {
                if i >= positions_tracked {
                    continue;
                }
                // Reads that never reach this position cannot carry adapter
                // there, so they must not dilute the percentage.
                let denom = reaching.get(i).copied().unwrap_or(0);
                if denom == 0 {
                    continue;
                }
                positions += 1;
                let row = &self.adapter_counts[i * stride..(i + 1) * stride];
                for (v, &count) in values.iter_mut().zip(row) {
                    *v += count as f64 * 100.0 / denom as f64;
                }
            }
            if positions == 0 {
                continue;
            }
            for v in values.iter_mut() {
                *v /= positions as f64;
                if *v > max_percent {
                    max_percent = *v;
                }
            }
            rows.push(AdapterRow { group: *g, values });
        }

        (rows, adapter_content::adapter_status(max_percent))
    }

    fn build_long_length(&self) -> LongLengthSummary {
        let sorted = self.length_hist.iter_sorted();
        let labels = [
            "1-9",
            "10-99",
            "100-999",
            "1k-9.9k",
            "10k-99k",
            "100k-999k",
            "1M-9.9M",
            "10M+",
        ];
        let mut bins = [0u64; 8];
        for (len, count) in &sorted {
            bins[long_len_bin(*len)] += count;
        }
        let total_bases = self.length_hist.total_bases();
        let mean = if self.total_reads == 0 {
            0.0
        } else {
            total_bases as f64 / self.total_reads as f64
        };

        LongLengthSummary {
            bins: labels
                .iter()
                .enumerate()
                .map(|(i, l)| ((*l).to_string(), bins[i]))
                .collect(),
            mean,
            median: self.length_hist.quantile(0.5),
            n50: self.length_hist.nxx(0.5),
            n90: self.length_hist.nxx(0.9),
            min: self.length_hist.min_len(),
            max: self.length_hist.max_len(),
            total_bases,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Statuses {
    pub basic: Status,
    pub per_base_qual: Status,
    pub per_tile: Status,
    pub per_seq_qual: Status,
    pub per_base_content: Status,
    pub per_seq_gc: Status,
    pub per_base_n: Status,
    pub length_dist: Status,
    pub duplication: Status,
    pub overrepresented: Status,
    pub adapter_content: Status,
    pub per_seq_n: Status,
    pub kmer_content: Status,
}

pub struct FinalMetrics {
    pub basic: BasicStats,
    pub groups: Vec<BaseGroup>,
    pub per_base_qual: Vec<PerBaseQualRow>,
    pub per_tile: Vec<PerTileRow>,
    pub per_tile_available: bool,
    pub per_seq_qual: Vec<PerSeqQualRow>,
    pub per_base_content: Vec<PerBaseContentRow>,
    pub per_seq_gc: Vec<PerSeqGcRow>,
    pub per_base_n: Vec<PerBaseNRow>,
    pub length_dist: Vec<LengthDistRow>,
    pub duplication: Vec<DuplicationRow>,
    pub duplication_summary: Option<DuplicationSummary>,
    pub overrepresented: Vec<OverrepRow>,
    pub adapter_content: Vec<AdapterRow>,
    pub adapter_names: Vec<String>,
    pub per_seq_n: Vec<PerSeqNRow>,
    pub long_length: Option<LongLengthSummary>,
    pub kmer_rows: Vec<KmerRow>,
    pub statuses: Statuses,
}

fn long_len_bin(len: u64) -> usize {
    match len {
        0..=9 => 0,
        10..=99 => 1,
        100..=999 => 2,
        1_000..=9_999 => 3,
        10_000..=99_999 => 4,
        100_000..=999_999 => 5,
        1_000_000..=9_999_999 => 6,
        _ => 7,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ctx(mode: Mode) -> FinalizeContext {
        FinalizeContext {
            encoding: Encoding::Sanger,
            file_name: "test.fastq".into(),
            sample_name: "test".into(),
            mode,
            nogroup: false,
            kmer_enabled: false,
            dup_window_reads: 0,
            dup_window_complete: true,
            ignored_modules: Vec::new(),
        }
    }

    fn agg(mode: Mode) -> Agg {
        let mut cfg = AggConfig::defaults_for(mode);
        cfg.dup_window = 0;
        Agg::new(cfg)
    }

    fn q(score: u8, len: usize) -> Vec<u8> {
        vec![33 + score; len]
    }

    #[test]
    fn gc_rich_library_is_not_a_per_base_bias() {
        // 67 % GC, but each read is paired with its complement so A==T and
        // G==C at every position. FastQC calls this a pass; the old
        // distance-from-25 % rule called it a fail.
        let mut a = agg(Mode::Short);
        let seq = b"GCGCATGCGCAT";
        let complement = b"CGCGTACGCGTA";
        let qual = q(35, seq.len());
        for i in 0..1000u64 {
            let read = if i % 2 == 0 { seq } else { complement };
            a.update_read(&ReadView::new(read, &qual), 33, i);
        }
        let m = a.finalize(&ctx(Mode::Short));
        assert_eq!(m.statuses.per_base_content, Status::Pass);
    }

    #[test]
    fn strand_bias_is_still_caught() {
        let mut a = agg(Mode::Short);
        let seq = b"AAAAAAAAAAAA";
        let qual = q(35, seq.len());
        for i in 0..100u64 {
            a.update_read(&ReadView::new(seq, &qual), 33, i);
        }
        let m = a.finalize(&ctx(Mode::Short));
        assert_eq!(m.statuses.per_base_content, Status::Fail);
    }

    #[test]
    fn unique_library_passes_duplication() {
        let mut a = agg(Mode::Short);
        let qual = q(35, 40);
        for i in 0..5000u64 {
            let seq = format!("{:A>40}", i);
            a.update_read(&ReadView::new(seq.as_bytes(), &qual), 33, i);
        }
        let m = a.finalize(&ctx(Mode::Short));
        assert_eq!(m.statuses.duplication, Status::Pass);
        let s = m.duplication_summary.unwrap();
        assert!((s.total_deduplicated_percent - 100.0).abs() < 1e-9);
    }

    #[test]
    fn duplicated_library_fails() {
        let mut a = agg(Mode::Short);
        let qual = q(35, 12);
        for i in 0..1000u64 {
            a.update_read(&ReadView::new(b"ACGTACGTACGT", &qual), 33, i);
        }
        let m = a.finalize(&ctx(Mode::Short));
        assert_eq!(m.statuses.duplication, Status::Fail);
        assert_eq!(m.statuses.overrepresented, Status::Fail);
        assert_eq!(m.overrepresented.len(), 1);
        assert!((m.overrepresented[0].percent - 100.0).abs() < 1e-9);
    }

    #[test]
    fn empty_reads_fail_length_module() {
        let mut a = agg(Mode::Short);
        let qual = q(35, 10);
        a.update_read(&ReadView::new(b"ACGTACGTAC", &qual), 33, 0);
        a.update_read(&ReadView::new(b"", b""), 33, 1);
        let m = a.finalize(&ctx(Mode::Short));
        assert_eq!(m.statuses.length_dist, Status::Fail);
        assert_eq!(m.basic.total_sequences, 2);
    }

    #[test]
    fn variable_lengths_warn() {
        let mut a = agg(Mode::Short);
        a.update_read(&ReadView::new(b"ACGTACGTAC", &q(35, 10)), 33, 0);
        a.update_read(&ReadView::new(b"ACGTACGT", &q(35, 8)), 33, 1);
        let m = a.finalize(&ctx(Mode::Short));
        assert_eq!(m.statuses.length_dist, Status::Warn);
    }

    #[test]
    fn per_base_quality_uses_quartiles() {
        let mut a = agg(Mode::Short);
        for i in 0..101u64 {
            a.update_read(&ReadView::new(b"ACGT", &q(30, 4)), 33, i);
        }
        let m = a.finalize(&ctx(Mode::Short));
        assert_eq!(m.statuses.per_base_qual, Status::Pass);

        let mut b = agg(Mode::Short);
        for i in 0..101u64 {
            b.update_read(&ReadView::new(b"ACGT", &q(4, 4)), 33, i);
        }
        let m2 = b.finalize(&ctx(Mode::Short));
        assert_eq!(m2.statuses.per_base_qual, Status::Fail);
    }

    #[test]
    fn merge_matches_single_threaded_result() {
        let qual = q(30, 12);
        let seqs: Vec<Vec<u8>> = (0..200u64)
            .map(|i| format!("{:A>12}", i % 50).into_bytes())
            .collect();

        let mut single = agg(Mode::Short);
        for (i, s) in seqs.iter().enumerate() {
            single.update_read(&ReadView::new(s, &qual), 33, i as u64);
        }

        let mut first = agg(Mode::Short);
        let mut second = agg(Mode::Short);
        for (i, s) in seqs.iter().enumerate() {
            if i % 2 == 0 {
                first.update_read(&ReadView::new(s, &qual), 33, i as u64);
            } else {
                second.update_read(&ReadView::new(s, &qual), 33, i as u64);
            }
        }
        first.merge(&second);

        let a = single.finalize(&ctx(Mode::Short));
        let b = first.finalize(&ctx(Mode::Short));
        assert_eq!(a.basic.total_sequences, b.basic.total_sequences);
        assert_eq!(a.basic.gc_percent, b.basic.gc_percent);
        assert_eq!(
            a.duplication_summary.unwrap().total_deduplicated_percent,
            b.duplication_summary.unwrap().total_deduplicated_percent
        );
        assert_eq!(a.overrepresented.len(), b.overrepresented.len());
    }

    #[test]
    fn long_mode_quality_uses_error_probabilities() {
        let mut a = agg(Mode::Long);
        // Half the bases at Q30, half at Q10. The arithmetic mean would be 20;
        // the error-probability mean is about 13.
        let mut qual = vec![33 + 30u8; 50];
        qual.extend(vec![33 + 10u8; 50]);
        let seq = vec![b'A'; 100];
        a.update_read(&ReadView::new(&seq, &qual), 33, 0);
        let m = a.finalize(&ctx(Mode::Long));
        let reported = m.per_seq_qual[0].mean_q;
        assert!(
            (12..=14).contains(&reported),
            "expected ~13, got {reported}"
        );
    }

    #[test]
    fn long_mode_n50_is_exact() {
        let mut a = agg(Mode::Long);
        for (i, len) in [1000usize, 5000, 5000].iter().enumerate() {
            let seq = vec![b'A'; *len];
            let qual = q(15, *len);
            a.update_read(&ReadView::new(&seq, &qual), 33, i as u64);
        }
        let m = a.finalize(&ctx(Mode::Long));
        let ll = m.long_length.unwrap();
        assert_eq!(ll.n50, 5000);
        assert_eq!(ll.max, 5000);
        assert_eq!(ll.min, 1000);
        assert_eq!(ll.total_bases, 11_000);
    }

    #[test]
    fn adapter_readthrough_is_detected() {
        let mut a = agg(Mode::Short);
        let mut seq = vec![b'C'; 30];
        seq.extend_from_slice(b"AGATCGGAAGAGCACACG");
        let qual = q(35, seq.len());
        for i in 0..100u64 {
            a.update_read(&ReadView::new(&seq, &qual), 33, i);
        }
        let m = a.finalize(&ctx(Mode::Short));
        assert_eq!(m.statuses.adapter_content, Status::Fail);
        let last = m.adapter_content.last().unwrap();
        assert!(last.values[0] > 90.0);
    }

    #[test]
    fn base_groups_bound_row_count() {
        let mut a = agg(Mode::Short);
        let seq = vec![b'A'; 2000];
        let qual = q(35, 2000);
        a.update_read(&ReadView::new(&seq, &qual), 33, 0);
        let m = a.finalize(&ctx(Mode::Short));
        assert!(
            m.per_base_qual.len() < 100,
            "expected grouped rows, got {}",
            m.per_base_qual.len()
        );
        assert_eq!(m.per_base_qual[0].group.label(), "1");
        assert_eq!(m.per_base_qual[9].group.label(), "10-49");
    }

    #[test]
    fn dup_window_bounds_tracking() {
        let mut cfg = AggConfig::defaults_for(Mode::Short);
        cfg.dup_window = 10;
        cfg.per_tile_enabled = false;
        let mut a = Agg::new(cfg);
        let qual = q(35, 4);
        for i in 0..100u64 {
            a.update_read(&ReadView::new(b"ACGT", &qual), 33, i);
        }
        let m = a.finalize(&ctx(Mode::Short));
        let s = m.duplication_summary.unwrap();
        assert_eq!(s.observed_reads, 10);
        assert_eq!(m.basic.total_sequences, 100);
    }
}
