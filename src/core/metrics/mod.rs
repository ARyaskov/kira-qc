use crate::core::fastq::ReadView;
use crate::core::model::{
    Encoding, FinalizeContext, MAX_Q, Mode, QualHist, Status, quantile_from_hist,
};
use crate::simd;
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

pub use adapter_content::{ADAPTERS, AdapterRow};
pub use basic::BasicStats;
pub use duplication::{DupLevel, DuplicationRow, SpaceSaving as DupSpaceSaving};
pub use kmer_content::KmerRow;
pub use length_dist::LengthDistRow;
pub use overrepresented::{OverrepRow, SpaceSavingSeq};
pub use per_base_content::PerBaseContentRow;
pub use per_base_n::PerBaseNRow;
pub use per_base_qual::PerBaseQualRow;
pub use per_seq_gc::PerSeqGcRow;
pub use per_seq_n::PerSeqNRow;
pub use per_seq_qual::PerSeqQualRow;

#[derive(Clone, Debug)]
pub struct BaseCounts {
    pub a: u64,
    pub c: u64,
    pub g: u64,
    pub t: u64,
    pub n: u64,
}

/// Optional per-worker timing breakdown, gated on `KIRA_STATS=1`.
#[derive(Clone, Debug, Default)]
pub struct UpdateTimings {
    pub kmer: Duration,
    pub kmer_updates: u64,
}

impl BaseCounts {
    fn zero() -> Self {
        Self {
            a: 0,
            c: 0,
            g: 0,
            t: 0,
            n: 0,
        }
    }

    fn add_assign(&mut self, other: &BaseCounts) {
        self.a += other.a;
        self.c += other.c;
        self.g += other.g;
        self.t += other.t;
        self.n += other.n;
    }
}

#[derive(Clone, Debug)]
pub struct Agg {
    pub mode: Mode,
    pub total_reads: u64,
    pub total_bases: u64,
    pub gc_bases: u64,
    pub n_bases: u64,
    pub min_len: u32,
    pub max_len: u32,
    pub per_pos_qual: Vec<QualHist>,
    pub per_pos_base: Vec<BaseCounts>,
    pub per_seq_mean_q_hist: Vec<u64>,
    pub per_seq_gc_hist: [u64; 101],
    pub length_hist: Vec<u64>,
    pub reads_mean_q_lt_20: u64,
    pub dup_space: DupSpaceSaving,
    pub overrep_space: SpaceSavingSeq,
    pub adapter_counts: Vec<[u64; ADAPTERS.len()]>,
    pub per_seq_n_hist: [u64; 101],
    pub reads_n_gt10: u64,
    pub reads_n_gt20: u64,
    pub adapter_reads_any: [u64; ADAPTERS.len()],
    pub long_len_bins: [u64; 8],
    pub kmer_bins: Vec<kmer_content::BinCounts>,
    pub kmer_bin_counts: [u64; kmer_content::BINS],
    pub kmer_total: u64,
}

impl Agg {
    pub fn new(mode: Mode) -> Self {
        Self {
            mode,
            total_reads: 0,
            total_bases: 0,
            gc_bases: 0,
            n_bases: 0,
            min_len: u32::MAX,
            max_len: 0,
            per_pos_qual: Vec::new(),
            per_pos_base: Vec::new(),
            per_seq_mean_q_hist: vec![0u64; MAX_Q + 1],
            per_seq_gc_hist: [0u64; 101],
            length_hist: Vec::new(),
            reads_mean_q_lt_20: 0,
            dup_space: DupSpaceSaving::new(),
            overrep_space: SpaceSavingSeq::new(),
            adapter_counts: Vec::new(),
            per_seq_n_hist: [0u64; 101],
            reads_n_gt10: 0,
            reads_n_gt20: 0,
            adapter_reads_any: [0u64; ADAPTERS.len()],
            long_len_bins: [0u64; 8],
            kmer_bins: if mode == Mode::Short {
                (0..kmer_content::BINS)
                    .map(|_| kmer_content::BinCounts::new())
                    .collect()
            } else {
                Vec::new()
            },
            kmer_bin_counts: [0u64; kmer_content::BINS],
            kmer_total: 0,
        }
    }

    pub fn update_read(&mut self, read: &ReadView<'_>, phred_offset: u8) {
        self.update_read_with_timings(read, phred_offset, None);
    }

    pub fn update_read_timed(
        &mut self,
        read: &ReadView<'_>,
        phred_offset: u8,
        timing: &mut UpdateTimings,
    ) {
        self.update_read_with_timings(read, phred_offset, Some(timing));
    }

    fn update_read_with_timings(
        &mut self,
        read: &ReadView<'_>,
        phred_offset: u8,
        timing: Option<&mut UpdateTimings>,
    ) {
        let len = read.seq.len();
        if len == 0 {
            return;
        }

        self.total_reads += 1;
        self.total_bases += len as u64;
        if (len as u32) > self.max_len {
            self.max_len = len as u32;
        }
        if (len as u32) < self.min_len {
            self.min_len = len as u32;
        }

        let (_a, c, g, _t, n) = simd::count_bases(read.seq);
        let gc = c as u64 + g as u64;
        let n_count = n as u64;
        self.gc_bases += gc;
        self.n_bases += n_count;

        match self.mode {
            Mode::Short => self.update_short(read, len, timing),
            Mode::Long => self.update_long(read, len, n_count),
        }

        let sum_q = simd::sum_qual(read.qual, phred_offset) as u64;
        let mean_q = (sum_q + (len as u64 / 2)) / len as u64;
        if mean_q < 20 {
            self.reads_mean_q_lt_20 += 1;
        }
        let mean_q_bin = (mean_q as usize).min(MAX_Q);
        self.per_seq_mean_q_hist[mean_q_bin] += 1;

        let gc_percent = ((gc * 100) + (len as u64 / 2)) / len as u64;
        let gc_bin = (gc_percent as usize).min(100);
        self.per_seq_gc_hist[gc_bin] += 1;

        if self.mode == Mode::Short {
            self.update_per_position(read, phred_offset, len);
        }
    }

    fn update_short(
        &mut self,
        read: &ReadView<'_>,
        len: usize,
        timing: Option<&mut UpdateTimings>,
    ) {
        if self.per_pos_qual.len() < len {
            self.per_pos_qual.resize(len, [0u64; MAX_Q + 1]);
        }
        if self.per_pos_base.len() < len {
            self.per_pos_base.resize(len, BaseCounts::zero());
        }
        if self.length_hist.len() <= len {
            self.length_hist.resize(len + 1, 0);
        }
        self.length_hist[len] += 1;

        let key = duplication::hash_seq(read.seq);
        self.dup_space.add(key, 1);

        let key2 = overrepresented::hash_seq(read.seq);
        self.overrep_space.add(key2, read.seq, 1);

        if self.adapter_counts.len() < len {
            self.adapter_counts.resize(len, [0u64; ADAPTERS.len()]);
        }
        adapter_content::scan(read.seq, &mut self.adapter_counts);

        #[cfg(not(feature = "no-kmer"))]
        if len >= kmer_content::K {
            kmer_content::update_kmers(
                read.seq,
                len,
                &mut self.kmer_bins,
                &mut self.kmer_bin_counts,
                &mut self.kmer_total,
                timing,
            );
        }
        #[cfg(feature = "no-kmer")]
        let _ = timing;
    }

    fn update_long(&mut self, read: &ReadView<'_>, len: usize, n_count: u64) {
        let bin = long_len_bin(len as u64);
        self.long_len_bins[bin] += 1;

        let n_percent = ((n_count * 100) + (len as u64 / 2)) / len as u64;
        let n_bin = n_percent.min(100) as usize;
        self.per_seq_n_hist[n_bin] += 1;
        if n_percent > 20 {
            self.reads_n_gt20 += 1;
        } else if n_percent > 10 {
            self.reads_n_gt10 += 1;
        }

        let mut hits = [false; ADAPTERS.len()];
        adapter_content::scan_any(read.seq, &mut hits);
        for (i, hit) in hits.iter().enumerate() {
            if *hit {
                self.adapter_reads_any[i] += 1;
            }
        }
    }

    fn update_per_position(&mut self, read: &ReadView<'_>, phred_offset: u8, len: usize) {
        // LUT slot: 0=A 1=C 2=G 3=T 4=N 5=other.
        const fn build_lut() -> [u8; 256] {
            let mut lut = [5u8; 256];
            lut[b'A' as usize] = 0;
            lut[b'C' as usize] = 1;
            lut[b'G' as usize] = 2;
            lut[b'T' as usize] = 3;
            lut[b'N' as usize] = 4;
            lut
        }
        static BASE_LUT: [u8; 256] = build_lut();

        let seq = &read.seq[..len];
        let qual = &read.qual[..len];
        let per_pos_base = &mut self.per_pos_base[..len];
        let per_pos_qual = &mut self.per_pos_qual[..len];

        for i in 0..len {
            let folded = seq[i] & 0xDF;
            let slot = BASE_LUT[folded as usize];
            let bc = &mut per_pos_base[i];
            match slot {
                0 => bc.a += 1,
                1 => bc.c += 1,
                2 => bc.g += 1,
                3 => bc.t += 1,
                4 => bc.n += 1,
                _ => {}
            }

            let q_raw = qual[i].saturating_sub(phred_offset) as usize;
            let q_bin = q_raw.min(MAX_Q);
            per_pos_qual[i][q_bin] += 1;
        }
    }

    pub fn merge(&mut self, other: &Agg) {
        self.total_reads += other.total_reads;
        self.total_bases += other.total_bases;
        self.gc_bases += other.gc_bases;
        self.n_bases += other.n_bases;
        if other.min_len < self.min_len {
            self.min_len = other.min_len;
        }
        if other.max_len > self.max_len {
            self.max_len = other.max_len;
        }
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
        self.reads_mean_q_lt_20 += other.reads_mean_q_lt_20;
        self.reads_n_gt10 += other.reads_n_gt10;
        self.reads_n_gt20 += other.reads_n_gt20;

        match self.mode {
            Mode::Short => {
                if self.length_hist.len() < other.length_hist.len() {
                    self.length_hist.resize(other.length_hist.len(), 0);
                }
                for i in 0..other.length_hist.len() {
                    self.length_hist[i] += other.length_hist[i];
                }
                self.dup_space.merge(&other.dup_space);
                self.overrep_space.merge(&other.overrep_space);
                if self.adapter_counts.len() < other.adapter_counts.len() {
                    self.adapter_counts
                        .resize(other.adapter_counts.len(), [0u64; ADAPTERS.len()]);
                }
                for (target, row) in self.adapter_counts.iter_mut().zip(&other.adapter_counts) {
                    for (t, &s) in target.iter_mut().zip(row.iter()) {
                        *t += s;
                    }
                }
                #[cfg(not(feature = "no-kmer"))]
                {
                    for b in 0..kmer_content::BINS {
                        self.kmer_bins[b].merge(&other.kmer_bins[b]);
                        self.kmer_bin_counts[b] += other.kmer_bin_counts[b];
                    }
                    self.kmer_total += other.kmer_total;
                }
            }
            Mode::Long => {
                for i in 0..self.long_len_bins.len() {
                    self.long_len_bins[i] += other.long_len_bins[i];
                }
                for i in 0..ADAPTERS.len() {
                    self.adapter_reads_any[i] += other.adapter_reads_any[i];
                }
            }
        }
    }

    pub fn finalize(&self, ctx: &FinalizeContext) -> FinalMetrics {
        let min_len = if self.total_reads == 0 {
            0
        } else {
            self.min_len
        };
        let max_len = if self.total_reads == 0 {
            0
        } else {
            self.max_len
        };
        let gc_percent = ((self.gc_bases * 100) + (self.total_bases / 2))
            .checked_div(self.total_bases)
            .unwrap_or(0) as u32;

        let encoding_str = match ctx.encoding {
            Encoding::Sanger => "Sanger / Illumina 1.9",
            Encoding::Illumina15 => "Illumina 1.5",
        };

        let basic = BasicStats {
            file_type: "Conventional base calls",
            encoding: encoding_str,
            total_sequences: self.total_reads,
            filtered_sequences: 0,
            min_len,
            max_len,
            gc_percent,
        };

        let mut per_base_qual = Vec::new();
        if ctx.mode == Mode::Short {
            per_base_qual.reserve(self.per_pos_qual.len());
            for (i, hist) in self.per_pos_qual.iter().enumerate() {
                let mut total: u64 = 0;
                let mut sum: u64 = 0;
                for (q, &c) in hist.iter().enumerate() {
                    total += c;
                    sum += c * q as u64;
                }
                let mean = if total == 0 {
                    0.0
                } else {
                    sum as f64 / total as f64
                };
                let median = quantile_from_hist(hist, 0.5);
                let lq = quantile_from_hist(hist, 0.25);
                let uq = quantile_from_hist(hist, 0.75);
                let p10 = quantile_from_hist(hist, 0.10);
                let p90 = quantile_from_hist(hist, 0.90);
                per_base_qual.push(PerBaseQualRow {
                    base: i + 1,
                    mean,
                    median,
                    lower_quartile: lq,
                    upper_quartile: uq,
                    p10,
                    p90,
                });
            }
        }

        let mut per_seq_qual = Vec::new();
        for (q, &count) in self.per_seq_mean_q_hist.iter().enumerate() {
            if count > 0 {
                per_seq_qual.push(PerSeqQualRow {
                    mean_q: q as u8,
                    count,
                });
            }
        }

        let mut per_base_content = Vec::new();
        let mut max_deviation: f64 = 0.0;
        if ctx.mode == Mode::Short {
            per_base_content.reserve(self.per_pos_base.len());
            for (i, bc) in self.per_pos_base.iter().enumerate() {
                let denom = bc.a + bc.c + bc.g + bc.t;
                let (g, a, t, c) = if denom == 0 {
                    (0.0, 0.0, 0.0, 0.0)
                } else {
                    let d = denom as f64;
                    (
                        bc.g as f64 * 100.0 / d,
                        bc.a as f64 * 100.0 / d,
                        bc.t as f64 * 100.0 / d,
                        bc.c as f64 * 100.0 / d,
                    )
                };
                if denom > 0 {
                    let devs = [
                        (g - 25.0).abs(),
                        (a - 25.0).abs(),
                        (t - 25.0).abs(),
                        (c - 25.0).abs(),
                    ];
                    for d in devs {
                        if d > max_deviation {
                            max_deviation = d;
                        }
                    }
                }
                per_base_content.push(PerBaseContentRow {
                    base: i + 1,
                    g,
                    a,
                    t,
                    c,
                });
            }
        }

        let mut per_seq_gc = Vec::new();
        for (gc, &count) in self.per_seq_gc_hist.iter().enumerate() {
            if count > 0 {
                per_seq_gc.push(PerSeqGcRow {
                    gc: gc as u8,
                    count,
                });
            }
        }

        let mut per_base_n = Vec::new();
        let mut max_n_percent: f64 = 0.0;
        if ctx.mode == Mode::Short {
            per_base_n.reserve(self.per_pos_base.len());
            for (i, bc) in self.per_pos_base.iter().enumerate() {
                let total = bc.a + bc.c + bc.g + bc.t + bc.n;
                let n_percent = if total == 0 {
                    0.0
                } else {
                    bc.n as f64 * 100.0 / total as f64
                };
                if n_percent > max_n_percent {
                    max_n_percent = n_percent;
                }
                per_base_n.push(PerBaseNRow {
                    base: i + 1,
                    n_percent,
                });
            }
        }

        let mut length_dist = Vec::new();
        let mut long_length = None;
        if ctx.mode == Mode::Short {
            for (len, &count) in self.length_hist.iter().enumerate() {
                if len > 0 && count > 0 {
                    length_dist.push(LengthDistRow { length: len, count });
                }
            }
        } else {
            long_length = Some(build_long_length(
                &self.long_len_bins,
                self.total_reads,
                self.total_bases,
                min_len,
                max_len,
            ));
        }

        let mut per_base_qual_status = Status::Pass;
        let mut per_seq_qual_status = Status::Pass;
        if ctx.mode == Mode::Short {
            for row in &per_base_qual {
                if row.median < 20 {
                    per_base_qual_status = Status::Fail;
                    break;
                }
                if row.median < 25 {
                    per_base_qual_status = Status::Warn;
                }
            }
            if self.total_reads > 0 {
                let low = self.reads_mean_q_lt_20 as f64 / self.total_reads as f64 * 100.0;
                if low > 20.0 {
                    per_seq_qual_status = Status::Fail;
                } else if low > 10.0 {
                    per_seq_qual_status = Status::Warn;
                }
            }
        } else {
            let median = quantile_from_hist(&self.per_seq_mean_q_hist, 0.5);
            if median < 7 {
                per_seq_qual_status = Status::Fail;
            } else if median < 10 {
                per_seq_qual_status = Status::Warn;
            }
        }

        let per_base_content_status = if ctx.mode == Mode::Short {
            if max_deviation > 20.0 {
                Status::Fail
            } else if max_deviation > 10.0 {
                Status::Warn
            } else {
                Status::Pass
            }
        } else {
            Status::Pass
        };

        let per_base_n_status = if ctx.mode == Mode::Short {
            if max_n_percent > 20.0 {
                Status::Fail
            } else if max_n_percent > 5.0 {
                Status::Warn
            } else {
                Status::Pass
            }
        } else {
            Status::Pass
        };

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

        let mut dup_counts = [0u64; 7];
        let mut tracked_total: u64 = 0;
        for e in self.dup_space.entries() {
            tracked_total += e.count;
            let idx = if e.count >= 7 {
                6
            } else {
                (e.count as usize).saturating_sub(1)
            };
            dup_counts[idx] += e.count;
        }
        if tracked_total > self.total_reads {
            tracked_total = self.total_reads;
        }
        let unique_extra = self.total_reads.saturating_sub(tracked_total);
        dup_counts[0] += unique_extra;

        let mut duplication = Vec::new();
        let mut overrep = Vec::new();
        let mut adapter_rows = Vec::new();
        #[cfg_attr(feature = "no-kmer", allow(unused_mut))]
        let mut kmer_rows = Vec::new();
        #[cfg_attr(feature = "no-kmer", allow(unused_mut))]
        let mut kmer_status = Status::Pass;
        let total_reads = self.total_reads.max(1);
        let mut duplication_status = Status::Pass;
        let mut overrep_status = Status::Pass;
        let mut adapter_status = Status::Pass;

        if ctx.mode == Mode::Short {
            let levels = [
                DupLevel::One,
                DupLevel::Two,
                DupLevel::Three,
                DupLevel::Four,
                DupLevel::Five,
                DupLevel::Six,
                DupLevel::SevenPlus,
            ];
            for i in 0..7 {
                let rel = dup_counts[i] as f64 / total_reads as f64;
                duplication.push(DuplicationRow {
                    level: levels[i],
                    relative: rel,
                });
            }

            let duplicated_reads = total_reads.saturating_sub(dup_counts[0]);
            let duplicated_pct = duplicated_reads as f64 * 100.0 / total_reads as f64;
            duplication_status = if duplicated_pct > 80.0 {
                Status::Fail
            } else if duplicated_pct > 50.0 {
                Status::Warn
            } else {
                Status::Pass
            };

            let mut warn_hit = false;
            for e in self.overrep_space.entries() {
                if e.count == 0 {
                    continue;
                }
                let pct = e.count as f64 * 100.0 / total_reads as f64;
                if pct >= 0.1 {
                    let seq = String::from_utf8_lossy(&e.seq).to_string();
                    let source = overrepresented::classify_source(&e.seq);
                    overrep.push(OverrepRow {
                        sequence: seq,
                        count: e.count,
                        percent: pct,
                        source,
                    });
                    overrep_status = Status::Fail;
                } else if pct >= 0.05 {
                    warn_hit = true;
                }
            }
            if overrep_status == Status::Pass && warn_hit {
                overrep_status = Status::Warn;
            }
            overrep.sort_by(|a, b| {
                b.count
                    .cmp(&a.count)
                    .then_with(|| a.sequence.cmp(&b.sequence))
            });

            for (i, row) in self.adapter_counts.iter().enumerate() {
                let mut values = [0.0f64; ADAPTERS.len()];
                for j in 0..ADAPTERS.len() {
                    let pct = row[j] as f64 * 100.0 / total_reads as f64;
                    values[j] = pct;
                    if pct > 10.0 {
                        adapter_status = Status::Fail;
                    } else if pct > 5.0 && adapter_status != Status::Fail {
                        adapter_status = Status::Warn;
                    }
                }
                adapter_rows.push(AdapterRow {
                    position: i + 1,
                    values,
                });
            }

            #[cfg(not(feature = "no-kmer"))]
            if self.kmer_total > 0 {
                let (rows, status) = kmer_content::build_rows(
                    &self.kmer_bins,
                    &self.kmer_bin_counts,
                    self.kmer_total,
                );
                kmer_rows = rows;
                kmer_status = status;
            }
        } else {
            let mut values = [0.0f64; ADAPTERS.len()];
            for (i, &count) in self.adapter_reads_any.iter().enumerate() {
                let pct = count as f64 * 100.0 / total_reads as f64;
                values[i] = pct;
                if pct > 10.0 {
                    adapter_status = Status::Fail;
                } else if pct > 5.0 && adapter_status != Status::Fail {
                    adapter_status = Status::Warn;
                }
            }
            adapter_rows.push(AdapterRow {
                position: 1,
                values,
            });
        }

        let statuses = Statuses {
            basic: Status::Pass,
            per_base_qual: per_base_qual_status,
            per_seq_qual: per_seq_qual_status,
            per_base_content: per_base_content_status,
            per_seq_gc: Status::Pass,
            per_base_n: per_base_n_status,
            length_dist: Status::Pass,
            duplication: duplication_status,
            overrepresented: overrep_status,
            adapter_content: adapter_status,
            per_seq_n: per_seq_n_status,
            kmer_content: kmer_status,
        };

        FinalMetrics {
            basic,
            per_base_qual,
            per_seq_qual,
            per_base_content,
            per_seq_gc,
            per_base_n,
            length_dist,
            duplication,
            overrepresented: overrep,
            adapter_content: adapter_rows,
            per_seq_n,
            long_length,
            kmer_rows,
            statuses,
        }
    }
}

pub struct Statuses {
    pub basic: Status,
    pub per_base_qual: Status,
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
    pub per_base_qual: Vec<PerBaseQualRow>,
    pub per_seq_qual: Vec<PerSeqQualRow>,
    pub per_base_content: Vec<PerBaseContentRow>,
    pub per_seq_gc: Vec<PerSeqGcRow>,
    pub per_base_n: Vec<PerBaseNRow>,
    pub length_dist: Vec<LengthDistRow>,
    pub duplication: Vec<DuplicationRow>,
    pub overrepresented: Vec<OverrepRow>,
    pub adapter_content: Vec<AdapterRow>,
    pub per_seq_n: Vec<PerSeqNRow>,
    pub long_length: Option<LongLengthSummary>,
    pub kmer_rows: Vec<KmerRow>,
    pub statuses: Statuses,
}

#[derive(Clone, Debug)]
pub struct LongLengthSummary {
    pub bins: [u64; 8],
    pub labels: [&'static str; 8],
    pub mean: f64,
    pub n50: u64,
    pub n90: u64,
    pub min: u32,
    pub max: u32,
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

fn build_long_length(
    bins: &[u64; 8],
    total_reads: u64,
    total_bases: u64,
    min: u32,
    max: u32,
) -> LongLengthSummary {
    let labels = [
        "1-9",
        "10-99",
        "100-999",
        "1k-9k",
        "10k-99k",
        "100k-999k",
        "1M-9M",
        "10M+",
    ];
    let mean = if total_reads == 0 {
        0.0
    } else {
        total_bases as f64 / total_reads as f64
    };
    let n50 = approx_nxx(bins, total_bases, 0.5);
    let n90 = approx_nxx(bins, total_bases, 0.9);
    LongLengthSummary {
        bins: *bins,
        labels,
        mean,
        n50,
        n90,
        min,
        max,
    }
}

fn approx_nxx(bins: &[u64; 8], total_bases: u64, frac: f64) -> u64 {
    let target = (total_bases as f64 * frac) as u64;
    let mut acc = 0u64;
    let mids: [u64; 8] = [5, 55, 550, 5_500, 55_000, 550_000, 5_500_000, 10_000_000];
    for i in (0..bins.len()).rev() {
        acc += bins[i] * mids[i];
        if acc >= target {
            return mids[i];
        }
    }
    0
}
