//! Exact k-mer counting with K=7 (4^K = 16384 buckets per position bin).

#[cfg(not(feature = "no-kmer"))]
mod real {
    use crate::core::metrics::UpdateTimings;
    use std::time::Instant;

    pub const K: usize = 7;
    pub const BINS: usize = 10;
    pub const KMER_SPACE: usize = 1 << (2 * K);
    const MAX_REPORT: usize = 50;

    #[derive(Clone, Debug)]
    pub struct KmerRow {
        pub sequence: String,
        pub count: u64,
        pub p_value: f64,
        pub obs_exp: f64,
        pub max_pos: u32,
    }

    #[derive(Clone, Debug)]
    pub struct BinCounts {
        cells: Vec<u32>,
    }

    impl Default for BinCounts {
        fn default() -> Self {
            Self::new()
        }
    }

    impl BinCounts {
        #[inline]
        pub fn new() -> Self {
            Self {
                cells: vec![0u32; KMER_SPACE],
            }
        }

        #[inline]
        pub fn add(&mut self, key: u16) {
            let c = &mut self.cells[key as usize];
            *c = c.saturating_add(1);
        }

        #[inline]
        pub fn get(&self, key: u16) -> u32 {
            self.cells[key as usize]
        }

        pub fn merge(&mut self, other: &BinCounts) {
            debug_assert_eq!(self.cells.len(), other.cells.len());
            for (a, &b) in self.cells.iter_mut().zip(other.cells.iter()) {
                *a = a.saturating_add(b);
            }
        }
    }

    /// 2-bit base encoding: first base of the k-mer in the high bits.
    #[inline]
    fn base_code(b: u8) -> Option<u64> {
        match b & 0xDF {
            b'A' => Some(0),
            b'C' => Some(1),
            b'G' => Some(2),
            b'T' => Some(3),
            _ => None,
        }
    }

    pub fn decode_kmer(mut key: u16) -> String {
        let mut buf = [b'A'; K];
        for i in (0..K).rev() {
            let bits = (key & 0x3) as u8;
            buf[i] = match bits {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            };
            key >>= 2;
        }
        // SAFETY: buf is ASCII-only by construction.
        unsafe { String::from_utf8_unchecked(buf.to_vec()) }
    }

    pub fn bin_mid_percent(bin: usize) -> u32 {
        if bin >= BINS - 1 {
            95
        } else {
            (bin as u32) * 10 + 5
        }
    }

    pub fn compute_pvalue(obs_exp: f64) -> f64 {
        if obs_exp <= 1.0 {
            return 1.0;
        }
        let score = (obs_exp - 1.0) * 2.0;
        let exp = score.min(50.0);
        10f64.powf(-exp)
    }

    pub fn select_top(rows: &mut Vec<KmerRow>) {
        rows.sort_by(|a, b| {
            b.obs_exp
                .partial_cmp(&a.obs_exp)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        if rows.len() > MAX_REPORT {
            rows.truncate(MAX_REPORT);
        }
    }

    #[inline(always)]
    fn next_bin_threshold(len: usize, bin: usize) -> usize {
        let t = (bin + 1) * 10;
        (t * len).div_ceil(100)
    }

    /// Single rolling-hash pass over the read; Ns reset the run.
    pub fn update_kmers(
        seq: &[u8],
        len: usize,
        counts: &mut [BinCounts],
        bin_counts: &mut [u64; BINS],
        total: &mut u64,
        timing: Option<&mut UpdateTimings>,
    ) {
        if len < K {
            return;
        }
        let t_total = timing.as_ref().map(|_| Instant::now());

        const MASK: u64 = (1u64 << (2 * K)) - 1;
        let mut rolling: u64 = 0;
        let mut valid_run: usize = 0;
        let mut bin: usize = 0;
        let mut next_threshold = next_bin_threshold(len, bin);

        let mut updates: u64 = 0;

        for (pos, &b) in seq.iter().enumerate() {
            match base_code(b) {
                Some(code) => {
                    rolling = ((rolling << 2) | code) & MASK;
                    if valid_run < K {
                        valid_run += 1;
                    }
                }
                None => {
                    valid_run = 0;
                    rolling = 0;
                    continue;
                }
            }
            if valid_run >= K {
                let start_pos_plus1 = pos + 2 - K;
                while bin + 1 < BINS && start_pos_plus1 >= next_threshold {
                    bin += 1;
                    next_threshold = next_bin_threshold(len, bin);
                }
                counts[bin].add(rolling as u16);
                bin_counts[bin] += 1;
                *total += 1;
                updates += 1;
            }
        }

        if let (Some(t0), Some(t)) = (t_total, timing) {
            t.kmer += t0.elapsed();
            t.kmer_updates += updates;
        }
    }

    /// Rows for k-mers whose obs/exp in any bin is ≥ 3.0 (FastQC heuristic).
    pub fn build_rows(
        counts: &[BinCounts],
        bin_counts: &[u64; BINS],
        total: u64,
    ) -> (Vec<KmerRow>, crate::core::model::Status) {
        use crate::core::model::Status;
        let mut rows = Vec::new();
        let mut status = Status::Pass;
        if total == 0 {
            return (rows, status);
        }

        for key in 0..(KMER_SPACE as u16) {
            let total_for_key: u64 = counts.iter().map(|c| c.get(key) as u64).sum();
            if total_for_key == 0 {
                continue;
            }
            let expected = total_for_key as f64 / total as f64;
            if expected == 0.0 {
                continue;
            }
            let mut max_obs = 0.0f64;
            let mut max_bin = 0usize;
            for (b, bc) in counts.iter().enumerate() {
                let denom = bin_counts[b] as f64;
                if denom == 0.0 {
                    continue;
                }
                let obs = bc.get(key) as f64 / denom;
                let obs_exp = obs / expected;
                if obs_exp > max_obs {
                    max_obs = obs_exp;
                    max_bin = b;
                }
            }
            if max_obs >= 3.0 {
                if max_obs >= 5.0 {
                    status = Status::Fail;
                } else if status != Status::Fail {
                    status = Status::Warn;
                }
                rows.push(KmerRow {
                    sequence: decode_kmer(key),
                    count: total_for_key,
                    p_value: compute_pvalue(max_obs),
                    obs_exp: max_obs,
                    max_pos: bin_mid_percent(max_bin),
                });
            }
        }
        select_top(&mut rows);
        (rows, status)
    }
}

#[cfg(not(feature = "no-kmer"))]
pub use real::*;

#[cfg(feature = "no-kmer")]
mod stub {
    use crate::core::model::Status;

    pub const K: usize = 7;
    pub const BINS: usize = 10;
    pub const KMER_SPACE: usize = 1 << (2 * K);

    #[derive(Clone, Debug)]
    pub struct KmerRow {
        pub sequence: String,
        pub count: u64,
        pub p_value: f64,
        pub obs_exp: f64,
        pub max_pos: u32,
    }

    #[derive(Clone, Debug, Default)]
    pub struct BinCounts;

    impl BinCounts {
        pub fn new() -> Self {
            Self
        }
        pub fn add(&mut self, _key: u16) {}
        pub fn get(&self, _key: u16) -> u32 {
            0
        }
        pub fn merge(&mut self, _other: &BinCounts) {}
    }

    pub fn decode_kmer(_key: u16) -> String {
        String::new()
    }

    pub fn bin_mid_percent(_bin: usize) -> u32 {
        0
    }

    pub fn compute_pvalue(_obs_exp: f64) -> f64 {
        1.0
    }

    pub fn select_top(_rows: &mut Vec<KmerRow>) {}

    pub fn update_kmers(
        _seq: &[u8],
        _len: usize,
        _counts: &mut [BinCounts],
        _bin_counts: &mut [u64; BINS],
        _total: &mut u64,
        _timing: Option<&mut crate::core::metrics::UpdateTimings>,
    ) {
    }

    pub fn build_rows(
        _counts: &[BinCounts],
        _bin_counts: &[u64; BINS],
        _total: u64,
    ) -> (Vec<KmerRow>, Status) {
        (Vec::new(), Status::Pass)
    }
}

#[cfg(feature = "no-kmer")]
pub use stub::*;
