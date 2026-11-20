//! Positional k-mer enrichment (K=7, 4^K = 16384 buckets per position bin).
//!
//! This module is a heuristic and is **off by default**, matching FastQC, which
//! stopped running its Kmer module by default in 0.11.x. Enable it with
//! `--kmers`.

#[cfg(not(feature = "no-kmer"))]
mod real {
    use crate::core::metrics::UpdateTimings;
    use std::time::Instant;

    pub const K: usize = 7;
    pub const BINS: usize = 10;
    pub const KMER_SPACE: usize = 1 << (2 * K);
    const MAX_REPORT: usize = 20;

    /// Report threshold on observed/expected ratio.
    const OBS_EXP_WARN: f64 = 3.0;
    const OBS_EXP_FAIL: f64 = 5.0;
    const SIGNIFICANCE: f64 = 0.05;

    /// Bonferroni-adjusted upper bound for enrichment under random placement
    /// of the observed occurrences across positions, conditional on their total.
    fn adjusted_enrichment_bound(obs: u64, trials: u64, frequency: f64) -> f64 {
        let observed = obs as f64 / trials as f64;
        if observed <= frequency {
            return 1.0;
        }
        let mut divergence = observed * (observed / frequency).ln();
        if observed < 1.0 {
            divergence += (1.0 - observed) * ((1.0 - observed) / (1.0 - frequency)).ln();
        }
        ((-(trials as f64) * divergence).exp() * (KMER_SPACE * BINS) as f64).min(1.0)
    }

    #[derive(Clone, Debug)]
    pub struct KmerRow {
        pub sequence: String,
        pub count: u64,
        pub p_value: f64,
        pub obs_exp: f64,
        /// Estimated base position of peak enrichment.
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

    /// Midpoint of a positional bin, as a percentage of read length.
    pub fn bin_mid_percent(bin: usize) -> u32 {
        if bin >= BINS - 1 {
            95
        } else {
            (bin as u32) * 10 + 5
        }
    }

    /// Bin midpoint expressed as a base position, given the observed read length.
    pub fn bin_mid_position(bin: usize, max_len: u64) -> u32 {
        if max_len == 0 {
            return 0;
        }
        let frac = (bin as f64 + 0.5) / BINS as f64;
        ((max_len as f64 * frac).round() as u32).max(1)
    }

    /// Complementary error function, Chebyshev approximation (Numerical Recipes
    /// `erfcc`); fractional error below 1.2e-7, which is far tighter than the
    /// two significant figures a p-value column shows.
    fn erfc(x: f64) -> f64 {
        let z = x.abs();
        let t = 2.0 / (2.0 + z);
        let ty = 4.0 * t - 2.0;

        // Published Chebyshev coefficients; kept verbatim from the reference
        // rather than rounded to f64 precision.
        #[allow(clippy::excessive_precision)]
        const COF: [f64; 28] = [
            -1.3026537197817094,
            6.419_697_923_564_902e-1,
            1.9476473204185836e-2,
            -9.561_514_786_808_631e-3,
            -9.46595344482036e-4,
            3.66839497852761e-4,
            4.2523324806907e-5,
            -2.0278578112534e-5,
            -1.624290004647e-6,
            1.303655835580e-6,
            1.5626441722e-8,
            -8.5238095915e-8,
            6.529054439e-9,
            5.059343495e-9,
            -9.91364156e-10,
            -2.27365122e-10,
            9.6467911e-11,
            2.394038e-12,
            -6.886027e-12,
            8.94487e-13,
            3.13092e-13,
            -1.12708e-13,
            3.81e-16,
            7.106e-15,
            -1.523e-15,
            -9.4e-17,
            1.21e-16,
            -2.8e-17,
        ];

        let mut d = 0.0f64;
        let mut dd = 0.0f64;
        for j in (1..COF.len()).rev() {
            let tmp = d;
            d = ty * d - dd + COF[j];
            dd = tmp;
        }
        let ans = t * (-z * z + 0.5 * (COF[0] + ty * d) - dd).exp();
        if x >= 0.0 { ans } else { 2.0 - ans }
    }

    /// Upper-tail probability of seeing at least `obs` successes in `n` trials
    /// with success probability `p`, via a normal approximation with continuity
    /// correction. Replaces the previous `10^(-2*(obs/exp - 1))` expression,
    /// which was not a probability of anything.
    pub fn binomial_upper_tail(obs: u64, n: u64, p: f64) -> f64 {
        if n == 0 || p <= 0.0 {
            return 1.0;
        }
        if p >= 1.0 {
            return 1.0;
        }
        let mean = n as f64 * p;
        let var = mean * (1.0 - p);
        if var <= 0.0 {
            return 1.0;
        }
        let z = (obs as f64 - mean - 0.5) / var.sqrt();
        if z <= 0.0 {
            return 1.0;
        }
        // Underflow to 0 is the conventional display for "astronomically
        // small"; clamping to the smallest subnormal would be misleading.
        let pv = 0.5 * erfc(z / std::f64::consts::SQRT_2);
        pv.clamp(0.0, 1.0)
    }

    pub fn select_top(rows: &mut Vec<KmerRow>) {
        rows.sort_by(|a, b| {
            b.obs_exp
                .partial_cmp(&a.obs_exp)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.sequence.cmp(&b.sequence))
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
        if len < K || counts.is_empty() {
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

    /// Rows for k-mers enriched in any positional bin.
    pub fn build_rows(
        counts: &[BinCounts],
        bin_counts: &[u64; BINS],
        total: u64,
        max_len: u64,
    ) -> (Vec<KmerRow>, crate::core::model::Status) {
        use crate::core::model::Status;
        let mut rows = Vec::new();
        let mut status = Status::Pass;
        if total == 0 || counts.is_empty() {
            return (rows, status);
        }

        for key in 0..(KMER_SPACE as u16) {
            let total_for_key: u64 = counts.iter().map(|c| c.get(key) as u64).sum();
            if total_for_key == 0 {
                continue;
            }
            let expected = total_for_key as f64 / total as f64;
            if expected <= 0.0 {
                continue;
            }
            let mut max_obs = 0.0f64;
            let mut max_bin = 0usize;
            let mut max_bin_obs = 0u64;
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
                    max_bin_obs = bc.get(key) as u64;
                }
            }
            if max_obs >= OBS_EXP_WARN {
                let p_value = adjusted_enrichment_bound(max_bin_obs, bin_counts[max_bin], expected);
                if p_value > SIGNIFICANCE {
                    continue;
                }
                if max_obs >= OBS_EXP_FAIL {
                    status = Status::Fail;
                } else if status != Status::Fail {
                    status = Status::Warn;
                }
                rows.push(KmerRow {
                    sequence: decode_kmer(key),
                    count: total_for_key,
                    p_value,
                    obs_exp: max_obs,
                    max_pos: bin_mid_position(max_bin, max_len),
                });
            }
        }
        select_top(&mut rows);
        (rows, status)
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn singleton_is_ignored_but_repeated_enrichment_fails() {
            use crate::core::model::Status;
            let mut bins: Vec<_> = (0..BINS).map(|_| BinCounts::new()).collect();
            bins[0].add(0);
            let (rows, status) = build_rows(&bins, &[1_000; BINS], 10_000, 100);
            assert!(rows.is_empty());
            assert_eq!(status, Status::Pass);
            for _ in 0..99 {
                bins[0].add(0);
            }
            let (rows, status) = build_rows(&bins, &[1_000; BINS], 10_000, 100);
            assert_eq!(status, Status::Fail);
            assert_eq!(rows[0].count, 100);
            assert!(rows[0].p_value < SIGNIFICANCE);
        }

        #[test]
        fn erfc_matches_known_values() {
            assert!((erfc(0.0) - 1.0).abs() < 1e-6);
            assert!((erfc(1.0) - 0.157_299_2).abs() < 1e-5);
            assert!((erfc(2.0) - 0.004_677_7).abs() < 1e-6);
        }

        #[test]
        fn p_value_is_a_probability() {
            // an observation right at the mean is not surprising
            assert!((binomial_upper_tail(100, 1000, 0.1) - 1.0).abs() < 1e-9);
            // a large excess is
            let p = binomial_upper_tail(300, 1000, 0.1);
            assert!(p < 1e-9, "p = {p}");
            // and every result stays inside [0, 1]
            for obs in [0u64, 1, 50, 500, 1000] {
                let p = binomial_upper_tail(obs, 1000, 0.1);
                assert!((0.0..=1.0).contains(&p), "obs={obs} p={p}");
            }
        }

        #[test]
        fn degenerate_inputs_do_not_panic() {
            assert_eq!(binomial_upper_tail(5, 0, 0.5), 1.0);
            assert_eq!(binomial_upper_tail(5, 10, 0.0), 1.0);
            assert_eq!(binomial_upper_tail(5, 10, 1.0), 1.0);
        }

        #[test]
        fn bin_position_is_in_read_coordinates() {
            assert_eq!(bin_mid_position(0, 100), 5);
            assert_eq!(bin_mid_position(9, 100), 95);
            assert_eq!(bin_mid_position(0, 0), 0);
        }
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

    pub fn bin_mid_position(_bin: usize, _max_len: u64) -> u32 {
        0
    }

    pub fn binomial_upper_tail(_obs: u64, _n: u64, _p: f64) -> f64 {
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
        _max_len: u64,
    ) -> (Vec<KmerRow>, Status) {
        (Vec::new(), Status::Pass)
    }
}

#[cfg(feature = "no-kmer")]
pub use stub::*;
