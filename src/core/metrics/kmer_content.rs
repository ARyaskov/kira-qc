#[cfg(not(feature = "no-kmer"))]
mod real {
    use crate::core::metrics::UpdateTimings;
    use crate::simd;
    use std::cmp::Reverse;
    use std::collections::{BinaryHeap, HashMap};
    use std::time::Instant;

    pub const K: usize = 7;
    pub const BINS: usize = 10;
    const CMS_DEPTH: usize = 4;
    const CMS_WIDTH: usize = 1 << 18;
    // Per-bin heavy hitters; sized for stability while keeping memory bounded.
    const HH_K: usize = 2000;
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
    pub struct Cms {
        data: Vec<u32>,
    }

    impl Cms {
        pub fn new() -> Self {
            Self {
                data: vec![0u32; CMS_DEPTH * CMS_WIDTH],
            }
        }

        pub fn add(&mut self, key: u64, weight: u32) {
            for d in 0..CMS_DEPTH {
                let idx = self.index(key, d);
                let slot = &mut self.data[d * CMS_WIDTH + idx];
                *slot = slot.saturating_add(weight);
            }
        }

        pub fn estimate(&self, key: u64) -> u32 {
            let mut min = u32::MAX;
            for d in 0..CMS_DEPTH {
                let idx = self.index(key, d);
                let v = self.data[d * CMS_WIDTH + idx];
                if v < min {
                    min = v;
                }
            }
            if min == u32::MAX { 0 } else { min }
        }

        pub fn merge(&mut self, other: &Cms) {
            for i in 0..self.data.len() {
                self.data[i] = self.data[i].saturating_add(other.data[i]);
            }
        }

        #[inline]
        fn index(&self, key: u64, depth: usize) -> usize {
            let mut x = key ^ ((depth as u64).wrapping_mul(0x9e3779b97f4a7c15));
            x ^= x >> 33;
            x = x.wrapping_mul(0xff51afd7ed558ccd);
            x ^= x >> 33;
            (x as usize) & (CMS_WIDTH - 1)
        }
    }

    #[derive(Clone, Debug)]
    struct Entry {
        key: u64,
        count: u64,
    }

    #[derive(Clone, Debug)]
    pub struct SpaceSaving {
        map: HashMap<u64, usize>,
        entries: Vec<Entry>,
        heap: BinaryHeap<(Reverse<u64>, u64, usize)>,
    }

    impl SpaceSaving {
        pub fn new() -> Self {
            Self {
                map: HashMap::with_capacity(HH_K),
                entries: Vec::with_capacity(HH_K),
                heap: BinaryHeap::with_capacity(HH_K),
            }
        }

        pub fn add(&mut self, key: u64, weight: u64) {
            if let Some(&idx) = self.map.get(&key) {
                let e = &mut self.entries[idx];
                e.count += weight;
                self.heap.push((Reverse(e.count), e.key, idx));
                return;
            }

            if self.entries.len() < HH_K {
                let idx = self.entries.len();
                self.entries.push(Entry { key, count: weight });
                self.map.insert(key, idx);
                self.heap.push((Reverse(weight), key, idx));
                return;
            }

            let (min_idx, min_count) = self.min_entry();
            let removed = self.entries[min_idx].key;
            self.map.remove(&removed);
            self.entries[min_idx] = Entry {
                key,
                count: min_count + weight,
            };
            self.map.insert(key, min_idx);
            self.heap.push((Reverse(min_count + weight), key, min_idx));
        }

        pub fn merge(&mut self, other: &SpaceSaving) {
            let mut items = other.entries.clone();
            items.sort_by_key(|e| e.key);
            for e in items {
                self.add(e.key, e.count);
            }
        }

        pub fn keys(&self) -> Vec<u64> {
            self.entries.iter().map(|e| e.key).collect()
        }

        fn min_entry(&mut self) -> (usize, u64) {
            loop {
                if let Some((Reverse(count), key, idx)) = self.heap.pop() {
                    let e = &self.entries[idx];
                    if e.key == key && e.count == count {
                        return (idx, count);
                    }
                } else {
                    return (0, self.entries[0].count);
                }
            }
        }
    }

    pub fn encode_kmer(seq: &[u8]) -> Option<u64> {
        if seq.len() != K {
            return None;
        }
        let mut v = 0u64;
        for &b in seq {
            let u = b & 0xDF;
            let bits = match u {
                b'A' => 0u64,
                b'C' => 1u64,
                b'G' => 2u64,
                b'T' => 3u64,
                _ => return None,
            };
            v = (v << 2) | bits;
        }
        Some(v)
    }

    pub fn decode_kmer(mut key: u64) -> String {
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
        String::from_utf8_lossy(&buf).to_string()
    }

    pub fn pos_bin(pos: usize, len: usize) -> usize {
        if len == 0 {
            return 0;
        }
        let pct = ((pos + 1) * 100) / len;
        let bin = pct / 10;
        if bin >= BINS { BINS - 1 } else { bin }
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
        ((t * len) + 99) / 100
    }

    pub fn update_kmers(
        seq: &[u8],
        len: usize,
        cms: &mut [Cms],
        hh: &mut [SpaceSaving],
        bin_counts: &mut [u64; BINS],
        total: &mut u64,
        mut timing: Option<&mut UpdateTimings>,
    ) {
        if len < K {
            return;
        }
        let t_total: Option<Instant> = timing.as_deref_mut().map(|_| Instant::now());
        let mask: u64 = (1u64 << (2 * K)) - 1;
        let mut pos: usize = 0;
        let mut bin = 0usize;
        let mut next_threshold = next_bin_threshold(len, bin);
        let mut carry_bits: u64 = 0;
        let mut carry_len: usize = 0;

        const BATCH: usize = 256;
        let mut batch_keys: [u16; BATCH] = [0u16; BATCH];
        let mut batch_bins: [u8; BATCH] = [0u8; BATCH];
        let mut batch_len: usize = 0;

        #[inline(always)]
        fn flush_batch(
            batch_len: &mut usize,
            keys: &mut [u16; BATCH],
            bins: &mut [u8; BATCH],
            cms: &mut [Cms],
            hh: &mut [SpaceSaving],
            bin_counts: &mut [u64; BINS],
            total: &mut u64,
            timing: Option<&mut UpdateTimings>,
        ) {
            let len = *batch_len;
            if len == 0 {
                return;
            }
            if let Some(t) = timing {
                for i in 0..len {
                    let bin = bins[i] as usize;
                    let key = keys[i] as u64;
                    let t0 = Instant::now();
                    cms[bin].add(key, 1);
                    t.kmer_cms += t0.elapsed();
                    let t1 = Instant::now();
                    hh[bin].add(key, 1);
                    t.kmer_hh += t1.elapsed();
                    bin_counts[bin] += 1;
                    *total += 1;
                }
            } else {
                for i in 0..len {
                    let bin = bins[i] as usize;
                    let key = keys[i] as u64;
                    cms[bin].add(key, 1);
                    hh[bin].add(key, 1);
                    bin_counts[bin] += 1;
                    *total += 1;
                }
            }
            *batch_len = 0;
        }

        while pos + 16 <= len {
            let t_encode: Option<Instant> = timing.as_deref_mut().map(|_| Instant::now());
            let (valid_mask, packed_codes) =
                simd::acgt_2bit_block_16(unsafe { seq.as_ptr().add(pos) });
            if let (Some(t), Some(t0)) = (timing.as_deref_mut(), t_encode) {
                t.kmer_encode += t0.elapsed();
            }

            let combined_len = carry_len + 16;
            let mut vbits: u32 = if carry_len == 0 {
                valid_mask as u32
            } else {
                ((valid_mask as u32) << carry_len) | ((1u32 << carry_len) - 1)
            };
            if combined_len < 32 {
                vbits &= (1u32 << combined_len) - 1;
            }

            let mut w = vbits;
            w &= w >> 1;
            w &= w >> 2;
            w &= w >> 3;
            w &= w >> 4;
            w &= w >> 5;
            w &= w >> 6;
            let max_start = combined_len - K;
            if max_start < 31 {
                w &= (1u32 << (max_start + 1)) - 1;
            }

            let stream_bits: u64 = carry_bits | ((packed_codes as u64) << (2 * carry_len));
            let base_start = pos - carry_len;

            if let Some(t) = timing.as_deref_mut() {
                let t0 = Instant::now();
                let mut m = w;
                while m != 0 {
                    let i = m.trailing_zeros() as usize;
                    let start_pos = base_start + i;
                    let start_pos_plus1 = start_pos + 1;
                    let t_bin = Instant::now();
                    while bin + 1 < BINS && start_pos_plus1 >= next_threshold {
                        bin += 1;
                        next_threshold = next_bin_threshold(len, bin);
                    }
                    t.kmer_binning += t_bin.elapsed();
                    let key = ((stream_bits >> (2 * i)) & mask) as u16;
                    batch_keys[batch_len] = key;
                    batch_bins[batch_len] = bin as u8;
                    batch_len += 1;
                    t.kmer_updates += 1;
                    if batch_len == BATCH {
                        flush_batch(
                            &mut batch_len,
                            &mut batch_keys,
                            &mut batch_bins,
                            cms,
                            hh,
                            bin_counts,
                            total,
                            Some(t),
                        );
                    }
                    m &= m - 1;
                }
                t.kmer_keygen += t0.elapsed();
            } else {
                let mut m = w;
                while m != 0 {
                    let i = m.trailing_zeros() as usize;
                    let start_pos = base_start + i;
                    let start_pos_plus1 = start_pos + 1;
                    while bin + 1 < BINS && start_pos_plus1 >= next_threshold {
                        bin += 1;
                        next_threshold = next_bin_threshold(len, bin);
                    }
                    let key = ((stream_bits >> (2 * i)) & mask) as u16;
                    batch_keys[batch_len] = key;
                    batch_bins[batch_len] = bin as u8;
                    batch_len += 1;
                    if batch_len == BATCH {
                        flush_batch(
                            &mut batch_len,
                            &mut batch_keys,
                            &mut batch_bins,
                            cms,
                            hh,
                            bin_counts,
                            total,
                            None,
                        );
                    }
                    m &= m - 1;
                }
            }

            let mut suffix = 0usize;
            for s in 0..combined_len.min(6) {
                let idx = combined_len - 1 - s;
                if ((vbits >> idx) & 1) != 0 {
                    suffix += 1;
                } else {
                    break;
                }
            }
            carry_len = suffix;
            if carry_len > 0 {
                let shift = 2 * (combined_len - carry_len);
                carry_bits = (stream_bits >> shift) & ((1u64 << (2 * carry_len)) - 1);
            } else {
                carry_bits = 0;
            }
            pos += 16;
        }

        let mut rolling = carry_bits;
        let mut valid_run = carry_len;
        let t_tail: Option<Instant> = timing.as_deref_mut().map(|_| Instant::now());
        while pos < len {
            let b = seq[pos] & 0xDF;
            let bits = match b {
                b'A' => 0u64,
                b'C' => 1u64,
                b'G' => 2u64,
                b'T' => 3u64,
                _ => {
                    valid_run = 0;
                    rolling = 0;
                    pos += 1;
                    continue;
                }
            };
            if valid_run < K {
                valid_run += 1;
            }
            rolling = ((rolling << 2) | bits) & mask;
            if valid_run >= K {
                let start_pos_plus1 = pos + 2 - K;
                if let Some(t) = timing.as_deref_mut() {
                    let t_bin = Instant::now();
                    while bin + 1 < BINS && start_pos_plus1 >= next_threshold {
                        bin += 1;
                        next_threshold = next_bin_threshold(len, bin);
                    }
                    t.kmer_binning += t_bin.elapsed();
                    t.kmer_updates += 1;
                } else {
                    while bin + 1 < BINS && start_pos_plus1 >= next_threshold {
                        bin += 1;
                        next_threshold = next_bin_threshold(len, bin);
                    }
                }
                batch_keys[batch_len] = rolling as u16;
                batch_bins[batch_len] = bin as u8;
                batch_len += 1;
                if batch_len == BATCH {
                    flush_batch(
                        &mut batch_len,
                        &mut batch_keys,
                        &mut batch_bins,
                        cms,
                        hh,
                        bin_counts,
                        total,
                        timing.as_deref_mut(),
                    );
                }
            }
            pos += 1;
        }
        if let (Some(t), Some(t0)) = (timing.as_deref_mut(), t_tail) {
            t.kmer_keygen += t0.elapsed();
        }
        flush_batch(
            &mut batch_len,
            &mut batch_keys,
            &mut batch_bins,
            cms,
            hh,
            bin_counts,
            total,
            timing.as_deref_mut(),
        );

        if let (Some(t), Some(t0)) = (timing.as_deref_mut(), t_total) {
            t.kmer += t0.elapsed();
        }
    }
}

#[cfg(not(feature = "no-kmer"))]
pub use real::*;

#[cfg(feature = "no-kmer")]
mod stub {
    pub const K: usize = 7;
    pub const BINS: usize = 10;

    #[derive(Clone, Debug)]
    pub struct KmerRow {
        pub sequence: String,
        pub count: u64,
        pub p_value: f64,
        pub obs_exp: f64,
        pub max_pos: u32,
    }

    #[derive(Clone, Debug)]
    pub struct Cms;

    impl Cms {
        pub fn new() -> Self {
            Self
        }
        pub fn add(&mut self, _key: u64, _weight: u32) {}
        pub fn estimate(&self, _key: u64) -> u32 {
            0
        }
        pub fn merge(&mut self, _other: &Cms) {}
    }

    #[derive(Clone, Debug)]
    pub struct SpaceSaving;

    impl SpaceSaving {
        pub fn new() -> Self {
            Self
        }
        pub fn add(&mut self, _key: u64, _weight: u64) {}
        pub fn merge(&mut self, _other: &SpaceSaving) {}
        pub fn keys(&self) -> Vec<u64> {
            Vec::new()
        }
    }

    pub fn decode_kmer(_key: u64) -> String {
        String::new()
    }

    pub fn pos_bin(_pos: usize, _len: usize) -> usize {
        0
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
        _cms: &mut [Cms],
        _hh: &mut [SpaceSaving],
        _bin_counts: &mut [u64; BINS],
        _total: &mut u64,
        _timing: Option<&mut crate::core::metrics::UpdateTimings>,
    ) {
    }
}

#[cfg(feature = "no-kmer")]
pub use stub::*;
