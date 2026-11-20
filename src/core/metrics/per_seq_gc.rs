use crate::core::model::Status;
use rustc_hash::FxHashMap;

pub struct PerSeqGcRow {
    pub gc: u8,
    pub count: f64,
}

pub(super) fn model_length(len: usize) -> usize {
    match len {
        0..=100 => len,
        101..=1000 => len / 100 * 100,
        _ => len / 1000 * 1000,
    }
}

#[derive(Clone, Debug, Default)]
pub(super) struct ShortGc(FxHashMap<(usize, usize), u64>);

impl ShortGc {
    pub fn add(&mut self, length: usize, gc: usize) {
        *self.0.entry((length, gc)).or_default() += 1;
    }

    pub fn merge(&mut self, other: &Self) {
        for (&key, &count) in &other.0 {
            *self.0.entry(key).or_default() += count;
        }
    }

    pub fn histogram(&self) -> [f64; 101] {
        let mut counts: Vec<_> = self.0.iter().collect();
        counts.sort_unstable_by_key(|&(key, _)| *key);
        let mut histogram = [0.0; 101];
        for (&(len, gc), &count) in counts {
            let len = len as i64;
            let gc = gc as i64;
            let low = ((100 * gc - 50).max(0) + len / 2) / len;
            let high = ((100 * gc + 50).min(100 * len) + len / 2) / len;
            for percent in low..=high {
                // Invert the rounded endpoints to count every GC count whose
                // interval includes this percentage, including shared edges.
                let first = (percent * len - 50 - len / 2 + 99).div_euclid(100).max(0);
                let last = ((percent + 1) * len - 1 + 50 - len / 2)
                    .div_euclid(100)
                    .min(len);
                histogram[percent as usize] += count as f64 / (last - first + 1) as f64;
            }
        }
        histogram
    }
}

pub(super) fn short_status(histogram: &[f64; 101]) -> Status {
    let total: f64 = histogram.iter().sum();
    if total <= 1.0 {
        return Status::Pass;
    }
    let mut peak = 0;
    for i in 1..histogram.len() {
        if histogram[i] > histogram[peak] {
            peak = i;
        }
    }
    let cutoff = histogram[peak] - histogram[peak] / 10.0;
    let start = (0..peak)
        .rev()
        .take_while(|&i| histogram[i] > cutoff)
        .last()
        .unwrap_or(peak);
    let end = (peak + 1..histogram.len())
        .take_while(|&i| histogram[i] > cutoff)
        .last()
        .unwrap_or(peak);
    let center = if start == 0 || end == histogram.len() - 1 {
        peak as f64
    } else {
        (start + end) as f64 / 2.0
    };
    let variance = histogram
        .iter()
        .enumerate()
        .map(|(i, &count)| (i as f64 - center).powi(2) * count)
        .sum::<f64>()
        / (total - 1.0);
    if variance == 0.0 {
        return Status::Pass;
    }
    let norm = total / (2.0 * std::f64::consts::PI * variance).sqrt();
    let deviation = histogram
        .iter()
        .enumerate()
        .map(|(i, &count)| {
            let expected = norm * (-(i as f64 - center).powi(2) / (2.0 * variance)).exp();
            (count - expected).abs()
        })
        .sum::<f64>()
        / total
        * 100.0;
    if deviation > 30.0 {
        Status::Fail
    } else if deviation > 15.0 {
        Status::Warn
    } else {
        Status::Pass
    }
}
