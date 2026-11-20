use crate::core::model::BaseGroup;

#[derive(Clone, Debug)]
pub struct PerBaseQualRow {
    pub group: BaseGroup,
    pub mean: f64,
    pub median: f64,
    pub lower_quartile: f64,
    pub upper_quartile: f64,
    pub p10: f64,
    pub p90: f64,
}

pub(super) fn position_percentiles(hist: &[u64], coverage: u64) -> [f64; 5] {
    let ranks = [10, 25, 50, 75, 90].map(|p| (coverage as u128 * p / 100) as u64);
    let mut values = [0.0; 5];
    let mut next = 0;
    let mut cumulative = 0;
    for (quality, count) in hist.iter().enumerate() {
        cumulative += count;
        while next < ranks.len() && cumulative >= ranks[next] {
            values[next] = quality as f64;
            next += 1;
        }
        if next == ranks.len() {
            break;
        }
    }
    values
}
