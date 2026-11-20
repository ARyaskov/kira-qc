#[derive(Clone, Debug)]
pub struct LengthDistRow {
    pub length: u64,
    pub count: u64,
}

/// Long-read length summary. N50/N90 are computed exactly from the full length
/// histogram rather than from log-bin midpoints.
#[derive(Clone, Debug)]
pub struct LongLengthSummary {
    pub bins: Vec<(String, u64)>,
    pub mean: f64,
    pub median: u64,
    pub n50: u64,
    pub n90: u64,
    pub min: u64,
    pub max: u64,
    pub total_bases: u64,
}
