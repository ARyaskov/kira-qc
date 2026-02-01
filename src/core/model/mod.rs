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
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Encoding {
    Sanger,
    Illumina15,
}

pub struct FinalizeContext {
    pub phred_offset: u8,
    pub encoding: Encoding,
    pub file_name: String,
    pub sample_name: String,
    pub mode: Mode,
}

pub const MAX_Q: usize = 93;

pub type QualHist = [u64; MAX_Q + 1];

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
