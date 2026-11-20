#[derive(Clone, Debug)]
pub struct BasicStats {
    pub file_type: &'static str,
    pub encoding: &'static str,
    pub total_sequences: u64,
    pub total_bases: u64,
    /// FastQC 0.12 calls this "Sequences flagged as poor quality". kira-qc does
    /// no filtering of its own, so it is always 0.
    pub poor_quality_sequences: u64,
    pub min_len: u64,
    pub max_len: u64,
    pub gc_percent: u32,
}

impl BasicStats {
    /// FastQC 0.12 renders the base count with a magnitude suffix.
    pub fn format_bases(bases: u64) -> String {
        if bases < 1_000 {
            format!("{} bp", bases)
        } else if bases < 1_000_000 {
            format!("{} kbp", bases / 1_000)
        } else if bases < 1_000_000_000 {
            format!("{} Mbp", bases / 1_000_000)
        } else {
            format!("{} Gbp", bases / 1_000_000_000)
        }
    }

    pub fn length_range(&self) -> String {
        if self.min_len == self.max_len {
            self.min_len.to_string()
        } else {
            format!("{}-{}", self.min_len, self.max_len)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn base_counts_are_formatted_like_fastqc() {
        assert_eq!(BasicStats::format_bases(999), "999 bp");
        assert_eq!(BasicStats::format_bases(101_000), "101 kbp");
        assert_eq!(BasicStats::format_bases(5_400_000), "5 Mbp");
        assert_eq!(BasicStats::format_bases(3_200_000_000), "3 Gbp");
    }
}
