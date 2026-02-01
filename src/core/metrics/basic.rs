pub struct BasicStats {
    pub file_type: &'static str,
    pub encoding: &'static str,
    pub total_sequences: u64,
    pub filtered_sequences: u64,
    pub min_len: u32,
    pub max_len: u32,
    pub gc_percent: u32,
}
