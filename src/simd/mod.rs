#[cfg(target_arch = "aarch64")]
mod aarch64_neon;
mod scalar;
#[cfg(target_arch = "x86_64")]
mod x86_avx2;

#[cfg(target_arch = "x86_64")]
pub const KMER_CHUNK: usize = 32;
#[cfg(target_arch = "aarch64")]
pub const KMER_CHUNK: usize = 16;
#[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
pub const KMER_CHUNK: usize = 16;

pub fn count_bases(seq: &[u8]) -> (u32, u32, u32, u32, u32) {
    #[cfg(target_arch = "x86_64")]
    unsafe {
        return x86_avx2::count_bases_avx2(seq);
    }
    #[cfg(target_arch = "aarch64")]
    unsafe {
        return aarch64_neon::count_bases_neon(seq);
    }
    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    {
        return scalar::count_bases(seq);
    }
}

pub fn sum_qual(qual: &[u8], offset: u8) -> u32 {
    #[cfg(target_arch = "x86_64")]
    unsafe {
        return x86_avx2::sum_qual_avx2(qual, offset);
    }
    #[cfg(target_arch = "aarch64")]
    unsafe {
        return aarch64_neon::sum_qual_neon(qual, offset);
    }
    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    {
        return scalar::sum_qual(qual, offset);
    }
}

pub fn prefix_scan(seq: &[u8], prefix: &[u8]) -> bool {
    if prefix.is_empty() || seq.len() < prefix.len() {
        return false;
    }
    #[cfg(target_arch = "x86_64")]
    unsafe {
        return x86_avx2::prefix_scan_avx2(seq, prefix);
    }
    #[cfg(target_arch = "aarch64")]
    unsafe {
        return aarch64_neon::prefix_scan_neon(seq, prefix);
    }
    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    {
        return scalar::prefix_scan(seq, prefix);
    }
}

pub fn encode_acgt_chunk(seq: &[u8], out: &mut [u8]) -> u32 {
    #[cfg(target_arch = "x86_64")]
    unsafe {
        return x86_avx2::encode_acgt_chunk_avx2(seq, out);
    }
    #[cfg(target_arch = "aarch64")]
    unsafe {
        return aarch64_neon::encode_acgt_chunk_neon(seq, out);
    }
    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    {
        return scalar::encode_acgt_chunk_scalar(seq, out);
    }
}

pub fn acgt_2bit_encode_block(input: &[u8; 16]) -> (u16, [u8; 16]) {
    #[cfg(target_arch = "x86_64")]
    unsafe {
        return x86_avx2::acgt_2bit_encode_block_avx2(input);
    }
    #[cfg(target_arch = "aarch64")]
    unsafe {
        return aarch64_neon::acgt_2bit_encode_block_neon(input);
    }
    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    {
        return scalar::acgt_2bit_encode_block_scalar(input);
    }
}

pub fn acgt_2bit_block_16(input_ptr: *const u8) -> (u16, u32) {
    #[cfg(target_arch = "x86_64")]
    unsafe {
        return x86_avx2::acgt_2bit_block_16_avx2(input_ptr);
    }
    #[cfg(target_arch = "aarch64")]
    unsafe {
        return aarch64_neon::acgt_2bit_block_16_neon(input_ptr);
    }
    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    {
        return scalar::acgt_2bit_block_16_scalar(input_ptr);
    }
}
