#[cfg(target_arch = "aarch64")]
mod aarch64_neon;
mod scalar;
#[cfg(target_arch = "x86_64")]
mod x86_avx2;
#[cfg(target_arch = "x86_64")]
use std::sync::OnceLock;

#[cfg(target_arch = "x86_64")]
#[inline]
fn avx2_available() -> bool {
    static AVX2: OnceLock<bool> = OnceLock::new();
    *AVX2.get_or_init(|| std::arch::is_x86_feature_detected!("avx2"))
}

/// (A, C, G, T, N) byte tally; case-folded via `& 0xDF`.
#[inline]
pub fn count_bases(seq: &[u8]) -> (u32, u32, u32, u32, u32) {
    #[cfg(target_arch = "x86_64")]
    {
        if avx2_available() {
            // SAFETY: guarded by runtime AVX2 feature detection.
            unsafe {
                return x86_avx2::count_bases_avx2(seq);
            }
        }
    }
    #[cfg(target_arch = "aarch64")]
    {
        // SAFETY: NEON is part of the aarch64 baseline ABI.
        unsafe {
            return aarch64_neon::count_bases_neon(seq);
        }
    }
    #[cfg_attr(
        any(target_arch = "x86_64", target_arch = "aarch64"),
        allow(unreachable_code)
    )]
    scalar::count_bases(seq)
}

/// Sum of `max(0, qual[i] - offset)`.
#[inline]
pub fn sum_qual(qual: &[u8], offset: u8) -> u32 {
    #[cfg(target_arch = "x86_64")]
    {
        if avx2_available() {
            // SAFETY: guarded by runtime AVX2 feature detection.
            unsafe {
                return x86_avx2::sum_qual_avx2(qual, offset);
            }
        }
    }
    #[cfg(target_arch = "aarch64")]
    {
        // SAFETY: NEON is part of the aarch64 baseline ABI.
        unsafe {
            return aarch64_neon::sum_qual_neon(qual, offset);
        }
    }
    #[cfg_attr(
        any(target_arch = "x86_64", target_arch = "aarch64"),
        allow(unreachable_code)
    )]
    scalar::sum_qual(qual, offset)
}
