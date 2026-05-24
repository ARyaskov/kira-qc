#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
pub unsafe fn count_bases_avx2(seq: &[u8]) -> (u32, u32, u32, u32, u32) {
    let mut a = 0u32;
    let mut c = 0u32;
    let mut g = 0u32;
    let mut t = 0u32;
    let mut n = 0u32;
    let mut i = 0usize;
    let len = seq.len();

    let upper_mask = _mm256_set1_epi8(0xDFu8 as i8);
    let va = _mm256_set1_epi8(b'A' as i8);
    let vc = _mm256_set1_epi8(b'C' as i8);
    let vg = _mm256_set1_epi8(b'G' as i8);
    let vt = _mm256_set1_epi8(b'T' as i8);
    let vn = _mm256_set1_epi8(b'N' as i8);

    // SAFETY: loop condition keeps `i+32 <= len`, so loads of 32 bytes stay in bounds.
    while i + 32 <= len {
        let ptr = unsafe { seq.as_ptr().add(i) as *const __m256i };
        let mut v = unsafe { _mm256_loadu_si256(ptr) };
        v = _mm256_and_si256(v, upper_mask);
        let ma = _mm256_movemask_epi8(_mm256_cmpeq_epi8(v, va)) as u32;
        let mc = _mm256_movemask_epi8(_mm256_cmpeq_epi8(v, vc)) as u32;
        let mg = _mm256_movemask_epi8(_mm256_cmpeq_epi8(v, vg)) as u32;
        let mt = _mm256_movemask_epi8(_mm256_cmpeq_epi8(v, vt)) as u32;
        let mn = _mm256_movemask_epi8(_mm256_cmpeq_epi8(v, vn)) as u32;
        a += ma.count_ones();
        c += mc.count_ones();
        g += mg.count_ones();
        t += mt.count_ones();
        n += mn.count_ones();
        i += 32;
    }

    for &b in &seq[i..] {
        match b & 0xDF {
            b'A' => a += 1,
            b'C' => c += 1,
            b'G' => g += 1,
            b'T' => t += 1,
            b'N' => n += 1,
            _ => {}
        }
    }

    (a, c, g, t, n)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
pub unsafe fn sum_qual_avx2(qual: &[u8], offset: u8) -> u32 {
    let mut sum: u64 = 0;
    let mut i = 0usize;
    let len = qual.len();

    let off = _mm256_set1_epi8(offset as i8);
    let zero = _mm256_setzero_si256();

    // SAFETY: loop condition keeps `i+32 <= len`; `tmp` is a 32-byte local for the SAD store.
    while i + 32 <= len {
        let ptr = unsafe { qual.as_ptr().add(i) as *const __m256i };
        let v = unsafe { _mm256_loadu_si256(ptr) };
        let q = _mm256_subs_epu8(v, off);
        let sad = _mm256_sad_epu8(q, zero);
        let mut tmp = [0u64; 4];
        unsafe { _mm256_storeu_si256(tmp.as_mut_ptr() as *mut __m256i, sad) };
        sum += tmp[0] + tmp[1] + tmp[2] + tmp[3];
        i += 32;
    }

    for &b in &qual[i..] {
        sum += b.saturating_sub(offset) as u64;
    }

    sum as u32
}
