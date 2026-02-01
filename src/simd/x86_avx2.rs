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
        let q = if b >= offset { b - offset } else { 0 };
        sum += q as u64;
    }

    sum as u32
}

#[target_feature(enable = "avx2")]
pub unsafe fn prefix_scan_avx2(seq: &[u8], prefix: &[u8]) -> bool {
    if prefix.is_empty() || seq.len() < prefix.len() {
        return false;
    }
    let len = seq.len();
    let plen = prefix.len();
    let first = prefix[0] as i8;
    let upper_mask = _mm256_set1_epi8(0xDFu8 as i8);
    let target = _mm256_set1_epi8(first);
    let mut i = 0usize;
    while i + 32 <= len {
        let ptr = unsafe { seq.as_ptr().add(i) as *const __m256i };
        let mut v = unsafe { _mm256_loadu_si256(ptr) };
        v = _mm256_and_si256(v, upper_mask);
        let mask = _mm256_movemask_epi8(_mm256_cmpeq_epi8(v, target)) as u32;
        if mask != 0 {
            let mut m = mask;
            while m != 0 {
                let bit = m.trailing_zeros() as usize;
                let idx = i + bit;
                if idx + plen <= len {
                    let mut ok = true;
                    for j in 1..plen {
                        if (seq[idx + j] & 0xDF) != prefix[j] {
                            ok = false;
                            break;
                        }
                    }
                    if ok {
                        return true;
                    }
                }
                m &= m - 1;
            }
        }
        i += 32;
    }
    while i + plen <= len {
        if (seq[i] & 0xDF) == prefix[0] {
            let mut ok = true;
            for j in 1..plen {
                if (seq[i + j] & 0xDF) != prefix[j] {
                    ok = false;
                    break;
                }
            }
            if ok {
                return true;
            }
        }
        i += 1;
    }
    false
}

#[target_feature(enable = "avx2")]
pub unsafe fn encode_acgt_chunk_avx2(seq: &[u8], out: &mut [u8]) -> u32 {
    debug_assert!(seq.len() >= 32);
    debug_assert!(out.len() >= 32);
    let upper_mask = _mm256_set1_epi8(0xDFu8 as i8);
    let va = _mm256_set1_epi8(b'A' as i8);
    let vc = _mm256_set1_epi8(b'C' as i8);
    let vg = _mm256_set1_epi8(b'G' as i8);
    let vt = _mm256_set1_epi8(b'T' as i8);
    let c1 = _mm256_set1_epi8(1);
    let c2 = _mm256_set1_epi8(2);
    let c3 = _mm256_set1_epi8(3);
    let ptr = seq.as_ptr() as *const __m256i;
    let mut v = unsafe { _mm256_loadu_si256(ptr) };
    v = _mm256_and_si256(v, upper_mask);
    let is_a = _mm256_cmpeq_epi8(v, va);
    let is_c = _mm256_cmpeq_epi8(v, vc);
    let is_g = _mm256_cmpeq_epi8(v, vg);
    let is_t = _mm256_cmpeq_epi8(v, vt);
    let valid = _mm256_or_si256(_mm256_or_si256(is_a, is_c), _mm256_or_si256(is_g, is_t));
    let code_c = _mm256_and_si256(is_c, c1);
    let code_g = _mm256_and_si256(is_g, c2);
    let code_t = _mm256_and_si256(is_t, c3);
    let code = _mm256_or_si256(_mm256_or_si256(code_c, code_g), code_t);
    unsafe { _mm256_storeu_si256(out.as_mut_ptr() as *mut __m256i, code) };
    _mm256_movemask_epi8(valid) as u32
}

#[target_feature(enable = "avx2")]
pub unsafe fn acgt_2bit_encode_block_avx2(input: &[u8; 16]) -> (u16, [u8; 16]) {
    let upper_mask = _mm_set1_epi8(0xDFu8 as i8);
    let va = _mm_set1_epi8(b'A' as i8);
    let vc = _mm_set1_epi8(b'C' as i8);
    let vg = _mm_set1_epi8(b'G' as i8);
    let vt = _mm_set1_epi8(b'T' as i8);
    let c1 = _mm_set1_epi8(1);
    let c2 = _mm_set1_epi8(2);
    let c3 = _mm_set1_epi8(3);
    let v = unsafe { _mm_loadu_si128(input.as_ptr() as *const __m128i) };
    let v = _mm_and_si128(v, upper_mask);
    let is_a = _mm_cmpeq_epi8(v, va);
    let is_c = _mm_cmpeq_epi8(v, vc);
    let is_g = _mm_cmpeq_epi8(v, vg);
    let is_t = _mm_cmpeq_epi8(v, vt);
    let valid = _mm_or_si128(_mm_or_si128(is_a, is_c), _mm_or_si128(is_g, is_t));
    let code_c = _mm_and_si128(is_c, c1);
    let code_g = _mm_and_si128(is_g, c2);
    let code_t = _mm_and_si128(is_t, c3);
    let code = _mm_or_si128(_mm_or_si128(code_c, code_g), code_t);
    let mut out = [0u8; 16];
    unsafe { _mm_storeu_si128(out.as_mut_ptr() as *mut __m128i, code) };
    let mask = _mm_movemask_epi8(valid) as u16;
    (mask, out)
}

#[target_feature(enable = "avx2")]
pub unsafe fn acgt_2bit_block_16_avx2(input_ptr: *const u8) -> (u16, u32) {
    let upper_mask = _mm_set1_epi8(0xDFu8 as i8);
    let va = _mm_set1_epi8(b'A' as i8);
    let vc = _mm_set1_epi8(b'C' as i8);
    let vg = _mm_set1_epi8(b'G' as i8);
    let vt = _mm_set1_epi8(b'T' as i8);
    let c1 = _mm_set1_epi8(1);
    let c2 = _mm_set1_epi8(2);
    let c3 = _mm_set1_epi8(3);
    let v = unsafe { _mm_loadu_si128(input_ptr as *const __m128i) };
    let v = _mm_and_si128(v, upper_mask);
    let is_a = _mm_cmpeq_epi8(v, va);
    let is_c = _mm_cmpeq_epi8(v, vc);
    let is_g = _mm_cmpeq_epi8(v, vg);
    let is_t = _mm_cmpeq_epi8(v, vt);
    let valid = _mm_or_si128(_mm_or_si128(is_a, is_c), _mm_or_si128(is_g, is_t));
    let code_c = _mm_and_si128(is_c, c1);
    let code_g = _mm_and_si128(is_g, c2);
    let code_t = _mm_and_si128(is_t, c3);
    let code = _mm_or_si128(_mm_or_si128(code_c, code_g), code_t);
    let mut tmp = [0u8; 16];
    unsafe { _mm_storeu_si128(tmp.as_mut_ptr() as *mut __m128i, code) };
    let mut packed: u32 = 0;
    for i in 0..16 {
        packed |= (tmp[i] as u32) << (2 * i);
    }
    let mask = _mm_movemask_epi8(valid) as u16;
    (mask, packed)
}
