#[cfg(target_arch = "aarch64")]
use std::arch::aarch64::*;

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
pub unsafe fn count_bases_neon(seq: &[u8]) -> (u32, u32, u32, u32, u32) {
    let mut a = 0u32;
    let mut c = 0u32;
    let mut g = 0u32;
    let mut t = 0u32;
    let mut n = 0u32;
    let mut i = 0usize;
    let len = seq.len();

    let upper_mask = vdupq_n_u8(0xDF);
    let va = vdupq_n_u8(b'A');
    let vc = vdupq_n_u8(b'C');
    let vg = vdupq_n_u8(b'G');
    let vt = vdupq_n_u8(b'T');
    let vn = vdupq_n_u8(b'N');
    let ones = vdupq_n_u8(1);

    // SAFETY: loop condition keeps `i+16 <= len`, so loads of 16 bytes stay in bounds.
    while i + 16 <= len {
        let ptr = unsafe { seq.as_ptr().add(i) };
        let v = unsafe { vld1q_u8(ptr) };
        let v = vandq_u8(v, upper_mask);

        let ma = vceqq_u8(v, va);
        let mc = vceqq_u8(v, vc);
        let mg = vceqq_u8(v, vg);
        let mt = vceqq_u8(v, vt);
        let mn = vceqq_u8(v, vn);

        a += vaddvq_u8(vandq_u8(ma, ones)) as u32;
        c += vaddvq_u8(vandq_u8(mc, ones)) as u32;
        g += vaddvq_u8(vandq_u8(mg, ones)) as u32;
        t += vaddvq_u8(vandq_u8(mt, ones)) as u32;
        n += vaddvq_u8(vandq_u8(mn, ones)) as u32;

        i += 16;
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

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
pub unsafe fn sum_qual_neon(qual: &[u8], offset: u8) -> u32 {
    let mut sum: u64 = 0;
    let mut i = 0usize;
    let len = qual.len();

    let off = vdupq_n_u8(offset);

    // SAFETY: loop condition keeps `i+16 <= len`, so loads of 16 bytes stay in bounds.
    while i + 16 <= len {
        let ptr = unsafe { qual.as_ptr().add(i) };
        let v = unsafe { vld1q_u8(ptr) };
        let q = vqsubq_u8(v, off);
        sum += vaddlvq_u8(q) as u64;
        i += 16;
    }

    for &b in &qual[i..] {
        sum += b.saturating_sub(offset) as u64;
    }

    sum as u32
}
