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

    while i + 16 <= len {
        let ptr = seq.as_ptr().add(i);
        let v = vld1q_u8(ptr);
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

    while i + 16 <= len {
        let ptr = qual.as_ptr().add(i);
        let v = vld1q_u8(ptr);
        let q = vqsubq_u8(v, off);
        sum += vaddlvq_u8(q) as u64;
        i += 16;
    }

    for &b in &qual[i..] {
        let q = if b >= offset { b - offset } else { 0 };
        sum += q as u64;
    }

    sum as u32
}

#[target_feature(enable = "neon")]
pub unsafe fn prefix_scan_neon(seq: &[u8], prefix: &[u8]) -> bool {
    if prefix.is_empty() || seq.len() < prefix.len() {
        return false;
    }
    let len = seq.len();
    let plen = prefix.len();
    let first = prefix[0];
    let upper_mask = vdupq_n_u8(0xDF);
    let target = vdupq_n_u8(first);
    let ones = vdupq_n_u8(1);
    let mut i = 0usize;
    while i + 16 <= len {
        let ptr = seq.as_ptr().add(i);
        let v = vandq_u8(vld1q_u8(ptr), upper_mask);
        let eq = vceqq_u8(v, target);
        let any = vaddvq_u8(vandq_u8(eq, ones));
        if any != 0 {
            for lane in 0..16 {
                let idx = i + lane;
                if idx + plen <= len {
                    if (seq[idx] & 0xDF) == prefix[0] {
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
                }
            }
        }
        i += 16;
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

#[target_feature(enable = "neon")]
pub unsafe fn encode_acgt_chunk_neon(seq: &[u8], out: &mut [u8]) -> u32 {
    debug_assert!(seq.len() >= 16);
    debug_assert!(out.len() >= 16);
    let upper_mask = vdupq_n_u8(0xDF);
    let va = vdupq_n_u8(b'A');
    let vc = vdupq_n_u8(b'C');
    let vg = vdupq_n_u8(b'G');
    let vt = vdupq_n_u8(b'T');
    let c1 = vdupq_n_u8(1);
    let c2 = vdupq_n_u8(2);
    let c3 = vdupq_n_u8(3);
    let v = vandq_u8(vld1q_u8(seq.as_ptr()), upper_mask);
    let is_a = vceqq_u8(v, va);
    let is_c = vceqq_u8(v, vc);
    let is_g = vceqq_u8(v, vg);
    let is_t = vceqq_u8(v, vt);
    let valid = vorrq_u8(vorrq_u8(is_a, is_c), vorrq_u8(is_g, is_t));
    let code_c = vandq_u8(is_c, c1);
    let code_g = vandq_u8(is_g, c2);
    let code_t = vandq_u8(is_t, c3);
    let code = vorrq_u8(vorrq_u8(code_c, code_g), code_t);
    vst1q_u8(out.as_mut_ptr(), code);
    let mut tmp = [0u8; 16];
    vst1q_u8(tmp.as_mut_ptr(), valid);
    let mut mask = 0u32;
    for i in 0..16 {
        if tmp[i] != 0 {
            mask |= 1u32 << i;
        }
    }
    mask
}

#[target_feature(enable = "neon")]
pub unsafe fn acgt_2bit_encode_block_neon(input: &[u8; 16]) -> (u16, [u8; 16]) {
    let upper_mask = vdupq_n_u8(0xDF);
    let va = vdupq_n_u8(b'A');
    let vc = vdupq_n_u8(b'C');
    let vg = vdupq_n_u8(b'G');
    let vt = vdupq_n_u8(b'T');
    let c1 = vdupq_n_u8(1);
    let c2 = vdupq_n_u8(2);
    let c3 = vdupq_n_u8(3);
    let v = vandq_u8(vld1q_u8(input.as_ptr()), upper_mask);
    let is_a = vceqq_u8(v, va);
    let is_c = vceqq_u8(v, vc);
    let is_g = vceqq_u8(v, vg);
    let is_t = vceqq_u8(v, vt);
    let valid = vorrq_u8(vorrq_u8(is_a, is_c), vorrq_u8(is_g, is_t));
    let code_c = vandq_u8(is_c, c1);
    let code_g = vandq_u8(is_g, c2);
    let code_t = vandq_u8(is_t, c3);
    let code = vorrq_u8(vorrq_u8(code_c, code_g), code_t);
    let mut out = [0u8; 16];
    vst1q_u8(out.as_mut_ptr(), code);
    let mut tmp = [0u8; 16];
    vst1q_u8(tmp.as_mut_ptr(), valid);
    let mut mask = 0u16;
    for i in 0..16 {
        if tmp[i] != 0 {
            mask |= 1u16 << i;
        }
    }
    (mask, out)
}

#[target_feature(enable = "neon")]
pub unsafe fn acgt_2bit_block_16_neon(input_ptr: *const u8) -> (u16, u32) {
    let upper_mask = vdupq_n_u8(0xDF);
    let va = vdupq_n_u8(b'A');
    let vc = vdupq_n_u8(b'C');
    let vg = vdupq_n_u8(b'G');
    let vt = vdupq_n_u8(b'T');
    let c1 = vdupq_n_u8(1);
    let c2 = vdupq_n_u8(2);
    let c3 = vdupq_n_u8(3);
    let v = vandq_u8(vld1q_u8(input_ptr), upper_mask);
    let is_a = vceqq_u8(v, va);
    let is_c = vceqq_u8(v, vc);
    let is_g = vceqq_u8(v, vg);
    let is_t = vceqq_u8(v, vt);
    let valid = vorrq_u8(vorrq_u8(is_a, is_c), vorrq_u8(is_g, is_t));
    let code_c = vandq_u8(is_c, c1);
    let code_g = vandq_u8(is_g, c2);
    let code_t = vandq_u8(is_t, c3);
    let code = vorrq_u8(vorrq_u8(code_c, code_g), code_t);
    let mut tmp = [0u8; 16];
    vst1q_u8(tmp.as_mut_ptr(), code);
    let mut packed: u32 = 0;
    for i in 0..16 {
        packed |= (tmp[i] as u32) << (2 * i);
    }
    let mut vtmp = [0u8; 16];
    vst1q_u8(vtmp.as_mut_ptr(), valid);
    let mut mask = 0u16;
    for i in 0..16 {
        if vtmp[i] != 0 {
            mask |= 1u16 << i;
        }
    }
    (mask, packed)
}
