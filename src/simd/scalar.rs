pub fn count_bases(seq: &[u8]) -> (u32, u32, u32, u32, u32) {
    let mut a = 0u32;
    let mut c = 0u32;
    let mut g = 0u32;
    let mut t = 0u32;
    let mut n = 0u32;
    for &b in seq {
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

pub fn sum_qual(qual: &[u8], offset: u8) -> u32 {
    let mut sum: u32 = 0;
    for &b in qual {
        let q = if b >= offset { b - offset } else { 0 };
        sum += q as u32;
    }
    sum
}

pub fn prefix_scan(seq: &[u8], prefix: &[u8]) -> bool {
    if prefix.is_empty() || seq.len() < prefix.len() {
        return false;
    }
    let first = prefix[0];
    let len = seq.len();
    let plen = prefix.len();
    let mut i = 0usize;
    while i + plen <= len {
        if (seq[i] & 0xDF) == first {
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

pub fn encode_acgt_chunk_scalar(seq: &[u8], out: &mut [u8]) -> u32 {
    let mut mask = 0u32;
    let n = out.len().min(seq.len());
    for i in 0..n {
        let b = seq[i] & 0xDF;
        match b {
            b'A' => {
                out[i] = 0;
                mask |= 1u32 << i;
            }
            b'C' => {
                out[i] = 1;
                mask |= 1u32 << i;
            }
            b'G' => {
                out[i] = 2;
                mask |= 1u32 << i;
            }
            b'T' => {
                out[i] = 3;
                mask |= 1u32 << i;
            }
            _ => {
                out[i] = 0;
            }
        }
    }
    mask
}

pub fn acgt_2bit_encode_block_scalar(input: &[u8; 16]) -> (u16, [u8; 16]) {
    let mut out = [0u8; 16];
    let mut mask = 0u16;
    for i in 0..16 {
        let b = input[i] & 0xDF;
        match b {
            b'A' => {
                out[i] = 0;
                mask |= 1u16 << i;
            }
            b'C' => {
                out[i] = 1;
                mask |= 1u16 << i;
            }
            b'G' => {
                out[i] = 2;
                mask |= 1u16 << i;
            }
            b'T' => {
                out[i] = 3;
                mask |= 1u16 << i;
            }
            _ => {
                out[i] = 0;
            }
        }
    }
    (mask, out)
}

pub fn acgt_2bit_block_16_scalar(input_ptr: *const u8) -> (u16, u32) {
    let mut packed: u32 = 0;
    let mut mask: u16 = 0;
    for i in 0..16 {
        let b = unsafe { *input_ptr.add(i) } & 0xDF;
        let code = match b {
            b'A' => 0u32,
            b'C' => 1u32,
            b'G' => 2u32,
            b'T' => 3u32,
            _ => {
                continue;
            }
        };
        packed |= code << (2 * i);
        mask |= 1u16 << i;
    }
    (mask, packed)
}
