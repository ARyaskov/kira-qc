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
        sum += b.saturating_sub(offset) as u32;
    }
    sum
}
