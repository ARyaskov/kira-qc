#![cfg(not(feature = "no-kmer"))]

use kira_qc::core::metrics::kmer_content::{
    BINS, BinCounts, K, KMER_SPACE, decode_kmer, update_kmers,
};

fn key_of(seq: &[u8; K]) -> u16 {
    let mut v: u64 = 0;
    for &b in seq {
        let bits: u64 = match b & 0xDF {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => panic!("test k-mer must use ACGT only"),
        };
        v = (v << 2) | bits;
    }
    v as u16
}

#[test]
fn encode_decode_round_trip_all_kmers() {
    for k in 0..(KMER_SPACE as u16) {
        let s = decode_kmer(k);
        assert_eq!(s.len(), K);
        let bytes = s.as_bytes();
        let mut buf = [b'A'; K];
        buf.copy_from_slice(bytes);
        assert_eq!(key_of(&buf), k, "round trip failed for key {k}");
    }
}

#[test]
fn update_kmers_counts_clean_read_exactly() {
    let seq = b"ACGTACGTACGTACGTACGT";
    let mut bins = (0..BINS).map(|_| BinCounts::new()).collect::<Vec<_>>();
    let mut bin_counts = [0u64; BINS];
    let mut total = 0u64;

    update_kmers(seq, seq.len(), &mut bins, &mut bin_counts, &mut total, None);
    assert_eq!(total, 14);
    assert_eq!(bin_counts.iter().sum::<u64>(), 14);

    let key = key_of(b"ACGTACG");
    let occurrences: u32 = bins.iter().map(|b| b.get(key)).sum();
    assert!(occurrences >= 1);
}

#[test]
fn update_kmers_skips_over_n() {
    let seq = b"ACGTACGNACGTACGT";
    let mut bins = (0..BINS).map(|_| BinCounts::new()).collect::<Vec<_>>();
    let mut bin_counts = [0u64; BINS];
    let mut total = 0u64;

    update_kmers(seq, seq.len(), &mut bins, &mut bin_counts, &mut total, None);
    assert_eq!(total, 3);

    let acgtacg = key_of(b"ACGTACG");
    assert_eq!(bins.iter().map(|b| b.get(acgtacg)).sum::<u32>(), 2);
    let cgtacgt = key_of(b"CGTACGT");
    assert_eq!(bins.iter().map(|b| b.get(cgtacgt)).sum::<u32>(), 1);
}

#[test]
fn update_kmers_skips_short_reads() {
    let seq = b"ACGTAC";
    let mut bins = (0..BINS).map(|_| BinCounts::new()).collect::<Vec<_>>();
    let mut bin_counts = [0u64; BINS];
    let mut total = 0u64;

    update_kmers(seq, seq.len(), &mut bins, &mut bin_counts, &mut total, None);
    assert_eq!(total, 0);
}

#[test]
fn bin_counts_merge_is_associative() {
    let mut a = BinCounts::new();
    let mut b = BinCounts::new();
    a.add(42);
    a.add(42);
    b.add(42);
    b.add(7);

    let mut merged = a.clone();
    merged.merge(&b);
    assert_eq!(merged.get(42), 3);
    assert_eq!(merged.get(7), 1);

    let mut other = b.clone();
    other.merge(&a);
    assert_eq!(other.get(42), 3);
    assert_eq!(other.get(7), 1);
}
