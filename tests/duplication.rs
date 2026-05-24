use kira_qc::core::metrics::duplication::{HEAP_COMPACT_THRESHOLD, SpaceSaving, hash_seq};
use std::collections::HashMap;

#[test]
fn add_records_first_occurrence() {
    let mut ss = SpaceSaving::new();
    ss.add(1, 1);
    ss.add(2, 1);
    ss.add(1, 1);

    let map: HashMap<u64, u64> = ss.entries().iter().map(|e| (e.key, e.count)).collect();
    assert_eq!(map.get(&1), Some(&2));
    assert_eq!(map.get(&2), Some(&1));
}

#[test]
fn case_fold_matches_uppercase_and_lowercase() {
    assert_eq!(hash_seq(b"acgt"), hash_seq(b"ACGT"));
    assert_eq!(hash_seq(b"AcGt"), hash_seq(b"ACGT"));
    assert_ne!(hash_seq(b"ACGT"), hash_seq(b"ACGA"));
}

#[test]
fn merge_combines_counts() {
    let mut a = SpaceSaving::new();
    let mut b = SpaceSaving::new();
    a.add(10, 3);
    a.add(20, 1);
    b.add(20, 2);
    b.add(30, 5);

    a.merge(&b);
    let counts: HashMap<u64, u64> = a.entries().iter().map(|e| (e.key, e.count)).collect();
    assert_eq!(counts.get(&10), Some(&3));
    assert_eq!(counts.get(&20), Some(&3));
    assert_eq!(counts.get(&30), Some(&5));
}

#[test]
fn heap_compacts_under_pressure() {
    let mut ss = SpaceSaving::new();
    for _ in 0..(HEAP_COMPACT_THRESHOLD + 100) {
        ss.add(42, 1);
    }
    assert!(ss.heap_len() <= HEAP_COMPACT_THRESHOLD);
    assert_eq!(ss.entries()[0].count, (HEAP_COMPACT_THRESHOLD + 100) as u64);
}
