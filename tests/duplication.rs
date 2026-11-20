use kira_qc::core::metrics::duplication::{
    DUP_LEVEL_LABELS, SeqTracker, dup_level_index, duplication_status, hash_seq, truncate,
};
use kira_qc::core::model::Status;

#[test]
fn distinct_reads_are_counted_exactly() {
    let mut t = SeqTracker::new();
    t.add(b"ACGTACGTAC", 0);
    t.add(b"ACGTACGTAC", 1);
    t.add(b"TTTTTTTTTT", 2);

    let summary = t.duplication_levels();
    assert_eq!(summary.observed_reads, 3);
    assert_eq!(summary.distinct_reads, 2);
    assert!((summary.total_deduplicated_percent - (2.0 / 3.0 * 100.0)).abs() < 1e-9);
}

/// The previous Space-Saving sketch reported 100 % duplication for any library
/// with more distinct sequences than its capacity. This is the regression test
/// for that: a fully unique library must come out at 100 % deduplicated.
#[test]
fn high_complexity_library_is_not_reported_as_duplicated() {
    let mut t = SeqTracker::new();
    for i in 0..300_000u64 {
        let seq = format!("{:A>40}", i);
        t.add(seq.as_bytes(), i);
    }
    let summary = t.duplication_levels();
    assert_eq!(summary.distinct_reads, 300_000);
    assert!((summary.total_deduplicated_percent - 100.0).abs() < 1e-9);
    assert_eq!(
        duplication_status(summary.total_deduplicated_percent),
        Status::Pass
    );
}

#[test]
fn merging_is_order_independent() {
    let mut a = SeqTracker::new();
    let mut b = SeqTracker::new();
    for i in 0..500u64 {
        let seq = format!("{:A>30}", i % 120);
        if i % 2 == 0 {
            a.add(seq.as_bytes(), i);
        } else {
            b.add(seq.as_bytes(), i);
        }
    }

    let mut ab = a.clone();
    ab.merge(&b);
    let mut ba = b.clone();
    ba.merge(&a);

    assert_eq!(ab.entries_sorted(), ba.entries_sorted());
    assert_eq!(ab.observed_reads(), 500);
    assert_eq!(ab.distinct(), 120);
}

#[test]
fn merged_result_matches_single_pass() {
    let seqs: Vec<String> = (0..1000u64).map(|i| format!("{:A>30}", i % 70)).collect();

    let mut single = SeqTracker::new();
    for (i, s) in seqs.iter().enumerate() {
        single.add(s.as_bytes(), i as u64);
    }

    let mut parts = [SeqTracker::new(), SeqTracker::new(), SeqTracker::new()];
    for (i, s) in seqs.iter().enumerate() {
        parts[i % 3].add(s.as_bytes(), i as u64);
    }
    let mut merged = parts[0].clone();
    merged.merge(&parts[1]);
    merged.merge(&parts[2]);

    assert_eq!(single.entries_sorted(), merged.entries_sorted());
    assert_eq!(
        single.duplication_levels().total_deduplicated_percent,
        merged.duplication_levels().total_deduplicated_percent
    );
}

#[test]
fn duplication_levels_use_fastqc_buckets() {
    assert_eq!(DUP_LEVEL_LABELS.len(), 16);
    assert_eq!(DUP_LEVEL_LABELS[0], "1");
    assert_eq!(DUP_LEVEL_LABELS[9], ">10");
    assert_eq!(DUP_LEVEL_LABELS[15], ">10k+");
    assert_eq!(dup_level_index(1), 0);
    assert_eq!(dup_level_index(10), 9);
    assert_eq!(dup_level_index(100_000), 15);
}

#[test]
fn status_thresholds_match_fastqc() {
    assert_eq!(duplication_status(100.0), Status::Pass);
    assert_eq!(duplication_status(85.0), Status::Pass);
    assert_eq!(duplication_status(70.0), Status::Warn);
    assert_eq!(duplication_status(30.0), Status::Fail);
}

#[test]
fn long_sequences_are_truncated_like_fastqc() {
    let long = vec![b'A'; 300];
    assert_eq!(truncate(&long).len(), 50);
    let medium = vec![b'A'; 75];
    assert_eq!(truncate(&medium).len(), 75);
}

#[test]
fn hashing_folds_case() {
    assert_eq!(hash_seq(b"acgt"), hash_seq(b"ACGT"));
    assert_ne!(hash_seq(b"ACGT"), hash_seq(b"ACGA"));
}
