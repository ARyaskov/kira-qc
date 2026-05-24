use kira_qc::core::metrics::adapter_content::{ADAPTERS, scan, scan_any};

#[test]
fn cumulative_marks_from_start_to_end() {
    let mut seq = vec![b'A'; 30];
    seq.extend_from_slice(b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC");
    let mut counts = vec![[0u64; ADAPTERS.len()]; seq.len()];
    scan(&seq, &mut counts);

    for c in &counts[..30] {
        assert_eq!(c[0], 0);
    }
    for c in &counts[30..] {
        assert_eq!(c[0], 1);
    }
}

#[test]
fn first_occurrence_per_adapter_only() {
    let mut seq = Vec::new();
    seq.extend_from_slice(b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC");
    seq.extend_from_slice(&[b'N'; 5]);
    seq.extend_from_slice(b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC");

    let mut counts = vec![[0u64; ADAPTERS.len()]; seq.len()];
    scan(&seq, &mut counts);
    for c in &counts {
        assert_eq!(c[0], 1);
    }
}

#[test]
fn scan_any_sets_correct_flags() {
    let seq = b"NNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNN";
    let mut hits = [false; ADAPTERS.len()];
    scan_any(seq, &mut hits);
    assert!(hits[0]);
    for &hit in &hits[1..] {
        assert!(!hit);
    }
}
