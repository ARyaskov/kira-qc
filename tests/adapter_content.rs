use kira_qc::core::metrics::adapter_content::{
    AdapterSet, DEFAULT_LONG_ADAPTERS, DEFAULT_SHORT_ADAPTERS, adapter_status, revcomp,
};
use kira_qc::core::model::{Mode, Status};

#[test]
fn default_list_matches_fastqc_012() {
    // FastQC 0.12 dropped SOLiD and added PolyA / PolyG, and matches 12 bp
    // prefixes rather than whole adapters.
    let names: Vec<&str> = DEFAULT_SHORT_ADAPTERS.iter().map(|(n, _)| *n).collect();
    assert!(names.contains(&"Illumina Universal Adapter"));
    assert!(names.contains(&"Nextera Transposase Sequence"));
    assert!(names.contains(&"PolyA"));
    assert!(names.contains(&"PolyG"));
    assert!(!names.iter().any(|n| n.contains("SOLiD")));
    for (name, probe) in DEFAULT_SHORT_ADAPTERS {
        assert_eq!(probe.len(), 12, "{name} probe should be a 12-mer");
    }
}

/// Regression test for the central adapter bug: real read-through leaves a
/// partial adapter, which a full-length 34 bp pattern never matches.
#[test]
fn partial_readthrough_is_detected() {
    let set = AdapterSet::default_for(Mode::Short);
    let stride = set.len();
    let mut seq = vec![b'C'; 50];
    seq.extend_from_slice(b"AGATCGGAAGAGCACACG"); // 18 of 34 adapter bases
    let mut counts = vec![0u64; seq.len() * stride];
    set.scan(&seq, &mut counts, stride);

    assert_eq!(counts[49 * stride], 0, "adapter reported before it starts");
    assert_eq!(counts[50 * stride], 1, "partial adapter not detected");
}

#[test]
fn matches_are_cumulative_to_read_end() {
    let set = AdapterSet::default_for(Mode::Short);
    let stride = set.len();
    let mut seq = vec![b'C'; 20];
    seq.extend_from_slice(b"CTGTCTCTTATA");
    let mut counts = vec![0u64; seq.len() * stride];
    set.scan(&seq, &mut counts, stride);

    let nextera = 3;
    for pos in 0..20 {
        assert_eq!(counts[pos * stride + nextera], 0);
    }
    for pos in 20..seq.len() {
        assert_eq!(counts[pos * stride + nextera], 1);
    }
}

#[test]
fn overlapping_adapters_are_reported_independently() {
    let set = AdapterSet::default_for(Mode::Short);
    let mut seq = Vec::new();
    seq.extend_from_slice(b"GGGGGGGGGGGG");
    seq.extend_from_slice(b"AGATCGGAAGAG");
    let mut hits = vec![false; set.len()];
    set.scan_any(&seq, &mut hits);
    assert!(hits[0], "universal adapter missed");
    assert!(hits[5], "PolyG missed");
}

#[test]
fn long_mode_uses_long_probes_on_both_strands() {
    for (name, probe) in DEFAULT_LONG_ADAPTERS {
        assert!(
            probe.len() >= 20,
            "{name}: long-read probes must be long enough to avoid chance hits"
        );
    }

    let set = AdapterSet::default_for(Mode::Long);
    let rc = revcomp(b"ATCTCTCTCTTTTCCTCCTCCTCC");
    let mut seq = vec![b'A'; 100];
    seq.extend_from_slice(&rc);
    let mut hits = vec![false; set.len()];
    set.scan_any(&seq, &mut hits);
    assert!(hits[3], "reverse-complemented SMRTbell adapter missed");
}

#[test]
fn custom_adapter_lists_load() {
    let dir = std::env::temp_dir().join("kira_qc_adapter_test");
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("adapters.txt");
    std::fs::write(&path, "# custom\nMy Adapter\tACGTACGTACGT\n").unwrap();

    let set = AdapterSet::from_file(&path, false).unwrap();
    assert_eq!(set.len(), 1);
    assert_eq!(set.names(), vec!["My Adapter".to_string()]);

    let mut hits = vec![false; 1];
    set.scan_any(b"TTTTACGTACGTACGTTTTT", &mut hits);
    assert!(hits[0]);

    std::fs::remove_file(&path).ok();
}

#[test]
fn status_thresholds_match_fastqc() {
    assert_eq!(adapter_status(1.0), Status::Pass);
    assert_eq!(adapter_status(7.0), Status::Warn);
    assert_eq!(adapter_status(15.0), Status::Fail);
}
