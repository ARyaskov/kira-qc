use aho_corasick::{AhoCorasick, AhoCorasickBuilder, MatchKind};
use std::sync::OnceLock;

pub const ADAPTERS: [&str; 5] = [
    "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC", // Illumina Universal Adapter
    "TGGAATTCTCGGGTGCCAAGG",              // Illumina Small RNA 3' Adapter
    "GTTCAGAGTTCTACAGTCCGACGATC",         // Illumina Small RNA 5' Adapter
    "CTGTCTCTTATACACATCT",                // Nextera Transposase Sequence
    "CGCCTTGGCCGTACAGCAG",                // SOLiD Small RNA Adapter
];

pub fn adapter_matcher() -> &'static AhoCorasick {
    static AC: OnceLock<AhoCorasick> = OnceLock::new();
    AC.get_or_init(|| {
        AhoCorasickBuilder::new()
            .ascii_case_insensitive(true)
            .match_kind(MatchKind::LeftmostFirst)
            .build(ADAPTERS)
            .expect("adapter automaton")
    })
}

/// FastQC semantics: credit every position from the first match of each adapter to read end.
pub fn scan(seq: &[u8], counts: &mut [[u64; ADAPTERS.len()]]) {
    if seq.is_empty() {
        return;
    }
    let ac = adapter_matcher();
    let mut first_start: [Option<usize>; ADAPTERS.len()] = [None; ADAPTERS.len()];
    let mut remaining = ADAPTERS.len();
    for mat in ac.find_iter(seq) {
        let idx = mat.pattern().as_usize();
        if first_start[idx].is_none() {
            first_start[idx] = Some(mat.start());
            remaining -= 1;
            if remaining == 0 {
                break;
            }
        }
    }
    let cap = counts.len().min(seq.len());
    for (idx, start) in first_start.iter().enumerate() {
        let Some(start) = *start else { continue };
        if start >= cap {
            continue;
        }
        for c in counts[start..cap].iter_mut() {
            c[idx] += 1;
        }
    }
}

pub fn scan_any(seq: &[u8], hits: &mut [bool; ADAPTERS.len()]) {
    if seq.is_empty() {
        return;
    }
    let ac = adapter_matcher();
    for mat in ac.find_iter(seq) {
        hits[mat.pattern().as_usize()] = true;
    }
}

#[derive(Clone, Debug)]
pub struct AdapterRow {
    pub position: usize,
    pub values: [f64; ADAPTERS.len()],
}
