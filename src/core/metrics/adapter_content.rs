use crate::simd;
use aho_corasick::{AhoCorasick, AhoCorasickBuilder};
use std::sync::OnceLock;

pub const ADAPTERS: [&str; 5] = [
    "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC", // Illumina Universal Adapter
    "TGGAATTCTCGGGTGCCAAGG",              // Illumina Small RNA 3' Adapter
    "GTTCAGAGTTCTACAGTCCGACGATC",         // Illumina Small RNA 5' Adapter
    "CTGTCTCTTATACACATCT",                // Nextera Transposase Sequence
    "CGCCTTGGCCGTACAGCAG",                // SOLiD Small RNA Adapter
];

const PREFIXES: [&[u8]; 5] = [
    b"AGATCGGA",
    b"TGGAATTC",
    b"GTTCAGAG",
    b"CTGTCTCT",
    b"CGCCTTGG",
];

pub fn adapter_matcher() -> &'static AhoCorasick {
    static AC: OnceLock<AhoCorasick> = OnceLock::new();
    AC.get_or_init(|| {
        AhoCorasickBuilder::new()
            .ascii_case_insensitive(true)
            .build(ADAPTERS)
            .expect("adapter automaton")
    })
}

pub fn scan(seq: &[u8], counts: &mut [[u64; ADAPTERS.len()]]) {
    if seq.is_empty() {
        return;
    }
    if !prefilter(seq) {
        return;
    }
    let ac = adapter_matcher();
    for mat in ac.find_iter(seq) {
        let pos = mat.start();
        if pos < counts.len() {
            let idx = mat.pattern().as_usize();
            counts[pos][idx] += 1;
        }
    }
}

pub fn scan_any(seq: &[u8], hits: &mut [bool; ADAPTERS.len()]) {
    if seq.is_empty() {
        return;
    }
    if !prefilter(seq) {
        return;
    }
    let ac = adapter_matcher();
    for mat in ac.find_iter(seq) {
        let idx = mat.pattern().as_usize();
        hits[idx] = true;
    }
}

fn prefilter(seq: &[u8]) -> bool {
    for p in PREFIXES {
        if simd::prefix_scan(seq, p) {
            return true;
        }
    }
    false
}

#[derive(Clone, Debug)]
pub struct AdapterRow {
    pub position: usize,
    pub values: [f64; ADAPTERS.len()],
}
