//! Source attribution for overrepresented sequences.
//!
//! The counting itself lives in [`crate::core::metrics::duplication`] - FastQC
//! derives both modules from one exact sequence table, and so do we. This file
//! only answers "what is this sequence?".

use crate::core::metrics::adapter_content::parse_named_sequences;
pub use crate::core::metrics::duplication::OverrepRow;
use anyhow::{Context, Result};
use std::path::Path;
use std::sync::{Arc, OnceLock};

/// Subset of FastQC's `contaminant_list.txt` covering the sources that actually
/// show up in Illumina, ONT and PacBio libraries.
pub const DEFAULT_CONTAMINANTS: [(&str, &str); 30] = [
    (
        "Illumina Single End Adapter 1",
        "GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "Illumina Single End PCR Primer 1",
        "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    ),
    (
        "Illumina Paired End Adapter 1",
        "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    ),
    (
        "Illumina Paired End Adapter 2",
        "GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG",
    ),
    (
        "Illumina Paired End PCR Primer 2",
        "CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT",
    ),
    (
        "Illumina TruSeq Adapter Read 1",
        "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
    ),
    (
        "Illumina TruSeq Adapter Read 2",
        "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    ),
    (
        "Illumina Universal Adapter",
        "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
    ),
    (
        "TruSeq Universal Adapter",
        "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    ),
    (
        "Nextera Transposase Sequence Read 1",
        "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
    ),
    (
        "Nextera Transposase Sequence Read 2",
        "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
    ),
    ("Nextera Mate Pair Adapter", "CTGTCTCTTATACACATCT"),
    ("Illumina Small RNA 3' Adapter", "TGGAATTCTCGGGTGCCAAGG"),
    (
        "Illumina Small RNA 5' Adapter",
        "GTTCAGAGTTCTACAGTCCGACGATC",
    ),
    (
        "Illumina Multiplexing Adapter 1",
        "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
    ),
    (
        "Illumina Multiplexing Adapter 2",
        "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    ),
    (
        "Illumina Multiplexing PCR Primer 2.01",
        "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
    ),
    (
        "Illumina DpnII expression Adapter 1",
        "ACAGGTTCAGAGTTCTACAGTCCGACGATC",
    ),
    (
        "Illumina NlaIII expression Adapter 1",
        "TCGGACTGTAGAACTCTGAAC",
    ),
    (
        "PhiX174 control spike-in",
        "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTTGATAAAGCAGGAATTACTACTGCTTGTTTACGAATT",
    ),
    (
        "ONT Ligation Adapter (LSK top)",
        "AATGTACTTCGTTCAGTTACGTATTGCT",
    ),
    (
        "ONT Ligation Adapter (LSK bottom)",
        "GCAATACGTAACTGAACGAAGT",
    ),
    (
        "ONT Rapid Adapter (RAD)",
        "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
    ),
    (
        "ONT Native Barcode Flank",
        "GGTGCTGAAGAAAGTTGTCGGTGTCTTTGTG",
    ),
    (
        "PacBio SMRTbell Adapter",
        "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT",
    ),
    (
        "PacBio SMRTbell Adapter (v3)",
        "ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT",
    ),
    ("Oligo dT primer", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"),
    ("PolyA", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
    ("PolyT", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"),
    ("PolyG", "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"),
];

/// Shortest exact overlap that counts as a hit.
const MIN_OVERLAP: usize = 20;

#[derive(Clone, Debug)]
pub struct ContaminantList {
    entries: Vec<(String, Vec<u8>)>,
}

impl ContaminantList {
    pub fn builtin() -> Arc<ContaminantList> {
        static C: OnceLock<Arc<ContaminantList>> = OnceLock::new();
        C.get_or_init(|| {
            Arc::new(ContaminantList {
                entries: DEFAULT_CONTAMINANTS
                    .iter()
                    .map(|(name, seq)| ((*name).to_string(), seq.as_bytes().to_vec()))
                    .collect(),
            })
        })
        .clone()
    }

    /// Load a FastQC-format contaminant list: `name<TAB>sequence`, `#` comments.
    pub fn from_file(path: &Path) -> Result<Arc<ContaminantList>> {
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("failed to read contaminant list {}", path.display()))?;
        let entries = parse_named_sequences(&text)
            .with_context(|| format!("failed to parse contaminant list {}", path.display()))?
            .into_iter()
            .map(|(name, seq)| (name, seq.into_bytes()))
            .collect();
        Ok(Arc::new(ContaminantList { entries }))
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// FastQC-style attribution: report the contaminant with the longest exact
    /// overlap, along with how much of it matched.
    pub fn classify(&self, seq: &[u8]) -> String {
        if seq.is_empty() {
            return "No Hit".to_string();
        }
        if let Some(base) = homopolymer_base(seq) {
            return format!("Poly-{}", base as char);
        }

        let mut best: Option<(&str, usize, usize)> = None;
        for (name, contaminant) in &self.entries {
            let overlap = longest_common_run(seq, contaminant);
            let required = MIN_OVERLAP.min(contaminant.len());
            if overlap < required {
                continue;
            }
            let better = match best {
                None => true,
                Some((_, best_overlap, _)) => overlap > best_overlap,
            };
            if better {
                best = Some((name, overlap, contaminant.len()));
            }
        }

        match best {
            Some((name, overlap, c_len)) => {
                let pct = overlap as f64 * 100.0 / c_len as f64;
                format!("{} ({:.0}% over {}bp)", name, pct, overlap)
            }
            None => "No Hit".to_string(),
        }
    }
}

/// Longest exact common substring between two sequences, case-folded.
fn longest_common_run(a: &[u8], b: &[u8]) -> usize {
    if a.is_empty() || b.is_empty() {
        return 0;
    }
    let mut best = 0usize;
    // Slide `b` across `a` and measure the longest run of equal bases at each
    // offset. Both inputs are at most a few dozen bases and this only runs for
    // the handful of sequences that get reported.
    let a_len = a.len() as isize;
    let b_len = b.len() as isize;
    for offset in -(b_len - 1)..a_len {
        let mut run = 0usize;
        let start = offset.max(0);
        let end = (offset + b_len).min(a_len);
        for i in start..end {
            let bi = i - offset;
            if (a[i as usize] & 0xDF) == (b[bi as usize] & 0xDF) {
                run += 1;
                if run > best {
                    best = run;
                }
            } else {
                run = 0;
            }
        }
    }
    best
}

/// Reads made almost entirely of one base, allowing a little noise.
fn homopolymer_base(seq: &[u8]) -> Option<u8> {
    if seq.len() < MIN_OVERLAP {
        return None;
    }
    for base in *b"ACGT" {
        let matches = seq.iter().filter(|&&b| (b & 0xDF) == base).count();
        if matches * 100 >= seq.len() * 90 {
            return Some(base);
        }
    }
    None
}

/// FastQC: warn above 0.1 % of reads, fail above 1 %.
pub fn overrep_status(max_percent: f64) -> crate::core::model::Status {
    use crate::core::model::Status;
    if max_percent > 1.0 {
        Status::Fail
    } else if max_percent > 0.1 {
        Status::Warn
    } else {
        Status::Pass
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn polya_is_recognised() {
        let seq = vec![b'A'; 40];
        assert_eq!(ContaminantList::builtin().classify(&seq), "Poly-A");
    }

    #[test]
    fn adapter_is_attributed_with_overlap() {
        let source = ContaminantList::builtin().classify(b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC");
        assert!(
            source.contains("Adapter") || source.contains("Universal"),
            "unexpected source: {source}"
        );
        assert!(source.contains("bp"), "overlap not reported: {source}");
    }

    #[test]
    fn random_sequence_has_no_hit() {
        let seq = b"CGTAACGTGCATTAGCCATGACTGCAAGTCCGATTACGGA";
        assert_eq!(ContaminantList::builtin().classify(seq), "No Hit");
    }

    #[test]
    fn longest_common_run_finds_the_overlap() {
        assert_eq!(longest_common_run(b"AAACGTACGTAAA", b"CGTACGT"), 7);
        assert_eq!(longest_common_run(b"AAAA", b"CCCC"), 0);
        assert_eq!(longest_common_run(b"acgtacgt", b"ACGTACGT"), 8);
    }

    #[test]
    fn ont_adapter_is_attributed() {
        let source = ContaminantList::builtin().classify(b"AATGTACTTCGTTCAGTTACGTATTGCTAAAA");
        assert!(source.contains("ONT"), "unexpected source: {source}");
    }

    #[test]
    fn status_thresholds_match_fastqc() {
        use crate::core::model::Status;
        assert_eq!(overrep_status(0.05), Status::Pass);
        assert_eq!(overrep_status(0.5), Status::Warn);
        assert_eq!(overrep_status(2.0), Status::Fail);
    }
}
