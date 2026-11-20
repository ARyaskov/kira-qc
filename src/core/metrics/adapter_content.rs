//! Adapter content.
//!
//! FastQC matches short adapter *prefixes* (12 bp), not whole adapter
//! sequences. Read-through typically leaves only the first 10-25 bases of the
//! adapter at the 3' end of a read, so matching a full 34 bp adapter - as this
//! module used to - never fires on real data and silently under-reports
//! contamination.

use crate::core::model::{BaseGroup, Mode, Status};
use aho_corasick::{AhoCorasick, AhoCorasickBuilder, MatchKind};
use anyhow::{Context, Result, bail};
use std::collections::HashSet;
use std::path::Path;
use std::sync::{Arc, OnceLock};

/// FastQC 0.12.x `adapter_list.txt`. SOLiD was dropped upstream; PolyA and
/// PolyG were added.
pub const DEFAULT_SHORT_ADAPTERS: [(&str, &str); 6] = [
    ("Illumina Universal Adapter", "AGATCGGAAGAG"),
    ("Illumina Small RNA 3' Adapter", "TGGAATTCTCGG"),
    ("Illumina Small RNA 5' Adapter", "GATCGTCGGACT"),
    ("Nextera Transposase Sequence", "CTGTCTCTTATA"),
    ("PolyA", "AAAAAAAAAAAA"),
    ("PolyG", "GGGGGGGGGGGG"),
];

/// Long reads need longer probes: a 12-mer occurs by chance roughly once every
/// 16 Mbp, which on 100 kb reads would produce a steady stream of false
/// positives. These probes are 20-28 bp and are matched on both strands.
pub const DEFAULT_LONG_ADAPTERS: [(&str, &str); 6] = [
    ("ONT Ligation Adapter (LSK)", "AATGTACTTCGTTCAGTTACGTATTGCT"),
    ("ONT Ligation Adapter (bottom)", "GCAATACGTAACTGAACGAAGT"),
    ("ONT Rapid Adapter (RAD)", "GTTTTCGCATTTATCGTGAAACGC"),
    ("PacBio SMRTbell Adapter", "ATCTCTCTCTTTTCCTCCTCCTCC"),
    ("PolyA", "AAAAAAAAAAAAAAAAAAAA"),
    ("PolyG", "GGGGGGGGGGGGGGGGGGGG"),
];

/// Upper bound on adapters, so per-read scratch can stay on the stack.
pub const MAX_ADAPTERS: usize = 32;

#[derive(Clone, Debug)]
pub struct Adapter {
    pub name: String,
    pub probe: String,
}

pub struct AdapterSet {
    adapters: Vec<Adapter>,
    ac: AhoCorasick,
    /// Maps automaton pattern index to adapter index; a reverse complement maps
    /// back onto the adapter it came from.
    pattern_to_adapter: Vec<usize>,
}

impl std::fmt::Debug for AdapterSet {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("AdapterSet")
            .field("adapters", &self.adapters)
            .finish()
    }
}

impl AdapterSet {
    pub fn new(adapters: Vec<Adapter>, with_revcomp: bool) -> Result<Self> {
        if adapters.is_empty() {
            bail!("adapter list is empty");
        }
        if adapters.len() > MAX_ADAPTERS {
            bail!(
                "too many adapters: {} (limit {})",
                adapters.len(),
                MAX_ADAPTERS
            );
        }
        let mut names = HashSet::with_capacity(adapters.len());
        for a in &adapters {
            if !names.insert(a.name.as_str()) {
                bail!("duplicate adapter name '{}'", a.name);
            }
            if a.probe.is_empty() {
                bail!("adapter '{}' has an empty sequence", a.name);
            }
            if !a
                .probe
                .bytes()
                .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N'))
            {
                bail!(
                    "adapter '{}' contains unsupported symbols (allowed: A, C, G, T, literal N): {}",
                    a.name,
                    a.probe
                );
            }
        }

        let mut patterns: Vec<String> = Vec::new();
        let mut pattern_to_adapter: Vec<usize> = Vec::new();
        for (i, a) in adapters.iter().enumerate() {
            patterns.push(a.probe.clone());
            pattern_to_adapter.push(i);
            if with_revcomp {
                let rc = String::from_utf8(revcomp(a.probe.as_bytes()))
                    .expect("revcomp of ASCII stays ASCII");
                if rc != a.probe {
                    patterns.push(rc);
                    pattern_to_adapter.push(i);
                }
            }
        }

        let ac = AhoCorasickBuilder::new()
            .ascii_case_insensitive(true)
            // `Standard` is required for overlapping search: with `LeftmostFirst`
            // an adapter starting inside another adapter's match is never reported.
            .match_kind(MatchKind::Standard)
            .build(&patterns)
            .context("failed to build adapter automaton")?;

        Ok(Self {
            adapters,
            ac,
            pattern_to_adapter,
        })
    }

    pub fn default_for(mode: Mode) -> Arc<AdapterSet> {
        match mode {
            Mode::Short => {
                static S: OnceLock<Arc<AdapterSet>> = OnceLock::new();
                S.get_or_init(|| {
                    Arc::new(
                        AdapterSet::new(to_adapters(&DEFAULT_SHORT_ADAPTERS), false)
                            .expect("built-in short adapter list"),
                    )
                })
                .clone()
            }
            Mode::Long => {
                static L: OnceLock<Arc<AdapterSet>> = OnceLock::new();
                L.get_or_init(|| {
                    Arc::new(
                        AdapterSet::new(to_adapters(&DEFAULT_LONG_ADAPTERS), true)
                            .expect("built-in long adapter list"),
                    )
                })
                .clone()
            }
        }
    }

    /// Load a FastQC-format adapter list: `name<TAB>sequence`, `#` comments.
    pub fn from_file(path: &Path, with_revcomp: bool) -> Result<Arc<AdapterSet>> {
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("failed to read adapter list {}", path.display()))?;
        let adapters = parse_named_sequences(&text)
            .with_context(|| format!("failed to parse adapter list {}", path.display()))?
            .into_iter()
            .map(|(name, probe)| Adapter { name, probe })
            .collect();
        Ok(Arc::new(AdapterSet::new(adapters, with_revcomp)?))
    }

    pub fn len(&self) -> usize {
        self.adapters.len()
    }

    pub fn is_empty(&self) -> bool {
        self.adapters.is_empty()
    }

    pub fn adapters(&self) -> &[Adapter] {
        &self.adapters
    }

    pub fn names(&self) -> Vec<String> {
        self.adapters.iter().map(|a| a.name.clone()).collect()
    }

    /// First match position of each adapter, or `None`.
    fn first_matches(&self, seq: &[u8], out: &mut [Option<usize>]) {
        let mut remaining = self.adapters.len();
        for mat in self.ac.find_overlapping_iter(seq) {
            let idx = self.pattern_to_adapter[mat.pattern().as_usize()];
            if out[idx].is_none() {
                out[idx] = Some(mat.start());
                remaining -= 1;
                if remaining == 0 {
                    break;
                }
            }
        }
    }

    /// FastQC semantics: once an adapter is seen at position `p`, every position
    /// from `p` to the end of the read counts as adapter-contaminated.
    ///
    /// `counts` is a flat `positions * stride` matrix.
    pub fn scan(&self, seq: &[u8], counts: &mut [u64], stride: usize) {
        if seq.is_empty() || stride == 0 {
            return;
        }
        let mut first_start = [None; MAX_ADAPTERS];
        let n = self.adapters.len();
        self.first_matches(seq, &mut first_start[..n]);

        let positions = counts.len() / stride;
        let cap = positions.min(seq.len());
        for (idx, start) in first_start[..n].iter().enumerate() {
            let Some(start) = *start else { continue };
            if start >= cap {
                continue;
            }
            for pos in start..cap {
                counts[pos * stride + idx] += 1;
            }
        }
    }

    /// Record only the boundaries of each contaminated interval.
    pub(super) fn scan_ranges(&self, seq: &[u8], deltas: &mut [i64], stride: usize) {
        let mut first_start = [None; MAX_ADAPTERS];
        let starts = &mut first_start[..self.adapters.len()];
        self.first_matches(seq, starts);
        for (adapter, start) in starts.iter().enumerate() {
            if let Some(start) = *start {
                deltas[start * stride + adapter] += 1;
                deltas[seq.len() * stride + adapter] -= 1;
            }
        }
    }

    /// Presence-only scan for long reads, where per-position curves are meaningless.
    pub fn scan_any(&self, seq: &[u8], hits: &mut [bool]) {
        if seq.is_empty() {
            return;
        }
        let mut first_start = [None; MAX_ADAPTERS];
        let n = self.adapters.len();
        self.first_matches(seq, &mut first_start[..n]);
        for (idx, start) in first_start[..n].iter().enumerate() {
            if start.is_some() {
                hits[idx] = true;
            }
        }
    }
}

fn to_adapters(list: &[(&str, &str)]) -> Vec<Adapter> {
    list.iter()
        .map(|(name, probe)| Adapter {
            name: (*name).to_string(),
            probe: (*probe).to_string(),
        })
        .collect()
}

/// Shared parser for FastQC's `name<whitespace>sequence` list files.
pub fn parse_named_sequences(text: &str) -> Result<Vec<(String, String)>> {
    let mut out = Vec::new();
    for (lineno, line) in text.lines().enumerate() {
        let line = line.trim_end_matches(['\r', '\n']);
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        // FastQC separates the name from the sequence with a tab (often several).
        let (name, seq) = match line.rsplit_once('\t') {
            Some((n, s)) => (n, s),
            None => match trimmed.rsplit_once(char::is_whitespace) {
                Some((n, s)) => (n, s),
                None => bail!("line {}: expected 'name<TAB>sequence'", lineno + 1),
            },
        };
        let name = name.trim();
        let seq = seq.trim().to_ascii_uppercase();
        if name.is_empty() || seq.is_empty() {
            bail!("line {}: expected 'name<TAB>sequence'", lineno + 1);
        }
        out.push((name.to_string(), seq));
    }
    if out.is_empty() {
        bail!("no entries found");
    }
    Ok(out)
}

pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b & 0xDF {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            other => other,
        })
        .collect()
}

#[derive(Clone, Debug)]
pub struct AdapterRow {
    pub group: BaseGroup,
    pub values: Vec<f64>,
}

/// FastQC: warn above 5 % at any position, fail above 10 %.
pub fn adapter_status(max_percent: f64) -> Status {
    if max_percent > 10.0 {
        Status::Fail
    } else if max_percent > 5.0 {
        Status::Warn
    } else {
        Status::Pass
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn short() -> Arc<AdapterSet> {
        AdapterSet::default_for(Mode::Short)
    }

    #[test]
    fn short_probes_match_partial_readthrough() {
        // 18 bases of the universal adapter: the whole 34 bp adapter is absent,
        // which is exactly the case the old full-length matcher missed.
        let set = short();
        let mut seq = vec![b'C'; 40];
        seq.extend_from_slice(b"AGATCGGAAGAGCACACG");
        let stride = set.len();
        let mut counts = vec![0u64; seq.len() * stride];
        set.scan(&seq, &mut counts, stride);
        assert_eq!(counts[39 * stride], 0);
        assert_eq!(counts[40 * stride], 1);
        assert_eq!(counts[(seq.len() - 1) * stride], 1);
    }

    #[test]
    fn overlapping_adapters_are_all_found() {
        let set = short();
        let mut seq = Vec::new();
        seq.extend_from_slice(b"AAAAAAAAAAAA");
        seq.extend_from_slice(b"CTGTCTCTTATA");
        let mut hits = vec![false; set.len()];
        set.scan_any(&seq, &mut hits);
        assert!(hits[4], "PolyA not detected");
        assert!(hits[3], "Nextera not detected");
    }

    #[test]
    fn cumulative_from_first_match_to_read_end() {
        let set = short();
        let stride = set.len();
        let mut seq = vec![b'C'; 30];
        seq.extend_from_slice(b"AGATCGGAAGAG");
        let mut counts = vec![0u64; seq.len() * stride];
        set.scan(&seq, &mut counts, stride);
        for pos in 0..30 {
            assert_eq!(counts[pos * stride], 0);
        }
        for pos in 30..seq.len() {
            assert_eq!(counts[pos * stride], 1);
        }
    }

    #[test]
    fn only_first_occurrence_counts() {
        let set = short();
        let stride = set.len();
        let mut seq = Vec::new();
        seq.extend_from_slice(b"AGATCGGAAGAG");
        seq.extend_from_slice(&[b'N'; 5]);
        seq.extend_from_slice(b"AGATCGGAAGAG");
        let mut counts = vec![0u64; seq.len() * stride];
        set.scan(&seq, &mut counts, stride);
        for pos in 0..seq.len() {
            assert_eq!(counts[pos * stride], 1);
        }
    }

    #[test]
    fn interval_boundaries_match_dense_counts_for_mixed_lengths() {
        let set = short();
        let stride = set.len();
        let sequences: [&[u8]; 4] = [
            b"AGATCGGAAGAGCTGTCTCTTATA",
            b"CCAGATCGGAAGAG",
            b"GGGGGGGGGGGGAC",
            b"AAAAAAAAAAAA",
        ];
        let length = sequences.iter().map(|seq| seq.len()).max().unwrap();
        let mut dense = vec![0; length * stride];
        let mut deltas = vec![0; (length + 1) * stride];
        for seq in sequences {
            set.scan(seq, &mut dense, stride);
            set.scan_ranges(seq, &mut deltas, stride);
        }
        let mut counts = vec![0; stride];
        for (position, row) in deltas.chunks_exact(stride).enumerate() {
            for (adapter, &delta) in row.iter().enumerate() {
                counts[adapter] += delta;
                if position < length {
                    assert_eq!(counts[adapter] as u64, dense[position * stride + adapter]);
                } else {
                    assert_eq!(counts[adapter], 0);
                }
            }
        }
    }

    #[test]
    fn long_adapters_match_both_strands() {
        let set = AdapterSet::default_for(Mode::Long);
        let fwd = b"TTTTTTAATGTACTTCGTTCAGTTACGTATTGCTTTTTTT";
        let mut hits = vec![false; set.len()];
        set.scan_any(fwd, &mut hits);
        assert!(hits[0]);

        let rc = revcomp(b"AATGTACTTCGTTCAGTTACGTATTGCT");
        let mut seq = vec![b'T'; 10];
        seq.extend_from_slice(&rc);
        let mut hits2 = vec![false; set.len()];
        set.scan_any(&seq, &mut hits2);
        assert!(hits2[0], "reverse-complemented ONT adapter missed");
    }

    #[test]
    fn revcomp_is_correct() {
        assert_eq!(revcomp(b"ACGT"), b"ACGT".to_vec());
        assert_eq!(revcomp(b"AAGG"), b"CCTT".to_vec());
    }

    #[test]
    fn case_insensitive() {
        let set = short();
        let mut hits = vec![false; set.len()];
        set.scan_any(b"cccagatcggaagagccc", &mut hits);
        assert!(hits[0]);
    }

    #[test]
    fn parses_fastqc_adapter_files() {
        let text = "# comment\nIllumina Universal Adapter\tAGATCGGAAGAG\nNextera\tCTGTCTCTTATA\n";
        let parsed = parse_named_sequences(text).unwrap();
        assert_eq!(parsed.len(), 2);
        assert_eq!(parsed[0].0, "Illumina Universal Adapter");
        assert_eq!(parsed[0].1, "AGATCGGAAGAG");
    }

    #[test]
    fn rejects_duplicate_names_in_adapter_files() {
        let file = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(file.path(), "shared\tAAAA\nshared\tCCCC\n").unwrap();
        let error = AdapterSet::from_file(file.path(), false).unwrap_err();
        assert!(
            error
                .to_string()
                .contains("duplicate adapter name 'shared'")
        );
    }

    #[test]
    fn rejects_unsupported_probe_symbols_in_both_modes() {
        for reverse_complement in [false, true] {
            for probe in ["ACGZ", "ACGR", "ACGU"] {
                let adapter = Adapter {
                    name: "probe".into(),
                    probe: probe.into(),
                };
                let error = AdapterSet::new(vec![adapter], reverse_complement).unwrap_err();
                assert!(error.to_string().contains("unsupported symbols"));
            }
            let adapter = Adapter {
                name: "probe".into(),
                probe: "acgn".into(),
            };
            let set = AdapterSet::new(vec![adapter], reverse_complement).unwrap();
            let mut hits = [false];
            set.scan_any(b"ACGA", &mut hits);
            assert!(!hits[0]);
            set.scan_any(b"ACGN", &mut hits);
            assert!(hits[0]);
        }
    }

    #[test]
    fn rejects_invalid_adapter_files() {
        assert!(parse_named_sequences("# only comments\n").is_err());
        let bad = vec![Adapter {
            name: "x".into(),
            probe: "ACGT123".into(),
        }];
        assert!(AdapterSet::new(bad, false).is_err());
    }
}
