use rustc_hash::FxHashMap;
use std::cmp::Reverse;
use std::collections::BinaryHeap;

const OVERREP_K: usize = 200_000;
const MAX_SEQ_LEN: usize = 150;
const HEAP_COMPACT_THRESHOLD: usize = 4 * OVERREP_K;

#[derive(Clone, Debug)]
pub struct OverrepRow {
    pub sequence: String,
    pub count: u64,
    pub percent: f64,
    pub source: &'static str,
}

#[derive(Clone, Debug)]
pub struct Entry {
    pub key: u64,
    pub count: u64,
    pub seq: Vec<u8>,
}

#[derive(Clone, Debug)]
pub struct SpaceSavingSeq {
    map: FxHashMap<u64, u32>,
    entries: Vec<Entry>,
    heap: BinaryHeap<(Reverse<u64>, u64, u32)>,
}

impl Default for SpaceSavingSeq {
    fn default() -> Self {
        Self::new()
    }
}

impl SpaceSavingSeq {
    pub fn new() -> Self {
        Self {
            map: FxHashMap::with_capacity_and_hasher(OVERREP_K, Default::default()),
            entries: Vec::with_capacity(OVERREP_K),
            heap: BinaryHeap::with_capacity(OVERREP_K),
        }
    }

    pub fn add(&mut self, key: u64, seq: &[u8], weight: u64) {
        if let Some(&idx) = self.map.get(&key) {
            let e = &mut self.entries[idx as usize];
            e.count += weight;
            self.heap.push((Reverse(e.count), e.key, idx));
            self.maybe_compact_heap();
            return;
        }

        if self.entries.len() < OVERREP_K {
            let idx = self.entries.len() as u32;
            self.entries.push(Entry {
                key,
                count: weight,
                seq: trim_seq(seq),
            });
            self.map.insert(key, idx);
            self.heap.push((Reverse(weight), key, idx));
            return;
        }

        let (min_idx, min_count) = self.min_entry();
        let removed = self.entries[min_idx as usize].key;
        self.map.remove(&removed);
        let new_count = min_count + weight;
        self.entries[min_idx as usize] = Entry {
            key,
            count: new_count,
            seq: trim_seq(seq),
        };
        self.map.insert(key, min_idx);
        self.heap.push((Reverse(new_count), key, min_idx));
        self.maybe_compact_heap();
    }

    pub fn merge(&mut self, other: &SpaceSavingSeq) {
        let mut idx: Vec<usize> = (0..other.entries.len()).collect();
        idx.sort_unstable_by_key(|&i| other.entries[i].key);
        for i in idx {
            let e = &other.entries[i];
            self.add(e.key, &e.seq, e.count);
        }
    }

    pub fn entries(&self) -> &[Entry] {
        &self.entries
    }

    fn min_entry(&mut self) -> (u32, u64) {
        loop {
            if let Some((Reverse(count), key, idx)) = self.heap.pop() {
                let e = &self.entries[idx as usize];
                if e.key == key && e.count == count {
                    return (idx, count);
                }
            } else {
                return (0, self.entries[0].count);
            }
        }
    }

    fn maybe_compact_heap(&mut self) {
        if self.heap.len() < HEAP_COMPACT_THRESHOLD {
            return;
        }
        self.heap.clear();
        self.heap.reserve(self.entries.len());
        for (i, e) in self.entries.iter().enumerate() {
            self.heap.push((Reverse(e.count), e.key, i as u32));
        }
    }
}

fn trim_seq(seq: &[u8]) -> Vec<u8> {
    if seq.len() <= MAX_SEQ_LEN {
        return seq.to_vec();
    }
    seq[..MAX_SEQ_LEN].to_vec()
}

#[inline]
pub fn hash_seq(seq: &[u8]) -> u64 {
    let mut buf: [u8; 256] = [0; 256];
    if seq.len() <= buf.len() {
        for (i, &b) in seq.iter().enumerate() {
            buf[i] = b & 0xDF;
        }
        return xxhash_rust::xxh3::xxh3_64(&buf[..seq.len()]);
    }
    let folded: Vec<u8> = seq.iter().map(|&b| b & 0xDF).collect();
    xxhash_rust::xxh3::xxh3_64(&folded)
}

pub fn classify_source(seq: &[u8]) -> &'static str {
    if is_poly(seq, b'A') {
        return "Poly-A";
    }
    if is_poly(seq, b'T') {
        return "Poly-T";
    }
    if contains_adapter(seq) {
        return "Adapter";
    }
    "No Hit"
}

fn is_poly(seq: &[u8], base: u8) -> bool {
    if seq.len() < 20 {
        return false;
    }
    seq.iter().all(|&b| (b & 0xDF) == base)
}

fn contains_adapter(seq: &[u8]) -> bool {
    use aho_corasick::{AhoCorasick, AhoCorasickBuilder};
    use std::sync::OnceLock;

    static AC: OnceLock<AhoCorasick> = OnceLock::new();
    let ac = AC.get_or_init(|| {
        AhoCorasickBuilder::new()
            .ascii_case_insensitive(true)
            .build([
                "AGATCGGAAGAG", // Illumina Universal Adapter prefix
                "TGGAATTCTCGG", // Small RNA 3'
                "ATCTCGTATGCC", // Small RNA 5'
                "CTGTCTCTTATA", // Nextera
                "CGCCTTGGCCGT", // SOLiD
            ])
            .expect("classify adapter automaton")
    });
    ac.is_match(seq)
}
