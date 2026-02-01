use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};

// Chosen to capture frequent contaminants without unbounded memory.
const OVERREP_K: usize = 200_000;
const MAX_SEQ_LEN: usize = 150;

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
    pub error: u64,
    pub seq: Vec<u8>,
}

#[derive(Clone, Debug)]
pub struct SpaceSavingSeq {
    map: HashMap<u64, usize>,
    entries: Vec<Entry>,
    heap: BinaryHeap<(Reverse<u64>, u64, usize)>,
}

impl SpaceSavingSeq {
    pub fn new() -> Self {
        Self {
            map: HashMap::with_capacity(OVERREP_K),
            entries: Vec::with_capacity(OVERREP_K),
            heap: BinaryHeap::with_capacity(OVERREP_K),
        }
    }

    pub fn add(&mut self, key: u64, seq: &[u8], weight: u64) {
        if let Some(&idx) = self.map.get(&key) {
            let e = &mut self.entries[idx];
            e.count += weight;
            self.heap.push((Reverse(e.count), e.key, idx));
            return;
        }

        if self.entries.len() < OVERREP_K {
            let idx = self.entries.len();
            self.entries.push(Entry {
                key,
                count: weight,
                error: 0,
                seq: trim_seq(seq),
            });
            self.map.insert(key, idx);
            self.heap.push((Reverse(weight), key, idx));
            return;
        }

        let (min_idx, min_count) = self.min_entry();
        let removed = self.entries[min_idx].key;
        self.map.remove(&removed);
        self.entries[min_idx] = Entry {
            key,
            count: min_count + weight,
            error: min_count,
            seq: trim_seq(seq),
        };
        self.map.insert(key, min_idx);
        self.heap.push((Reverse(min_count + weight), key, min_idx));
    }

    pub fn merge(&mut self, other: &SpaceSavingSeq) {
        let mut items = other.entries.clone();
        items.sort_by_key(|e| e.key);
        for e in items {
            self.add(e.key, &e.seq, e.count);
        }
    }

    pub fn entries(&self) -> &[Entry] {
        &self.entries
    }

    fn min_entry(&mut self) -> (usize, u64) {
        loop {
            if let Some((Reverse(count), key, idx)) = self.heap.pop() {
                let e = &self.entries[idx];
                if e.key == key && e.count == count {
                    return (idx, count);
                }
            } else {
                return (0, self.entries[0].count);
            }
        }
    }
}

fn trim_seq(seq: &[u8]) -> Vec<u8> {
    if seq.len() <= MAX_SEQ_LEN {
        return seq.to_vec();
    }
    seq[..MAX_SEQ_LEN].to_vec()
}

pub fn hash_seq(seq: &[u8]) -> u64 {
    const FNV_OFFSET: u64 = 0xcbf29ce484222325;
    const FNV_PRIME: u64 = 0x100000001b3;
    let mut h = FNV_OFFSET;
    for &b in seq {
        h ^= (b & 0xDF) as u64;
        h = h.wrapping_mul(FNV_PRIME);
    }
    h
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
    for &b in seq {
        let u = b & 0xDF;
        if u != base {
            return false;
        }
    }
    true
}

fn contains_adapter(seq: &[u8]) -> bool {
    const ADAPTERS: [&[u8]; 5] = [
        b"AGATCGGAAGAG", // Illumina Universal Adapter prefix
        b"TGGAATTCTCGG", // Small RNA 3'
        b"ATCTCGTATGCC", // Small RNA 5'
        b"CTGTCTCTTATA", // Nextera
        b"CGCCTTGGCCGT", // SOLiD
    ];
    for a in ADAPTERS {
        if find_subseq(seq, a) {
            return true;
        }
    }
    false
}

fn find_subseq(hay: &[u8], needle: &[u8]) -> bool {
    if needle.is_empty() || hay.len() < needle.len() {
        return false;
    }
    for i in 0..=hay.len() - needle.len() {
        let mut ok = true;
        for j in 0..needle.len() {
            if (hay[i + j] & 0xDF) != needle[j] {
                ok = false;
                break;
            }
        }
        if ok {
            return true;
        }
    }
    false
}
