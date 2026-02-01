use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};

// Chosen to balance accuracy and memory; ~200k entries keeps top duplicates without unbounded growth.
const DUP_K: usize = 200_000;

#[derive(Clone, Debug)]
pub struct DuplicationRow {
    pub level: DupLevel,
    pub relative: f64,
}

#[derive(Clone, Copy, Debug)]
pub enum DupLevel {
    One,
    Two,
    Three,
    Four,
    Five,
    Six,
    SevenPlus,
}

impl DupLevel {
    pub fn as_str(self) -> &'static str {
        match self {
            DupLevel::One => "1",
            DupLevel::Two => "2",
            DupLevel::Three => "3",
            DupLevel::Four => "4",
            DupLevel::Five => "5",
            DupLevel::Six => "6",
            DupLevel::SevenPlus => "7+",
        }
    }
}

#[derive(Clone, Debug)]
pub struct Entry {
    pub key: u64,
    pub count: u64,
    pub error: u64,
}

#[derive(Clone, Debug)]
pub struct SpaceSaving {
    map: HashMap<u64, usize>,
    entries: Vec<Entry>,
    heap: BinaryHeap<(Reverse<u64>, u64, usize)>,
}

impl SpaceSaving {
    pub fn new() -> Self {
        Self {
            map: HashMap::with_capacity(DUP_K),
            entries: Vec::with_capacity(DUP_K),
            heap: BinaryHeap::with_capacity(DUP_K),
        }
    }

    pub fn add(&mut self, key: u64, weight: u64) {
        if let Some(&idx) = self.map.get(&key) {
            let e = &mut self.entries[idx];
            e.count += weight;
            self.heap.push((Reverse(e.count), e.key, idx));
            return;
        }

        if self.entries.len() < DUP_K {
            let idx = self.entries.len();
            self.entries.push(Entry {
                key,
                count: weight,
                error: 0,
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
        };
        self.map.insert(key, min_idx);
        self.heap.push((Reverse(min_count + weight), key, min_idx));
    }

    pub fn merge(&mut self, other: &SpaceSaving) {
        let mut items = other.entries.clone();
        items.sort_by_key(|e| e.key);
        for e in items {
            self.add(e.key, e.count);
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
