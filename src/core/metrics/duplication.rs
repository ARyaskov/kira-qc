use rustc_hash::FxHashMap;
use std::cmp::Reverse;
use std::collections::BinaryHeap;

const DUP_K: usize = 200_000;
pub const HEAP_COMPACT_THRESHOLD: usize = 4 * DUP_K;

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
}

#[derive(Clone, Debug)]
pub struct SpaceSaving {
    map: FxHashMap<u64, u32>,
    entries: Vec<Entry>,
    heap: BinaryHeap<(Reverse<u64>, u64, u32)>,
}

impl Default for SpaceSaving {
    fn default() -> Self {
        Self::new()
    }
}

impl SpaceSaving {
    pub fn new() -> Self {
        Self {
            map: FxHashMap::with_capacity_and_hasher(DUP_K, Default::default()),
            entries: Vec::with_capacity(DUP_K),
            heap: BinaryHeap::with_capacity(DUP_K),
        }
    }

    pub fn add(&mut self, key: u64, weight: u64) {
        if let Some(&idx) = self.map.get(&key) {
            let e = &mut self.entries[idx as usize];
            e.count += weight;
            self.heap.push((Reverse(e.count), e.key, idx));
            self.maybe_compact_heap();
            return;
        }

        if self.entries.len() < DUP_K {
            let idx = self.entries.len() as u32;
            self.entries.push(Entry { key, count: weight });
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
        };
        self.map.insert(key, min_idx);
        self.heap.push((Reverse(new_count), key, min_idx));
        self.maybe_compact_heap();
    }

    pub fn merge(&mut self, other: &SpaceSaving) {
        let mut items: Vec<&Entry> = other.entries.iter().collect();
        items.sort_unstable_by_key(|e| e.key);
        for e in items {
            self.add(e.key, e.count);
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

    /// Rebuild the heap to cap peak memory at O(DUP_K).
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

    pub fn heap_len(&self) -> usize {
        self.heap.len()
    }
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
