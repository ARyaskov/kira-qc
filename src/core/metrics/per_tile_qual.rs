//! Per tile sequence quality.
//!
//! Reports, for every flowcell tile and base position, how far the mean quality
//! deviates from the mean across all tiles. Localised drops point at a physical
//! problem on the flowcell (bubbles, smudges, a bad edge) rather than at the
//! library.

use crate::core::model::{BaseGroup, Status};
use rustc_hash::FxHashMap;

/// Beyond this many distinct tiles the module is dropped rather than allowed to
/// grow without bound; patterned flowcells stay well below it.
pub const MAX_TILES: usize = 4096;

#[derive(Clone, Debug, Default)]
struct TileAcc {
    sum: Vec<u64>,
    count: Vec<u64>,
}

impl TileAcc {
    fn ensure(&mut self, len: usize) {
        if self.sum.len() < len {
            self.sum.resize(len, 0);
            self.count.resize(len, 0);
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct PerTileQual {
    tiles: FxHashMap<u32, TileAcc>,
    /// Set when the tile budget was exceeded; the module is then not reported.
    overflowed: bool,
    /// Set once any read carried a parseable tile.
    saw_tile: bool,
}

#[derive(Clone, Debug)]
pub struct PerTileRow {
    pub tile: u32,
    pub group: BaseGroup,
    pub mean_deviation: f64,
}

impl PerTileQual {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn is_available(&self) -> bool {
        self.saw_tile && !self.overflowed && !self.tiles.is_empty()
    }

    pub fn overflowed(&self) -> bool {
        self.overflowed
    }

    pub fn tile_count(&self) -> usize {
        self.tiles.len()
    }

    /// Accumulate one read's qualities under its tile.
    pub fn add(&mut self, tile: u32, qual: &[u8], phred_offset: u8) {
        if self.overflowed || qual.is_empty() {
            return;
        }
        self.saw_tile = true;
        if !self.tiles.contains_key(&tile) && self.tiles.len() >= MAX_TILES {
            self.overflowed = true;
            self.tiles.clear();
            return;
        }
        let acc = self.tiles.entry(tile).or_default();
        acc.ensure(qual.len());
        for (i, &q) in qual.iter().enumerate() {
            acc.sum[i] += q.saturating_sub(phred_offset) as u64;
            acc.count[i] += 1;
        }
    }

    pub fn merge(&mut self, other: &PerTileQual) {
        if other.saw_tile {
            self.saw_tile = true;
        }
        if self.overflowed || other.overflowed {
            self.overflowed = true;
            self.tiles.clear();
            return;
        }
        for (&tile, acc) in &other.tiles {
            if !self.tiles.contains_key(&tile) && self.tiles.len() >= MAX_TILES {
                self.overflowed = true;
                self.tiles.clear();
                return;
            }
            let target = self.tiles.entry(tile).or_default();
            target.ensure(acc.sum.len());
            for i in 0..acc.sum.len() {
                target.sum[i] += acc.sum[i];
                target.count[i] += acc.count[i];
            }
        }
    }

    /// Deviation of each tile from the all-tile mean, collapsed into base groups.
    pub fn rows(&self, groups: &[BaseGroup]) -> (Vec<PerTileRow>, Status) {
        let mut rows = Vec::new();
        let mut status = Status::Pass;
        if !self.is_available() || groups.is_empty() {
            return (rows, status);
        }

        let mut tile_ids: Vec<u32> = self.tiles.keys().copied().collect();
        tile_ids.sort_unstable();
        let mut group_means = vec![0.0; groups.len()];
        for tile in tile_ids {
            let acc = &self.tiles[&tile];
            for (i, &group) in groups.iter().enumerate() {
                let mut sum = 0.0;
                let mut covered = 0;
                for position in group.range().take_while(|&pos| pos < acc.sum.len()) {
                    let count = acc.count[position];
                    if count > 0 {
                        sum += acc.sum[position] as f64 / count as f64;
                        covered += 1;
                    }
                }
                let mean = if covered == 0 {
                    0.0
                } else {
                    sum / covered as f64
                };
                group_means[i] += mean;
                rows.push(PerTileRow {
                    tile,
                    group,
                    mean_deviation: mean,
                });
            }
        }
        for (row, total) in rows.iter_mut().zip(group_means.iter().cycle()) {
            row.mean_deviation -= total / self.tiles.len() as f64;
            if row.mean_deviation.abs() > 10.0 {
                status = Status::Fail;
            } else if row.mean_deviation.abs() > 5.0 && status != Status::Fail {
                status = Status::Warn;
            }
        }

        (rows, status)
    }
}

/// Extract the tile number from an Illumina read header.
///
/// Casava 1.8+: `@instrument:run:flowcell:lane:tile:x:y ...` (tile is field 4).
/// Older Casava: `@instrument:lane:tile:x:y#index/1` (tile is field 2).
/// Anything else (ONT, PacBio, SRA-renamed reads) has no tile.
pub fn parse_tile(header: &[u8]) -> Option<u32> {
    // Stop at the first whitespace: the description after it may contain colons.
    let end = header
        .iter()
        .position(|&b| b == b' ' || b == b'\t')
        .unwrap_or(header.len());
    let name = &header[..end];
    let name = name.strip_prefix(b"@").unwrap_or(name);

    let mut fields = name.split(|&b| b == b':');
    let legacy_tile = fields.nth(2)?;
    let modern_tile = fields.nth(1)?;
    let field = match (fields.next(), fields.next()) {
        (None, None) => legacy_tile,
        (Some(_), Some(_)) => modern_tile,
        _ => return None,
    };

    if field.is_empty() || !field.iter().all(|b| b.is_ascii_digit()) {
        return None;
    }
    std::str::from_utf8(field).ok()?.parse::<u32>().ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_casava_18_header() {
        let h = b"@M01234:55:000000000-A1B2C:1:1101:15589:1338 1:N:0:ATCACG";
        assert_eq!(parse_tile(h), Some(1101));
    }

    #[test]
    fn parses_legacy_header() {
        let h = b"@HWUSI-EAS100R:6:73:941:1973#0/1";
        assert_eq!(parse_tile(h), Some(73));
    }

    #[test]
    fn rejects_nanopore_header() {
        let h =
            b"@a9ef233c-5db5-4199-9623-63c7dc060a94 runid=d841 sampleid=rbclA read=19309 ch=104";
        assert_eq!(parse_tile(h), None);
    }

    #[test]
    fn rejects_non_numeric_tile() {
        let h = b"@a:b:c:d:tile:f:g";
        assert_eq!(parse_tile(h), None);
    }

    #[test]
    fn parses_extended_headers_and_rejects_six_field_names() {
        assert_eq!(parse_tile(b"@a:b:c:d:101:42:43:extra"), Some(101));
        assert_eq!(parse_tile(b"@a:b:100:d:101:42"), None);
        assert_eq!(parse_tile(b"@a:b:c:d:4294967296:42:43"), None);
    }

    #[test]
    fn deviation_is_relative_to_all_tiles() {
        let mut p = PerTileQual::new();
        // tile 1 is 10 points worse than tile 2 at every position
        p.add(1, &[b'!' + 20, b'!' + 20], 33);
        p.add(2, &[b'!' + 30, b'!' + 30], 33);
        let groups = vec![
            BaseGroup { start: 1, end: 1 },
            BaseGroup { start: 2, end: 2 },
        ];
        let (rows, status) = p.rows(&groups);
        assert_eq!(rows.len(), 4);
        assert!((rows[0].mean_deviation + 5.0).abs() < 1e-9);
        assert!((rows[2].mean_deviation - 5.0).abs() < 1e-9);
        assert_eq!(status, Status::Pass);
        p.add(1, &[b'!' + 14, b'!' + 14], 33);
        assert_eq!(p.rows(&groups).1, Status::Warn);
    }

    #[test]
    fn a_badly_lagging_tile_fails() {
        let mut p = PerTileQual::new();
        p.add(1, &[b'!' + 10], 33);
        p.add(2, &[b'!' + 35], 33);
        p.add(2, &[b'!' + 35], 33);
        let groups = vec![BaseGroup { start: 1, end: 1 }];
        let (_, status) = p.rows(&groups);
        assert_eq!(status, Status::Fail);
    }

    #[test]
    fn merge_is_order_independent() {
        let mut a = PerTileQual::new();
        a.add(1, &[b'!' + 30], 33);
        let mut b = PerTileQual::new();
        b.add(2, &[b'!' + 30], 33);

        let groups = vec![BaseGroup { start: 1, end: 1 }];
        let mut ab = a.clone();
        ab.merge(&b);
        let mut ba = b.clone();
        ba.merge(&a);
        let ra = ab.rows(&groups).0;
        let rb = ba.rows(&groups).0;
        assert_eq!(ra.len(), rb.len());
        for (x, y) in ra.iter().zip(rb.iter()) {
            assert_eq!(x.tile, y.tile);
            assert!((x.mean_deviation - y.mean_deviation).abs() < 1e-12);
        }
    }

    #[test]
    fn module_disabled_without_tiles() {
        let p = PerTileQual::new();
        assert!(!p.is_available());
    }
}
