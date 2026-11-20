/// Zero-copy view of one FASTQ record.
///
/// The header itself is dropped during packing; the only thing needed
/// downstream is the flowcell tile, which the producer parses once.
#[derive(Clone, Copy, Debug)]
pub struct ReadView<'a> {
    pub seq: &'a [u8],
    pub qual: &'a [u8],
    pub tile: Option<u32>,
}

impl<'a> ReadView<'a> {
    pub fn new(seq: &'a [u8], qual: &'a [u8]) -> Self {
        Self {
            seq,
            qual,
            tile: None,
        }
    }
}
