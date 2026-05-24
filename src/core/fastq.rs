/// Zero-copy view of one FASTQ record (header is dropped during packing).
#[derive(Clone, Copy, Debug)]
pub struct ReadView<'a> {
    pub seq: &'a [u8],
    pub qual: &'a [u8],
}
