#[derive(Clone, Copy, Debug)]
pub struct ReadView<'a> {
    pub id: &'a [u8],
    pub seq: &'a [u8],
    pub qual: &'a [u8],
}

#[derive(Clone, Debug)]
pub struct OwnedRead {
    pub id: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}

impl OwnedRead {
    pub fn from_record(record: kira_fastq::FastqRecord<'_>) -> Self {
        Self {
            id: record.header().to_vec(),
            seq: record.seq().to_vec(),
            qual: record.qual().to_vec(),
        }
    }

    pub fn byte_len(&self) -> usize {
        self.id.len() + self.seq.len() + self.qual.len()
    }

    pub fn as_view(&self) -> ReadView<'_> {
        ReadView {
            id: &self.id,
            seq: &self.seq,
            qual: &self.qual,
        }
    }
}
