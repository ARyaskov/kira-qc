use memchr::memchr;
use std::io::BufRead;

#[derive(Clone, Copy, Debug)]
pub struct ReadView<'a> {
    pub id: &'a [u8],
    pub seq: &'a [u8],
    pub qual: &'a [u8],
}

#[derive(Debug)]
pub struct FastqError {
    pub message: String,
    pub byte_offset: usize,
}

impl std::fmt::Display for FastqError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} at offset {}", self.message, self.byte_offset)
    }
}

impl std::error::Error for FastqError {}

fn next_line<'a>(data: &'a [u8], start: usize) -> Option<(&'a [u8], usize)> {
    if start >= data.len() {
        return None;
    }
    let rel = memchr(b'\n', &data[start..])?;
    let end = start + rel;
    let mut line = &data[start..end];
    if line.ends_with(b"\r") {
        line = &line[..line.len() - 1];
    }
    Some((line, end + 1))
}

pub fn parse_chunk<'a>(chunk: &'a [u8], reads: &mut Vec<ReadView<'a>>) -> Result<(), FastqError> {
    reads.clear();
    let mut pos: usize = 0;
    while pos < chunk.len() {
        let record_start = pos;
        let (id, p1) = next_line(chunk, pos).ok_or_else(|| FastqError {
            message: "incomplete FASTQ record".to_string(),
            byte_offset: record_start,
        })?;
        let (seq, p2) = next_line(chunk, p1).ok_or_else(|| FastqError {
            message: "incomplete FASTQ record".to_string(),
            byte_offset: record_start,
        })?;
        let (plus, p3) = next_line(chunk, p2).ok_or_else(|| FastqError {
            message: "incomplete FASTQ record".to_string(),
            byte_offset: record_start,
        })?;
        let (qual, p4) = next_line(chunk, p3).ok_or_else(|| FastqError {
            message: "incomplete FASTQ record".to_string(),
            byte_offset: record_start,
        })?;

        if !id.starts_with(b"@") {
            return Err(FastqError {
                message: "invalid FASTQ id line".to_string(),
                byte_offset: record_start,
            });
        }
        if !plus.starts_with(b"+") {
            return Err(FastqError {
                message: "invalid FASTQ plus line".to_string(),
                byte_offset: record_start,
            });
        }
        if seq.len() != qual.len() {
            return Err(FastqError {
                message: "sequence/quality length mismatch".to_string(),
                byte_offset: record_start,
            });
        }

        reads.push(ReadView { id, seq, qual });
        pos = p4;
    }
    Ok(())
}

pub fn scan_qual_range(data: &[u8], max_reads: usize) -> Result<(u8, u8, usize), FastqError> {
    let mut pos: usize = 0;
    let mut reads: usize = 0;
    let mut min_q: u8 = u8::MAX;
    let mut max_q: u8 = 0;

    while pos < data.len() && reads < max_reads {
        let record_start = pos;
        let (_, p1) = next_line(data, pos).ok_or_else(|| FastqError {
            message: "incomplete FASTQ record".to_string(),
            byte_offset: record_start,
        })?;
        let (_, p2) = next_line(data, p1).ok_or_else(|| FastqError {
            message: "incomplete FASTQ record".to_string(),
            byte_offset: record_start,
        })?;
        let (_, p3) = next_line(data, p2).ok_or_else(|| FastqError {
            message: "incomplete FASTQ record".to_string(),
            byte_offset: record_start,
        })?;
        let (qual, p4) = next_line(data, p3).ok_or_else(|| FastqError {
            message: "incomplete FASTQ record".to_string(),
            byte_offset: record_start,
        })?;
        for &b in qual {
            if b < min_q {
                min_q = b;
            }
            if b > max_q {
                max_q = b;
            }
        }
        pos = p4;
        reads += 1;
    }

    if reads == 0 {
        return Err(FastqError {
            message: "no reads found".to_string(),
            byte_offset: 0,
        });
    }
    Ok((min_q, max_q, reads))
}

pub fn scan_qual_range_reader<R: BufRead>(
    mut reader: R,
    max_reads: usize,
) -> Result<(u8, u8, usize), FastqError> {
    let mut reads: usize = 0;
    let mut min_q: u8 = u8::MAX;
    let mut max_q: u8 = 0;
    let mut line = Vec::new();

    while reads < max_reads {
        line.clear();
        if reader
            .read_until(b'\n', &mut line)
            .map_err(|e| FastqError {
                message: format!("read error: {e}"),
                byte_offset: 0,
            })?
            == 0
        {
            break;
        }
        line.clear();
        if reader
            .read_until(b'\n', &mut line)
            .map_err(|e| FastqError {
                message: format!("read error: {e}"),
                byte_offset: 0,
            })?
            == 0
        {
            return Err(FastqError {
                message: "incomplete FASTQ record".to_string(),
                byte_offset: 0,
            });
        }
        line.clear();
        if reader
            .read_until(b'\n', &mut line)
            .map_err(|e| FastqError {
                message: format!("read error: {e}"),
                byte_offset: 0,
            })?
            == 0
        {
            return Err(FastqError {
                message: "incomplete FASTQ record".to_string(),
                byte_offset: 0,
            });
        }
        line.clear();
        if reader
            .read_until(b'\n', &mut line)
            .map_err(|e| FastqError {
                message: format!("read error: {e}"),
                byte_offset: 0,
            })?
            == 0
        {
            return Err(FastqError {
                message: "incomplete FASTQ record".to_string(),
                byte_offset: 0,
            });
        }
        if line.ends_with(b"\n") {
            line.pop();
            if line.ends_with(b"\r") {
                line.pop();
            }
        }
        for &b in &line {
            if b < min_q {
                min_q = b;
            }
            if b > max_q {
                max_q = b;
            }
        }
        reads += 1;
    }

    if reads == 0 {
        return Err(FastqError {
            message: "no reads found".to_string(),
            byte_offset: 0,
        });
    }
    Ok((min_q, max_q, reads))
}
