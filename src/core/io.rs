use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use gzp::deflate::{Bgzf, Mgzip};
use gzp::par::decompress::ParDecompressBuilder;
use memmap2::Mmap;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;
use std::sync::Arc;
use std::time::{Duration, Instant};

pub struct MmapSource {
    mmap: Mmap,
}

impl MmapSource {
    pub fn open(path: &Path) -> Result<Self> {
        let file =
            File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
        // SAFETY: read-only file mapping.
        let mmap = unsafe { Mmap::map(&file) }.with_context(|| "mmap failed")?;
        Ok(Self { mmap })
    }

    pub fn bytes(&self) -> &[u8] {
        &self.mmap
    }

    pub fn len(&self) -> usize {
        self.mmap.len()
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum InputKind {
    Plain,
    Gzip,
}

#[derive(Clone, Debug)]
pub enum ChunkData {
    MmapRange { start: usize, end: usize },
    Owned(Vec<u8>),
}

#[derive(Clone, Debug)]
pub struct Chunk {
    pub index: usize,
    pub data: ChunkData,
    pub timing: ChunkTiming,
}

pub const CHUNK_SIZE: usize = 16 * 1024 * 1024;
const GZIP_READ_BUF: usize = 8 * 1024 * 1024;

#[derive(Clone, Copy, Debug, Default)]
pub struct ChunkTiming {
    pub bytes: usize,
    pub decompress: Duration,
    pub align: Duration,
}

pub struct MmapChunker {
    data: Arc<MmapSource>,
    pos: usize,
    chunk_size: usize,
    index: usize,
}

impl MmapChunker {
    pub fn new(data: Arc<MmapSource>, chunk_size: usize) -> Self {
        Self {
            data,
            pos: 0,
            chunk_size,
            index: 0,
        }
    }

    pub fn next_chunk(&mut self) -> Option<Chunk> {
        let bytes = self.data.bytes();
        let len = bytes.len();
        if self.pos >= len {
            return None;
        }
        let t_align = Instant::now();
        let start = self.pos;
        let mut end = (start + self.chunk_size).min(len);
        let mut i = start;
        let mut lines: u32 = 0;
        while i < len {
            if bytes[i] == b'\n' {
                lines += 1;
                if lines > 0 && i + 1 >= end && lines % 4 == 0 {
                    i += 1;
                    break;
                }
            }
            i += 1;
            if i >= end && lines > 0 && lines % 4 == 0 {
                break;
            }
        }
        end = i.min(len);
        self.pos = end;
        let chunk = Chunk {
            index: self.index,
            data: ChunkData::MmapRange { start, end },
            timing: ChunkTiming {
                bytes: end - start,
                decompress: Duration::ZERO,
                align: t_align.elapsed(),
            },
        };
        self.index += 1;
        Some(chunk)
    }
}

pub struct GzipChunker {
    decoder: Box<dyn Read + Send>,
    buffer: Vec<u8>,
    read_buf: Vec<u8>,
    chunk_size: usize,
    index: usize,
    eof: bool,
    scan_pos: usize,
    nl_count: u32,
    last_cut: usize,
    total_out: usize,
    acc_decompress: Duration,
    acc_align: Duration,
}

impl GzipChunker {
    pub fn open(path: &Path, chunk_size: usize, threads: usize) -> Result<Self> {
        let decoder = open_gzip_reader(path, threads)?;
        Ok(Self {
            decoder,
            buffer: Vec::with_capacity(chunk_size + (chunk_size / 4)),
            read_buf: vec![0u8; GZIP_READ_BUF],
            chunk_size,
            index: 0,
            eof: false,
            scan_pos: 0,
            nl_count: 0,
            last_cut: 0,
            total_out: 0,
            acc_decompress: Duration::ZERO,
            acc_align: Duration::ZERO,
        })
    }

    pub fn next_chunk(&mut self) -> Result<Option<Chunk>> {
        loop {
            if (self.buffer.len() >= self.chunk_size || self.eof) && self.last_cut > 0 {
                let tail = self.buffer.split_off(self.last_cut);
                let chunk_bytes = std::mem::take(&mut self.buffer);
                self.buffer = tail;
                self.scan_pos = 0;
                self.nl_count = 0;
                self.last_cut = 0;
                let bytes_len = chunk_bytes.len();
                let chunk = Chunk {
                    index: self.index,
                    data: ChunkData::Owned(chunk_bytes),
                    timing: ChunkTiming {
                        bytes: bytes_len,
                        decompress: self.acc_decompress,
                        align: self.acc_align,
                    },
                };
                self.acc_decompress = Duration::ZERO;
                self.acc_align = Duration::ZERO;
                self.index += 1;
                return Ok(Some(chunk));
            }

            if self.eof {
                if self.buffer.is_empty() {
                    return Ok(None);
                }
                return Err(anyhow::anyhow!(
                    "incomplete FASTQ record at gzip offset {}",
                    self.total_out.saturating_sub(self.buffer.len())
                ));
            }

            let t_read = Instant::now();
            let n = self.decoder.read(&mut self.read_buf).with_context(|| {
                format!(
                    "gzip decompression error at chunk {} (offset {})",
                    self.index, self.total_out
                )
            })?;
            self.acc_decompress += t_read.elapsed();
            if n == 0 {
                self.eof = true;
                continue;
            }
            self.buffer.extend_from_slice(&self.read_buf[..n]);
            self.total_out += n;

            let t_align = Instant::now();
            while self.scan_pos < self.buffer.len() {
                if self.buffer[self.scan_pos] == b'\n' {
                    self.nl_count += 1;
                    if self.nl_count % 4 == 0 {
                        self.last_cut = self.scan_pos + 1;
                    }
                }
                self.scan_pos += 1;
            }
            self.acc_align += t_align.elapsed();
        }
    }
}

pub enum InputSource {
    Mmap { chunker: MmapChunker },
    Gzip { chunker: GzipChunker },
}

impl InputSource {
    pub fn open(path: &Path, threads: usize) -> Result<(Self, Option<Arc<MmapSource>>, InputKind)> {
        let kind = detect_input_kind(path)?;
        match kind {
            InputKind::Plain => {
                let source = Arc::new(MmapSource::open(path)?);
                let chunker = MmapChunker::new(Arc::clone(&source), CHUNK_SIZE);
                Ok((InputSource::Mmap { chunker }, Some(source), kind))
            }
            InputKind::Gzip => {
                let chunker = GzipChunker::open(path, CHUNK_SIZE, threads)?;
                Ok((InputSource::Gzip { chunker }, None, kind))
            }
        }
    }

    pub fn next_chunk(&mut self) -> Result<Option<Chunk>> {
        match self {
            InputSource::Mmap { chunker } => Ok(chunker.next_chunk()),
            InputSource::Gzip { chunker } => chunker.next_chunk(),
        }
    }
}

pub fn detect_input_kind(path: &Path) -> Result<InputKind> {
    if let Some(ext) = path.extension().and_then(|s| s.to_str()) {
        let ext = ext.to_ascii_lowercase();
        if ext == "gz" {
            return Ok(InputKind::Gzip);
        }
    }
    let mut file =
        File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let mut magic = [0u8; 2];
    let n = file
        .read(&mut magic)
        .with_context(|| "failed to read magic bytes")?;
    if n == 2 && magic == [0x1f, 0x8b] {
        Ok(InputKind::Gzip)
    } else {
        Ok(InputKind::Plain)
    }
}

#[derive(Clone, Copy, Debug)]
enum GzipVariant {
    Standard,
    Mgzip,
    Bgzf,
}

fn detect_gzip_variant(path: &Path) -> Result<GzipVariant> {
    let mut file =
        File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let mut header = [0u8; 20];
    let n = file
        .read(&mut header)
        .with_context(|| "failed to read gzip header")?;
    if n < 14 {
        return Ok(GzipVariant::Standard);
    }
    if header[0] != 0x1f || header[1] != 0x8b {
        return Ok(GzipVariant::Standard);
    }
    if header[3] & 4 == 0 {
        return Ok(GzipVariant::Standard);
    }
    if header[12] == b'B' && header[13] == b'C' {
        return Ok(GzipVariant::Bgzf);
    }
    if header[12] == b'I' && header[13] == b'G' {
        return Ok(GzipVariant::Mgzip);
    }
    Ok(GzipVariant::Standard)
}

pub fn open_gzip_reader(path: &Path, threads: usize) -> Result<Box<dyn Read + Send>> {
    let variant = detect_gzip_variant(path)?;
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let reader = BufReader::new(file);
    let reader: Box<dyn Read + Send> = match variant {
        GzipVariant::Bgzf => {
            if threads > 1 {
                Box::new(
                    ParDecompressBuilder::<Bgzf>::new()
                        .num_threads(threads)
                        .unwrap()
                        .from_reader(reader),
                )
            } else {
                Box::new(MultiGzDecoder::new(reader))
            }
        }
        GzipVariant::Mgzip => {
            if threads > 1 {
                Box::new(
                    ParDecompressBuilder::<Mgzip>::new()
                        .num_threads(threads)
                        .unwrap()
                        .from_reader(reader),
                )
            } else {
                Box::new(MultiGzDecoder::new(reader))
            }
        }
        GzipVariant::Standard => Box::new(MultiGzDecoder::new(reader)),
    };
    Ok(reader)
}
