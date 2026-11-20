use crate::core::fastq::ReadView;
use crate::core::metrics::per_tile_qual::parse_tile;
use crate::core::metrics::{AdapterSet, Agg, AggConfig, ContaminantList, UpdateTimings};
use crate::core::model::{Encoding, FinalizeContext, Mode};
use crate::core::stats::{StatsLog, fmt_dur};
use anyhow::{Context, Result, anyhow, bail};
use crossbeam_channel as channel;
use kira_fastq::FastqReader;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::thread;
use std::time::{Duration, Instant};

/// Quality bytes examined before committing to an encoding.
const AUTO_DETECT_READS: usize = 100_000;

/// Chunk payload target. Smaller chunks keep peak memory proportional to the
/// thread count instead of exploding with it.
const TARGET_CHUNK_BYTES: usize = 2 * 1024 * 1024;

/// Chunks allowed to sit in the channel between producer and workers.
const CHANNEL_SLACK: usize = 4;

/// Sentinel for "this read had no parseable tile".
const NO_TILE: u32 = u32::MAX;

#[derive(Clone, Copy, Debug)]
pub enum PhredOffsetConfig {
    Auto,
    Fixed(u8),
}

#[derive(Clone, Debug)]
pub enum Input {
    Path(PathBuf),
    Stdin,
}

impl Input {
    pub fn display_name(&self) -> String {
        match self {
            Input::Path(p) => p
                .file_name()
                .and_then(|s| s.to_str())
                .unwrap_or("input")
                .to_string(),
            Input::Stdin => "stdin".to_string(),
        }
    }
}

pub struct RunConfig {
    pub input: Input,
    pub sample_name: String,
    pub threads: usize,
    pub phred_offset: PhredOffsetConfig,
    pub mode: Mode,
    pub kmer_enabled: bool,
    pub dup_window: u64,
    pub nogroup: bool,
    pub adapters: Arc<AdapterSet>,
    pub contaminants: Arc<ContaminantList>,
    pub ignored_modules: Vec<String>,
}

pub struct RunOutput {
    pub agg: Agg,
    pub ctx: FinalizeContext,
}

/// Packed buffer of records: `seq||qual` per record, indexed by `offsets`.
struct PackedChunk {
    bytes: Vec<u8>,
    offsets: Vec<RecordOffset>,
    tiles: Vec<u32>,
    /// Global index of the first read in this chunk.
    first_read_index: u64,
    /// Determined once from the first chunk and carried along, so workers never
    /// need to synchronise and the input is never read twice.
    phred_offset: u8,
}

impl PackedChunk {
    fn empty(first_read_index: u64, phred_offset: u8) -> Self {
        Self {
            bytes: Vec::with_capacity(TARGET_CHUNK_BYTES + 4096),
            offsets: Vec::with_capacity(TARGET_CHUNK_BYTES / 256),
            tiles: Vec::new(),
            first_read_index,
            phred_offset,
        }
    }
}

#[derive(Clone, Copy)]
struct RecordOffset {
    seq_start: u32,
    len: u32,
}

#[derive(Clone, Debug, Default)]
struct ProducerStats {
    chunks: u64,
    bytes: u64,
    reads: u64,
    parse: Duration,
}

#[derive(Clone, Debug, Default)]
struct WorkerStats {
    chunks: u64,
    bytes: u64,
    reads: u64,
    process: Duration,
    kmer: Duration,
    kmer_updates: u64,
}

pub fn run(cfg: RunConfig) -> Result<RunOutput> {
    let log = StatsLog::from_env();
    let t_total = Instant::now();

    let file_name = match &cfg.input {
        Input::Path(p) => p
            .file_name()
            .and_then(|s| s.to_str())
            .map(str::to_string)
            .context("failed to determine input filename")?,
        Input::Stdin => "stdin".to_string(),
    };

    let threads = cfg.threads.max(1);
    let (chunk_tx, chunk_rx) = channel::bounded::<PackedChunk>(CHANNEL_SLACK);
    let (result_tx, result_rx) = channel::bounded::<(Agg, WorkerStats)>(threads);
    let (err_tx, err_rx) = channel::bounded::<anyhow::Error>(1);
    let (meta_tx, meta_rx) = channel::bounded::<ProducerMeta>(1);

    let producer_input = cfg.input.clone();
    let producer_err = err_tx.clone();
    let phred_cfg = cfg.phred_offset;
    let mode = cfg.mode;
    let producer = thread::spawn(move || {
        run_producer(
            producer_input,
            phred_cfg,
            mode,
            chunk_tx,
            producer_err,
            meta_tx,
        );
    });

    let agg_cfg = AggConfig {
        mode: cfg.mode,
        kmer_enabled: cfg.kmer_enabled && cfg.mode == Mode::Short,
        dup_window: cfg.dup_window,
        per_tile_enabled: cfg.mode == Mode::Short,
        adapters: cfg.adapters.clone(),
        contaminants: cfg.contaminants.clone(),
    };
    let kmer_enabled = agg_cfg.kmer_enabled;

    let mut workers = Vec::with_capacity(threads);
    let t_workers = Instant::now();
    for _ in 0..threads {
        let rx = chunk_rx.clone();
        let tx = result_tx.clone();
        let stats_enabled = log.enabled();
        let worker_cfg = agg_cfg.clone();
        workers.push(thread::spawn(move || {
            run_worker(rx, tx, worker_cfg, stats_enabled);
        }));
    }
    drop(chunk_rx);
    drop(result_tx);
    drop(err_tx);
    log.stage("engine.spawn_workers", t_workers);

    let t_collect = Instant::now();
    let mut final_agg = Agg::new(agg_cfg.clone());
    let mut worker_stats_total = WorkerStats::default();
    let mut wait_time = Duration::ZERO;
    let mut merge_time = Duration::ZERO;
    let mut first_error: Option<anyhow::Error> = None;

    let mut err_open = true;
    let mut results_open = true;
    while results_open {
        let t_wait = Instant::now();
        let msg = if err_open {
            channel::select! {
                recv(err_rx) -> err => match err {
                    Ok(err) => {
                        // Record the failure but keep draining so the worker
                        // threads can finish and be joined cleanly.
                        if first_error.is_none() {
                            first_error = Some(err);
                        }
                        err_open = false;
                        continue;
                    }
                    Err(_) => {
                        err_open = false;
                        continue;
                    }
                },
                recv(result_rx) -> msg => msg,
            }
        } else {
            result_rx.recv()
        };
        wait_time += t_wait.elapsed();

        match msg {
            Ok((agg, ws)) => {
                let tm = Instant::now();
                final_agg.merge(&agg);
                merge_time += tm.elapsed();
                accumulate_worker_stats(&mut worker_stats_total, &ws);
            }
            Err(_) => {
                results_open = false;
            }
        }
    }
    log.stage("engine.collect_merge", t_collect);

    // A panicking worker used to be silently dropped, producing a report with
    // missing reads and no error.
    if producer.join().is_err() {
        bail!("FASTQ producer thread panicked");
    }
    let mut panicked = 0usize;
    for w in workers {
        if w.join().is_err() {
            panicked += 1;
        }
    }
    if panicked > 0 {
        bail!("{panicked} worker thread(s) panicked; results would be incomplete");
    }
    if let Some(err) = first_error {
        return Err(err);
    }

    let meta = meta_rx.recv().unwrap_or_default();
    if final_agg.total_reads == 0 {
        bail!("no FASTQ records found in input");
    }

    let encoding = if meta.phred_offset == 64 {
        Encoding::Illumina15
    } else {
        Encoding::Sanger
    };

    if log.enabled() {
        emit_stats(&meta.producer, &worker_stats_total, wait_time, merge_time);
    }
    log.stage("engine.total", t_total);

    let dup_window_complete = cfg.dup_window == 0 || final_agg.total_reads <= cfg.dup_window;

    let ctx = FinalizeContext {
        encoding,
        file_name,
        sample_name: cfg.sample_name.clone(),
        mode: cfg.mode,
        nogroup: cfg.nogroup,
        kmer_enabled,
        dup_window_reads: final_agg.seq_tracker.observed_reads(),
        dup_window_complete,
        ignored_modules: cfg.ignored_modules.clone(),
    };

    Ok(RunOutput {
        agg: final_agg,
        ctx,
    })
}

#[derive(Clone, Debug, Default)]
struct ProducerMeta {
    producer: ProducerStats,
    phred_offset: u8,
}

fn open_reader(input: &Input) -> Result<FastqReader> {
    match input {
        Input::Path(path) => FastqReader::from_path_auto(path)
            .map_err(|e| anyhow!("failed to open FASTQ input {}: {e:?}", path.display())),
        Input::Stdin => Ok(FastqReader::from_reader(BufReader::with_capacity(
            1 << 20,
            std::io::stdin(),
        ))),
    }
}

fn run_producer(
    input: Input,
    phred_cfg: PhredOffsetConfig,
    mode: Mode,
    chunk_tx: channel::Sender<PackedChunk>,
    err_tx: channel::Sender<anyhow::Error>,
    meta_tx: channel::Sender<ProducerMeta>,
) {
    let mut reader = match open_reader(&input) {
        Ok(reader) => reader,
        Err(e) => {
            let _ = err_tx.send(e);
            let _ = meta_tx.send(ProducerMeta::default());
            return;
        }
    };

    let mut stats = ProducerStats::default();

    // Encoding detection runs on the records of the first chunk instead of
    // re-opening the file. Re-opening meant a gzip input was inflated twice.
    let fixed_offset = match phred_cfg {
        PhredOffsetConfig::Fixed(v) => Some(v),
        PhredOffsetConfig::Auto => None,
    };
    let mut detect = PhredDetector::new(mode);
    let mut phred_offset = fixed_offset.unwrap_or(33);
    let mut offset_decided = fixed_offset.is_some();

    // Tile parsing stops as soon as the first chunk shows there are no tiles.
    let mut parse_tiles = mode == Mode::Short;
    let mut tiles_seen = false;

    let mut read_index: u64 = 0;
    let mut current = PackedChunk::empty(0, phred_offset);
    let mut pending: Vec<PackedChunk> = Vec::new();

    loop {
        let t_next = Instant::now();
        let rec = match reader.next() {
            Ok(Some(rec)) => rec,
            Ok(None) => break,
            Err(e) => {
                let _ = err_tx.send(anyhow!("FASTQ parse/read error: {e:?}"));
                let _ = meta_tx.send(ProducerMeta {
                    producer: stats,
                    phred_offset,
                });
                return;
            }
        };
        stats.parse += t_next.elapsed();

        let seq = rec.seq();
        let qual = rec.qual();
        if seq.len() != qual.len() {
            let _ = err_tx.send(anyhow!(
                "seq/qual length mismatch at read {}: seq={} qual={}",
                read_index + 1,
                seq.len(),
                qual.len()
            ));
            let _ = meta_tx.send(ProducerMeta {
                producer: stats,
                phred_offset,
            });
            return;
        }

        if !offset_decided {
            detect.observe(qual);
        }

        let tile = if parse_tiles {
            match parse_tile(rec.header()) {
                Some(t) => {
                    tiles_seen = true;
                    t
                }
                None => NO_TILE,
            }
        } else {
            NO_TILE
        };

        let len = seq.len();
        let seq_start = current.bytes.len() as u32;
        current.bytes.extend_from_slice(seq);
        current.bytes.extend_from_slice(qual);
        current.offsets.push(RecordOffset {
            seq_start,
            len: len as u32,
        });
        current.tiles.push(tile);
        read_index += 1;

        if current.bytes.len() >= TARGET_CHUNK_BYTES {
            if !offset_decided && detect.is_confident() {
                phred_offset = detect.decide();
                offset_decided = true;
            }
            let next_index = read_index;
            let mut chunk =
                std::mem::replace(&mut current, PackedChunk::empty(next_index, phred_offset));
            stats.chunks += 1;
            stats.reads += chunk.offsets.len() as u64;
            stats.bytes += chunk.bytes.len() as u64;

            if offset_decided {
                for mut queued in pending.drain(..) {
                    queued.phred_offset = phred_offset;
                    if chunk_tx.send(queued).is_err() {
                        let _ = meta_tx.send(ProducerMeta {
                            producer: stats,
                            phred_offset,
                        });
                        return;
                    }
                }
                chunk.phred_offset = phred_offset;
                if chunk_tx.send(chunk).is_err() {
                    let _ = meta_tx.send(ProducerMeta {
                        producer: stats,
                        phred_offset,
                    });
                    return;
                }
            } else {
                // Hold the chunk back until enough quality bytes have been seen
                // to pick an encoding.
                pending.push(chunk);
            }

            if parse_tiles && !tiles_seen {
                parse_tiles = false;
            }
        }
    }

    if !offset_decided {
        phred_offset = detect.decide();
    }

    if !current.offsets.is_empty() {
        stats.chunks += 1;
        stats.reads += current.offsets.len() as u64;
        stats.bytes += current.bytes.len() as u64;
        pending.push(current);
    }

    for mut queued in pending.drain(..) {
        queued.phred_offset = phred_offset;
        if chunk_tx.send(queued).is_err() {
            break;
        }
    }

    let _ = meta_tx.send(ProducerMeta {
        producer: stats,
        phred_offset,
    });
}

/// Quality-encoding detector.
///
/// The old rule mapped "no quality byte below 64" straight to Phred+64, which
/// misreads any high-quality modern run - binned NovaSeq output and PacBio HiFi
/// in particular - and shifts every score by 31. Detection now requires
/// positive evidence for Phred+64 and defaults to Phred+33 otherwise.
struct PhredDetector {
    min_q: u8,
    max_q: u8,
    reads: usize,
    mode: Mode,
}

impl PhredDetector {
    fn new(mode: Mode) -> Self {
        Self {
            min_q: u8::MAX,
            max_q: 0,
            reads: 0,
            mode,
        }
    }

    fn observe(&mut self, qual: &[u8]) {
        for &b in qual {
            if b < self.min_q {
                self.min_q = b;
            }
            if b > self.max_q {
                self.max_q = b;
            }
        }
        self.reads += 1;
    }

    fn is_confident(&self) -> bool {
        // A single byte below the Phred+64 range settles the question.
        self.min_q < 59 || self.reads >= AUTO_DETECT_READS
    }

    fn decide(&self) -> u8 {
        decide_phred_offset(self.min_q, self.max_q, self.reads, self.mode)
    }
}

pub fn decide_phred_offset(min_q: u8, max_q: u8, reads: usize, mode: Mode) -> u8 {
    if reads == 0 {
        return 33;
    }
    // Anything below '(' cannot be Phred+64.
    if min_q < 59 {
        return 33;
    }
    // Phred+64 tops out at 'h' (104). Higher bytes mean Phred+33 with a high
    // score, which is normal for PacBio HiFi and for Q40+ Illumina data.
    if max_q > 104 {
        return 33;
    }
    // Long-read platforms never used Phred+64.
    if mode == Mode::Long {
        return 33;
    }
    // Require a score above 'J' (74): under Phred+33 that would be Q42+, which
    // no Illumina 1.8+ pipeline emits, so Phred+64 is the better reading.
    if min_q >= 64 && max_q > 74 {
        return 64;
    }
    33
}

fn run_worker(
    rx: channel::Receiver<PackedChunk>,
    tx: channel::Sender<(Agg, WorkerStats)>,
    cfg: AggConfig,
    stats_enabled: bool,
) {
    let mut agg = Agg::new(cfg);
    let mut wstats = WorkerStats::default();

    for chunk in rx.iter() {
        let t_start = Instant::now();
        wstats.chunks += 1;
        wstats.bytes += chunk.bytes.len() as u64;
        wstats.reads += chunk.offsets.len() as u64;
        let phred_offset = chunk.phred_offset;

        if stats_enabled {
            let mut timing = UpdateTimings::default();
            for (i, off) in chunk.offsets.iter().enumerate() {
                let view = view_from_offset(&chunk, *off, i);
                let index = chunk.first_read_index + i as u64;
                agg.update_read_timed(&view, phred_offset, index, &mut timing);
            }
            wstats.kmer += timing.kmer;
            wstats.kmer_updates += timing.kmer_updates;
        } else {
            for (i, off) in chunk.offsets.iter().enumerate() {
                let view = view_from_offset(&chunk, *off, i);
                let index = chunk.first_read_index + i as u64;
                agg.update_read(&view, phred_offset, index);
            }
        }

        wstats.process += t_start.elapsed();
    }

    let _ = tx.send((agg, wstats));
}

#[inline]
fn view_from_offset(chunk: &PackedChunk, off: RecordOffset, i: usize) -> ReadView<'_> {
    let start = off.seq_start as usize;
    let len = off.len as usize;
    let seq = &chunk.bytes[start..start + len];
    let qual = &chunk.bytes[start + len..start + 2 * len];
    let tile = match chunk.tiles.get(i) {
        Some(&t) if t != NO_TILE => Some(t),
        _ => None,
    };
    ReadView { seq, qual, tile }
}

fn accumulate_worker_stats(total: &mut WorkerStats, ws: &WorkerStats) {
    total.chunks += ws.chunks;
    total.bytes += ws.bytes;
    total.reads += ws.reads;
    total.process += ws.process;
    total.kmer += ws.kmer;
    total.kmer_updates += ws.kmer_updates;
}

fn emit_stats(
    prod: &ProducerStats,
    workers: &WorkerStats,
    wait_time: Duration,
    merge_time: Duration,
) {
    if prod.chunks > 0 {
        let avg = prod.bytes as f64 / prod.chunks as f64;
        eprintln!(
            "KIRA_STATS producer.chunks={} producer.avg_chunk_bytes={:.0} producer.bytes={} producer.reads={} producer.parse={}",
            prod.chunks,
            avg,
            prod.bytes,
            prod.reads,
            fmt_dur(prod.parse),
        );
    }
    eprintln!(
        "KIRA_STATS worker.chunks={} worker.bytes={} worker.reads={} worker.process={} worker.kmer={} worker.kmer_updates={}",
        workers.chunks,
        workers.bytes,
        workers.reads,
        fmt_dur(workers.process),
        fmt_dur(workers.kmer),
        workers.kmer_updates,
    );
    eprintln!(
        "KIRA_STATS reducer.wait={} reducer.merge_cost={}",
        fmt_dur(wait_time),
        fmt_dur(merge_time)
    );
}

/// Quick look at the first records of a file, used to warn about an obviously
/// wrong `--mode`.
pub fn peek_max_read_length(path: &Path, reads: usize) -> Option<usize> {
    let mut reader = FastqReader::from_path_auto(path).ok()?;
    let mut max = 0usize;
    for _ in 0..reads {
        match reader.next() {
            Ok(Some(rec)) => max = max.max(rec.seq().len()),
            _ => break,
        }
    }
    Some(max)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn low_quality_bytes_mean_phred33() {
        assert_eq!(decide_phred_offset(35, 74, 1000, Mode::Short), 33);
    }

    #[test]
    fn classic_illumina_15_is_detected() {
        // Phred+64 data: nothing below '@' (64), scores up to 'h'.
        assert_eq!(decide_phred_offset(64, 104, 1000, Mode::Short), 64);
    }

    #[test]
    fn high_quality_phred33_is_not_mistaken_for_phred64() {
        // Binned NovaSeq output: no byte below 'F' (70), max 'F'. Under the old
        // rule min_q >= 64 forced Phred+64 and shifted every score by 31.
        assert_eq!(decide_phred_offset(70, 70, 1000, Mode::Short), 33);
    }

    #[test]
    fn pacbio_hifi_stays_phred33() {
        // HiFi reaches Q60+, i.e. bytes well past the Phred+64 ceiling.
        assert_eq!(decide_phred_offset(70, 126, 1000, Mode::Short), 33);
    }

    #[test]
    fn long_mode_never_selects_phred64() {
        assert_eq!(decide_phred_offset(64, 104, 1000, Mode::Long), 33);
    }

    #[test]
    fn empty_input_defaults_to_phred33() {
        assert_eq!(decide_phred_offset(u8::MAX, 0, 0, Mode::Short), 33);
    }

    #[test]
    fn detector_short_circuits_on_low_quality() {
        let mut d = PhredDetector::new(Mode::Short);
        d.observe(b"!!!!");
        assert!(d.is_confident());
        assert_eq!(d.decide(), 33);
    }
}
