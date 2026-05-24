use crate::core::fastq::ReadView;
use crate::core::metrics::{Agg, UpdateTimings};
use crate::core::model::{Encoding, FinalizeContext, Mode};
use anyhow::{Context, Result, anyhow, bail};
use crossbeam_channel as channel;
use kira_fastq::FastqReader;
use std::path::PathBuf;
use std::sync::Arc;
use std::thread;
use std::time::{Duration, Instant};

const AUTO_DETECT_READS: usize = 50_000;
const TARGET_CHUNK_BYTES: usize = 16 * 1024 * 1024;

pub enum PhredOffsetConfig {
    Auto,
    Fixed(u8),
}

pub struct RunConfig {
    pub reads1: PathBuf,
    pub sample_name: String,
    pub threads: usize,
    pub phred_offset: PhredOffsetConfig,
    pub mode: Mode,
}

pub struct RunOutput {
    pub agg: Agg,
    pub ctx: FinalizeContext,
}

/// Packed buffer of records: `seq||qual` per record, indexed by `offsets`.
struct PackedChunk {
    bytes: Vec<u8>,
    offsets: Vec<RecordOffset>,
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
    let stats = stats_enabled();
    let t_total = Instant::now();

    let t_phred = Instant::now();
    let phred_offset = match cfg.phred_offset {
        PhredOffsetConfig::Fixed(v) => v,
        PhredOffsetConfig::Auto => detect_phred_offset(&cfg.reads1)
            .with_context(|| "failed to auto-detect phred offset")?,
    };
    log_stage(stats, "engine.phred_detect", t_phred);

    let encoding = if phred_offset == 64 {
        Encoding::Illumina15
    } else {
        Encoding::Sanger
    };

    let file_name = cfg
        .reads1
        .file_name()
        .and_then(|s| s.to_str())
        .map(str::to_string)
        .context("failed to determine input filename")?;

    let ctx = FinalizeContext {
        encoding,
        file_name,
        sample_name: cfg.sample_name.clone(),
        mode: cfg.mode,
    };

    let threads = cfg.threads.max(1);
    let (chunk_tx, chunk_rx) = channel::bounded::<PackedChunk>(threads * 2);
    let (result_tx, result_rx) = channel::bounded::<(Agg, WorkerStats)>(threads);
    let (err_tx, err_rx) = channel::bounded::<anyhow::Error>(1);
    let (prod_stats_tx, prod_stats_rx) = channel::bounded::<ProducerStats>(1);

    let producer_path = cfg.reads1.clone();
    let producer_err = err_tx.clone();
    let t_producer = Instant::now();
    let producer = thread::spawn(move || {
        run_producer(producer_path, chunk_tx, producer_err, prod_stats_tx);
    });
    log_stage(stats, "engine.spawn_producer", t_producer);

    let mode = cfg.mode;
    let stats_arc = Arc::new(stats);
    let mut workers = Vec::with_capacity(threads);
    let t_workers = Instant::now();
    for _ in 0..threads {
        let rx = chunk_rx.clone();
        let tx = result_tx.clone();
        let stats_enabled = *stats_arc;
        workers.push(thread::spawn(move || {
            run_worker(rx, tx, mode, phred_offset, stats_enabled);
        }));
    }
    drop(chunk_rx);
    drop(result_tx);
    drop(err_tx);
    log_stage(stats, "engine.spawn_workers", t_workers);

    let t_collect = Instant::now();
    let mut final_agg = Agg::new(cfg.mode);
    let mut worker_stats_total = WorkerStats::default();
    let mut wait_time = Duration::ZERO;
    let mut merge_time = Duration::ZERO;

    let mut err_open = true;
    let mut results_open = true;
    while results_open {
        let t_wait = Instant::now();
        let msg = if err_open {
            channel::select! {
                recv(err_rx) -> err => match err {
                    Ok(err) => return Err(err),
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
    log_stage(stats, "engine.collect_merge", t_collect);

    let _ = producer.join();
    for w in workers {
        let _ = w.join();
    }

    let prod_stats = prod_stats_rx.recv().unwrap_or_default();

    if stats {
        emit_stats(&prod_stats, &worker_stats_total, wait_time, merge_time);
    }
    log_stage(stats, "engine.total", t_total);

    Ok(RunOutput {
        agg: final_agg,
        ctx,
    })
}

fn run_producer(
    path: PathBuf,
    chunk_tx: channel::Sender<PackedChunk>,
    err_tx: channel::Sender<anyhow::Error>,
    stats_tx: channel::Sender<ProducerStats>,
) {
    let mut reader = match FastqReader::from_path_auto(&path) {
        Ok(reader) => reader,
        Err(e) => {
            let _ = err_tx.send(anyhow!("failed to open FASTQ input: {e:?}"));
            return;
        }
    };

    let mut stats = ProducerStats::default();
    let mut current = PackedChunk {
        bytes: Vec::with_capacity(TARGET_CHUNK_BYTES + 4096),
        offsets: Vec::with_capacity(TARGET_CHUNK_BYTES / 128),
    };

    loop {
        let t_next = Instant::now();
        let rec = match reader.next() {
            Ok(Some(rec)) => rec,
            Ok(None) => break,
            Err(e) => {
                let _ = err_tx.send(anyhow!("FASTQ parse/read error: {e:?}"));
                return;
            }
        };
        stats.parse += t_next.elapsed();

        let seq = rec.seq();
        let qual = rec.qual();
        if seq.len() != qual.len() {
            let _ = err_tx.send(anyhow!(
                "seq/qual length mismatch: seq={} qual={}",
                seq.len(),
                qual.len()
            ));
            return;
        }
        let len = seq.len();
        let seq_start = current.bytes.len() as u32;
        current.bytes.extend_from_slice(seq);
        current.bytes.extend_from_slice(qual);
        current.offsets.push(RecordOffset {
            seq_start,
            len: len as u32,
        });

        if current.bytes.len() >= TARGET_CHUNK_BYTES {
            let read_count = current.offsets.len() as u64;
            let chunk_bytes = current.bytes.len() as u64;
            let chunk = std::mem::replace(
                &mut current,
                PackedChunk {
                    bytes: Vec::with_capacity(TARGET_CHUNK_BYTES + 4096),
                    offsets: Vec::with_capacity(TARGET_CHUNK_BYTES / 128),
                },
            );
            if chunk_tx.send(chunk).is_err() {
                return;
            }
            stats.chunks += 1;
            stats.reads += read_count;
            stats.bytes += chunk_bytes;
        }
    }

    if !current.offsets.is_empty() {
        stats.chunks += 1;
        stats.reads += current.offsets.len() as u64;
        stats.bytes += current.bytes.len() as u64;
        let _ = chunk_tx.send(current);
    }

    let _ = stats_tx.send(stats);
}

fn run_worker(
    rx: channel::Receiver<PackedChunk>,
    tx: channel::Sender<(Agg, WorkerStats)>,
    mode: Mode,
    phred_offset: u8,
    stats_enabled: bool,
) {
    let mut agg = Agg::new(mode);
    let mut wstats = WorkerStats::default();

    for chunk in rx.iter() {
        let t_start = Instant::now();
        wstats.chunks += 1;
        wstats.bytes += chunk.bytes.len() as u64;
        wstats.reads += chunk.offsets.len() as u64;

        if stats_enabled {
            let mut timing = UpdateTimings::default();
            for off in &chunk.offsets {
                let view = view_from_offset(&chunk, *off);
                agg.update_read_timed(&view, phred_offset, &mut timing);
            }
            wstats.kmer += timing.kmer;
            wstats.kmer_updates += timing.kmer_updates;
        } else {
            for off in &chunk.offsets {
                let view = view_from_offset(&chunk, *off);
                agg.update_read(&view, phred_offset);
            }
        }

        wstats.process += t_start.elapsed();
    }

    let _ = tx.send((agg, wstats));
}

#[inline]
fn view_from_offset(chunk: &PackedChunk, off: RecordOffset) -> ReadView<'_> {
    let start = off.seq_start as usize;
    let len = off.len as usize;
    let seq = &chunk.bytes[start..start + len];
    let qual = &chunk.bytes[start + len..start + 2 * len];
    ReadView { seq, qual }
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

fn stats_enabled() -> bool {
    matches!(std::env::var("KIRA_STATS").as_deref(), Ok("1"))
}

fn log_stage(stats: bool, name: &str, t: Instant) {
    if stats {
        eprintln!("KIRA_STATS stage={} time={}", name, fmt_dur(t.elapsed()));
    }
}

fn fmt_dur(d: Duration) -> String {
    if d.as_secs_f64() < 1.0 {
        format!("{}ms", d.as_millis())
    } else {
        format!("{:.3}s", d.as_secs_f64())
    }
}

fn detect_phred_offset(path: &PathBuf) -> Result<u8> {
    let mut reader = FastqReader::from_path_auto(path)
        .map_err(|e| anyhow!("failed to open FASTQ for phred detection: {e:?}"))?;

    let mut reads: usize = 0;
    let mut min_q: u8 = u8::MAX;
    let mut max_q: u8 = 0;

    while reads < AUTO_DETECT_READS {
        let rec = match reader.next() {
            Ok(Some(rec)) => rec,
            Ok(None) => break,
            Err(e) => return Err(anyhow!("FASTQ parse/read error during phred detect: {e:?}")),
        };

        for &b in rec.qual() {
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
        bail!("input file is empty");
    }

    let offset = if min_q < 59 {
        33
    } else if min_q >= 64 {
        64
    } else if max_q <= 74 {
        33
    } else {
        64
    };

    Ok(offset)
}
