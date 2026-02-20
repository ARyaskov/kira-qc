use crate::core::fastq::{self, OwnedRead};
use crate::core::metrics::{Agg, UpdateTimings};
use crate::core::model::{Encoding, FinalizeContext, Mode};
use anyhow::{Context, Result, anyhow, bail};
use crossbeam_channel as channel;
use kira_fastq::FastqReader;
use std::path::PathBuf;
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
    pub out_dir: PathBuf,
    pub sample_name: String,
    pub threads: usize,
    pub phred_offset: PhredOffsetConfig,
    pub mode: Mode,
}

pub struct RunOutput {
    pub agg: Agg,
    pub ctx: FinalizeContext,
}

struct WorkChunk {
    index: usize,
    reads: Vec<OwnedRead>,
    bytes: usize,
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
    parse: Duration,
    metrics_core: Duration,
    adapters: Duration,
    heavyhitters: Duration,
    kmer: Duration,
    kmer_encode: Duration,
    kmer_keygen: Duration,
    kmer_binning: Duration,
    kmer_cms: Duration,
    kmer_hh: Duration,
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
        .map(|s| s.to_string())
        .context("failed to determine input filename")?;

    let ctx = FinalizeContext {
        phred_offset,
        encoding,
        file_name,
        sample_name: cfg.sample_name.clone(),
        mode: cfg.mode,
    };

    let (chunk_tx, chunk_rx) = channel::bounded::<WorkChunk>(cfg.threads * 2);
    let (result_tx, result_rx) = channel::unbounded::<(usize, Agg)>();
    let (total_tx, total_rx) = channel::bounded::<usize>(1);
    let (err_tx, err_rx) = channel::bounded::<anyhow::Error>(1);
    let (prod_stats_tx, prod_stats_rx) = channel::bounded::<ProducerStats>(1);
    let (worker_stats_tx, worker_stats_rx) = channel::unbounded::<WorkerStats>();

    let producer_path = cfg.reads1.clone();
    let producer_err = err_tx.clone();
    let t_producer = Instant::now();
    let producer = thread::spawn(move || {
        let mut reader = match FastqReader::from_path_auto(&producer_path) {
            Ok(reader) => reader,
            Err(e) => {
                let _ = producer_err.send(anyhow!("failed to open FASTQ input: {e:?}"));
                return;
            }
        };

        let mut stats = ProducerStats::default();
        let mut chunk_index = 0usize;
        let mut batch_reads = Vec::new();
        let mut batch_bytes = 0usize;

        loop {
            let t_next = Instant::now();
            let rec = match reader.next() {
                Ok(Some(rec)) => rec,
                Ok(None) => break,
                Err(e) => {
                    let _ = producer_err.send(anyhow!("FASTQ parse/read error: {e:?}"));
                    return;
                }
            };
            stats.parse += t_next.elapsed();

            let owned = fastq::OwnedRead::from_record(rec);
            batch_bytes += owned.byte_len();
            batch_reads.push(owned);

            if batch_bytes >= TARGET_CHUNK_BYTES {
                let read_count = batch_reads.len() as u64;
                let chunk_bytes = batch_bytes as u64;
                let chunk = WorkChunk {
                    index: chunk_index,
                    reads: std::mem::take(&mut batch_reads),
                    bytes: batch_bytes,
                };
                if chunk_tx.send(chunk).is_err() {
                    return;
                }
                stats.chunks += 1;
                stats.reads += read_count;
                stats.bytes += chunk_bytes;
                batch_bytes = 0;
                chunk_index += 1;
            }
        }

        if !batch_reads.is_empty() {
            let read_count = batch_reads.len() as u64;
            let chunk_bytes = batch_bytes as u64;
            let chunk = WorkChunk {
                index: chunk_index,
                reads: batch_reads,
                bytes: batch_bytes,
            };
            if chunk_tx.send(chunk).is_err() {
                return;
            }
            stats.chunks += 1;
            stats.reads += read_count;
            stats.bytes += chunk_bytes;
            chunk_index += 1;
        }

        let _ = total_tx.send(chunk_index);
        let _ = prod_stats_tx.send(stats);
    });
    log_stage(stats, "engine.spawn_producer", t_producer);

    let mut workers = Vec::with_capacity(cfg.threads);
    let t_workers = Instant::now();
    for _ in 0..cfg.threads {
        let rx = chunk_rx.clone();
        let tx = result_tx.clone();
        let stats_enabled = stats;
        let stats_tx = worker_stats_tx.clone();
        let mode = cfg.mode;
        workers.push(thread::spawn(move || {
            let mut wstats = WorkerStats::default();
            for chunk in rx.iter() {
                let mut agg = Agg::new(mode);
                let t_parse = Instant::now();
                for read in &chunk.reads {
                    let read_view = read.as_view();
                    if stats_enabled {
                        let mut ut = UpdateTimings::default();
                        agg.update_read_timed(&read_view, phred_offset, &mut ut);
                        wstats.metrics_core += ut.metrics_core;
                        wstats.adapters += ut.adapters;
                        wstats.heavyhitters += ut.heavyhitters;
                        wstats.kmer += ut.kmer;
                        wstats.kmer_encode += ut.kmer_encode;
                        wstats.kmer_keygen += ut.kmer_keygen;
                        wstats.kmer_binning += ut.kmer_binning;
                        wstats.kmer_cms += ut.kmer_cms;
                        wstats.kmer_hh += ut.kmer_hh;
                        wstats.kmer_updates += ut.kmer_updates;
                    } else {
                        agg.update_read(&read_view, phred_offset);
                    }
                }
                wstats.parse += t_parse.elapsed();
                wstats.chunks += 1;
                wstats.bytes += chunk.bytes as u64;
                wstats.reads += chunk.reads.len() as u64;

                if tx.send((chunk.index, agg)).is_err() {
                    break;
                }
            }

            if stats_enabled {
                let _ = stats_tx.send(wstats);
            }
        }));
    }
    log_stage(stats, "engine.spawn_workers", t_workers);
    drop(result_tx);
    drop(err_tx);
    drop(worker_stats_tx);

    let t_collect = Instant::now();
    let total_chunks = total_rx.recv().context("failed to receive chunk count")?;
    if total_chunks == 0 {
        return Err(anyhow!("input file is empty"));
    }

    let mut parts: Vec<Option<Agg>> = vec![None; total_chunks];
    let mut wait_time = Duration::ZERO;
    let mut err_open = true;
    for _ in 0..total_chunks {
        if err_open {
            let t_wait = Instant::now();
            channel::select! {
                recv(err_rx) -> err => {
                    match err {
                        Ok(err) => return Err(err),
                        Err(_) => {
                            err_open = false;
                            continue;
                        }
                    }
                }
                recv(result_rx) -> msg => {
                    wait_time += t_wait.elapsed();
                    let (index, agg) = msg.context("failed to receive chunk result")?;
                    if index >= parts.len() {
                        return Err(anyhow!("invalid chunk index {}", index));
                    }
                    parts[index] = Some(agg);
                }
            }
        } else {
            let t_wait = Instant::now();
            let (index, agg) = result_rx.recv().context("failed to receive chunk result")?;
            wait_time += t_wait.elapsed();
            if index >= parts.len() {
                return Err(anyhow!("invalid chunk index {}", index));
            }
            parts[index] = Some(agg);
        }
    }

    let mut final_agg = Agg::new(cfg.mode);
    let t_merge = Instant::now();
    for part in parts.into_iter().flatten() {
        final_agg.merge(&part);
    }
    let merge_time = t_merge.elapsed();
    log_stage(stats, "engine.merge", t_collect);

    let _ = producer.join();
    for worker in workers {
        let _ = worker.join();
    }

    let prod_stats = prod_stats_rx.recv().unwrap_or_default();
    let mut worker_stats = WorkerStats::default();
    for ws in worker_stats_rx.iter() {
        worker_stats.chunks += ws.chunks;
        worker_stats.bytes += ws.bytes;
        worker_stats.reads += ws.reads;
        worker_stats.parse += ws.parse;
        worker_stats.metrics_core += ws.metrics_core;
        worker_stats.adapters += ws.adapters;
        worker_stats.heavyhitters += ws.heavyhitters;
        worker_stats.kmer += ws.kmer;
        worker_stats.kmer_encode += ws.kmer_encode;
        worker_stats.kmer_keygen += ws.kmer_keygen;
        worker_stats.kmer_binning += ws.kmer_binning;
        worker_stats.kmer_cms += ws.kmer_cms;
        worker_stats.kmer_hh += ws.kmer_hh;
        worker_stats.kmer_updates += ws.kmer_updates;
    }

    if stats {
        if prod_stats.chunks > 0 {
            let avg = prod_stats.bytes as f64 / prod_stats.chunks as f64;
            eprintln!(
                "KIRA_STATS producer.chunks={} producer.avg_chunk_bytes={:.0} producer.bytes={} producer.reads={}",
                prod_stats.chunks, avg, prod_stats.bytes, prod_stats.reads
            );
        }
        eprintln!(
            "KIRA_STATS worker.chunks={} worker.bytes={} worker.reads={}",
            worker_stats.chunks, worker_stats.bytes, worker_stats.reads
        );
        eprintln!(
            "KIRA_STATS producer.fastq_read_parse={}",
            fmt_dur(prod_stats.parse)
        );
        let worker_total = worker_stats.parse
            + worker_stats.metrics_core
            + worker_stats.adapters
            + worker_stats.heavyhitters
            + worker_stats.kmer;
        eprintln!(
            "KIRA_STATS worker.parse={} worker.metrics_core={} worker.adapters={} worker.heavyhitters={} worker.kmer={} worker.total={}",
            fmt_dur(worker_stats.parse),
            fmt_dur(worker_stats.metrics_core),
            fmt_dur(worker_stats.adapters),
            fmt_dur(worker_stats.heavyhitters),
            fmt_dur(worker_stats.kmer),
            fmt_dur(worker_total)
        );
        eprintln!(
            "KIRA_STATS kmer.encode={} kmer.keygen={} kmer.binning={} kmer.cms={} kmer.hh={} kmer.updates={}",
            fmt_dur(worker_stats.kmer_encode),
            fmt_dur(worker_stats.kmer_keygen),
            fmt_dur(worker_stats.kmer_binning),
            fmt_dur(worker_stats.kmer_cms),
            fmt_dur(worker_stats.kmer_hh),
            worker_stats.kmer_updates
        );
        eprintln!(
            "KIRA_STATS reducer.wait={} reducer.merge_cost={}",
            fmt_dur(wait_time),
            fmt_dur(merge_time)
        );
    }

    log_stage(stats, "engine.total", t_total);

    Ok(RunOutput {
        agg: final_agg,
        ctx,
    })
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
            min_q = min_q.min(b);
            max_q = max_q.max(b);
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
