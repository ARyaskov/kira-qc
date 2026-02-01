use crate::core::fastq;
use crate::core::io::{InputKind, InputSource, MmapSource, detect_input_kind, open_gzip_reader};
use crate::core::metrics::{Agg, UpdateTimings};
use crate::core::model::{Encoding, FinalizeContext, Mode};
use anyhow::{Context, Result, anyhow, bail};
use crossbeam_channel as channel;
use std::io::BufReader;
use std::path::PathBuf;
use std::sync::Arc;
use std::thread;
use std::time::{Duration, Instant};

const AUTO_DETECT_READS: usize = 50_000;

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

#[derive(Clone, Debug, Default)]
struct ProducerStats {
    chunks: u64,
    bytes: u64,
    decompress: Duration,
    align: Duration,
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
    let kind = detect_input_kind(&cfg.reads1)?;
    let t_phred = Instant::now();
    let phred_offset = match cfg.phred_offset {
        PhredOffsetConfig::Fixed(v) => v,
        PhredOffsetConfig::Auto => detect_phred_offset(&cfg.reads1, kind, cfg.threads)
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

    let t_open = Instant::now();
    let (mut input, mmap_source_opt, _kind) = InputSource::open(&cfg.reads1, cfg.threads)?;
    log_stage(stats, "engine.input_open", t_open);
    let (chunk_tx, chunk_rx) = channel::bounded::<crate::core::io::Chunk>(cfg.threads * 2);
    let (result_tx, result_rx) = channel::unbounded::<(usize, Agg)>();
    let (total_tx, total_rx) = channel::bounded::<usize>(1);
    let (err_tx, err_rx) = channel::bounded::<anyhow::Error>(1);
    let (prod_stats_tx, prod_stats_rx) = channel::bounded::<ProducerStats>(1);
    let (worker_stats_tx, worker_stats_rx) = channel::unbounded::<WorkerStats>();

    let producer_err = err_tx.clone();
    let t_producer = Instant::now();
    let producer = thread::spawn(move || {
        let mut count = 0usize;
        let mut stats = ProducerStats::default();
        loop {
            match input.next_chunk() {
                Ok(Some(chunk)) => {
                    stats.chunks += 1;
                    stats.bytes += chunk.timing.bytes as u64;
                    stats.decompress += chunk.timing.decompress;
                    stats.align += chunk.timing.align;
                    if chunk_tx.send(chunk).is_err() {
                        return;
                    }
                    count += 1;
                }
                Ok(None) => break,
                Err(e) => {
                    let _ = producer_err.send(e);
                    return;
                }
            }
        }
        let _ = total_tx.send(count);
        let _ = prod_stats_tx.send(stats);
    });
    log_stage(stats, "engine.spawn_producer", t_producer);

    let mut workers = Vec::with_capacity(cfg.threads);
    let t_workers = Instant::now();
    for _ in 0..cfg.threads {
        let rx = chunk_rx.clone();
        let tx = result_tx.clone();
        let err = err_tx.clone();
        let worker_source = mmap_source_opt.as_ref().map(Arc::clone);
        let stats_enabled = stats;
        let stats_tx = worker_stats_tx.clone();
        workers.push(thread::spawn(move || {
            let mut wstats = WorkerStats::default();
            for chunk in rx.iter() {
                let slice = match &chunk.data {
                    crate::core::io::ChunkData::MmapRange { start, end } => {
                        let source = match &worker_source {
                            Some(s) => s,
                            None => {
                                let _ = err
                                    .send(anyhow!("mmap source missing for chunk {}", chunk.index));
                                break;
                            }
                        };
                        &source.bytes()[*start..*end]
                    }
                    crate::core::io::ChunkData::Owned(data) => data.as_slice(),
                };
                let mut reads = Vec::new();
                let t_parse = Instant::now();
                if let Err(e) = fastq::parse_chunk(slice, &mut reads) {
                    let msg = format!(
                        "FASTQ parse error in chunk {} at offset {}",
                        chunk.index, e.byte_offset
                    );
                    let _ = err.send(anyhow!(msg));
                    break;
                }
                wstats.parse += t_parse.elapsed();
                wstats.chunks += 1;
                wstats.bytes += chunk.timing.bytes as u64;
                wstats.reads += reads.len() as u64;
                let mut agg = Agg::new(cfg.mode);
                if stats_enabled {
                    let mut ut = UpdateTimings::default();
                    for read in &reads {
                        agg.update_read_timed(read, phred_offset, &mut ut);
                    }
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
                    for read in &reads {
                        agg.update_read(read, phred_offset);
                    }
                }
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
                "KIRA_STATS producer.chunks={} producer.avg_chunk_bytes={:.0} producer.bytes={}",
                prod_stats.chunks, avg, prod_stats.bytes
            );
        }
        eprintln!(
            "KIRA_STATS worker.chunks={} worker.bytes={} worker.reads={}",
            worker_stats.chunks, worker_stats.bytes, worker_stats.reads
        );
        let producer_total = prod_stats.decompress + prod_stats.align;
        eprintln!(
            "KIRA_STATS producer.decompress={} producer.chunk_align={} producer.total={}",
            fmt_dur(prod_stats.decompress),
            fmt_dur(prod_stats.align),
            fmt_dur(producer_total)
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

fn detect_phred_offset(path: &PathBuf, kind: InputKind, threads: usize) -> Result<u8> {
    let (min_q, max_q, _) = match kind {
        InputKind::Plain => {
            let source = MmapSource::open(path)?;
            if source.len() == 0 {
                bail!("input file is empty");
            }
            fastq::scan_qual_range(source.bytes(), AUTO_DETECT_READS).map_err(|e| anyhow!(e))?
        }
        InputKind::Gzip => {
            let reader = open_gzip_reader(path, threads)?;
            let reader = BufReader::new(reader);
            fastq::scan_qual_range_reader(reader, AUTO_DETECT_READS).map_err(|e| anyhow!(e))?
        }
    };
    // Heuristic: phred33 typically has low ASCII (<59); phred64 clusters higher.
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
