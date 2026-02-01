use crate::cli::args::{Cli, Commands, LatexExportArg, ModeArg, PhredOffsetArg, RunArgs};
use crate::core::engine::{self, PhredOffsetConfig, RunConfig};
use crate::core::model::Mode;
use crate::report;
use anyhow::{Context, Result, bail};
use clap::Parser;
use std::env;
use std::fs;
use std::time::{Duration, Instant};

pub fn entry() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Run(args) => run(args),
    }
}

fn run(args: RunArgs) -> Result<()> {
    let stats = stats_enabled();
    let t0 = Instant::now();

    stage(stats, "preflight", || {
        if args.reads1.as_os_str() == "-" {
            bail!("stdin is not supported in Stage 1; provide a FASTQ file path");
        }
        if !args.reads1.is_file() {
            bail!("input file not found: {}", args.reads1.display());
        }
        if args.threads == 0 {
            bail!("--threads must be >= 1");
        }
        Ok(())
    })?;

    let input_size = fs::metadata(&args.reads1).map(|m| m.len()).unwrap_or(0);

    let t_name = Instant::now();
    let sample_name = match args.sample_name {
        Some(s) => s,
        None => args
            .reads1
            .file_stem()
            .and_then(|s| s.to_str())
            .map(|s| s.to_string())
            .context("failed to determine sample name from input file")?,
    };
    stage_done(stats, "sample-name", t_name);

    let t_phred = Instant::now();
    let phred_offset = match args.phred_offset {
        PhredOffsetArg::Auto => PhredOffsetConfig::Auto,
        PhredOffsetArg::P33 => PhredOffsetConfig::Fixed(33),
        PhredOffsetArg::P64 => PhredOffsetConfig::Fixed(64),
    };
    stage_done(stats, "phred-config", t_phred);

    let t_mode = Instant::now();
    let mode = match args.mode {
        ModeArg::Short => Mode::Short,
        ModeArg::Long => Mode::Long,
    };
    stage_done(stats, "mode", t_mode);

    let t_out = Instant::now();
    let out_dir = args.out.join(format!("{}_fastqc", sample_name));
    fs::create_dir_all(&out_dir)
        .with_context(|| format!("failed to create output dir {}", out_dir.display()))?;
    stage_done(stats, "mkdir", t_out);

    let config = RunConfig {
        reads1: args.reads1.clone(),
        out_dir: out_dir.clone(),
        sample_name: sample_name.clone(),
        threads: args.threads,
        phred_offset,
        mode,
    };

    let t_engine = Instant::now();
    let output = engine::run(config)?;
    stage_done(stats, "engine", t_engine);
    if stats {
        eprintln!(
            "KIRA_STATS input={} bytes={} reads={} bases={}",
            args.reads1.display(),
            input_size,
            output.agg.total_reads,
            output.agg.total_bases
        );
    }

    let fastqc_path = out_dir.join("fastqc_data.txt");
    let summary_path = out_dir.join("summary.txt");
    let html_path = out_dir.join("fastqc_report.html");

    let t_fastqc = Instant::now();
    report::fastqc_txt::write(&fastqc_path, &output)
        .with_context(|| format!("failed to write {}", fastqc_path.display()))?;
    stage_done(stats, "fastqc_data", t_fastqc);
    if stats {
        let fastqc_size = fs::metadata(&fastqc_path).map(|m| m.len()).unwrap_or(0);
        eprintln!(
            "KIRA_STATS output fastqc_data={} bytes={}",
            fastqc_path.display(),
            fastqc_size
        );
    }

    let t_summary = Instant::now();
    report::summary_txt::write(&summary_path, &output)
        .with_context(|| format!("failed to write {}", summary_path.display()))?;
    stage_done(stats, "summary", t_summary);
    if stats {
        let summary_size = fs::metadata(&summary_path).map(|m| m.len()).unwrap_or(0);
        eprintln!(
            "KIRA_STATS output summary={} bytes={}",
            summary_path.display(),
            summary_size
        );
    }

    let t_html = Instant::now();
    report::html::write(&html_path, &output)
        .with_context(|| format!("failed to write {}", html_path.display()))?;
    stage_done(stats, "html", t_html);
    if stats {
        let html_size = fs::metadata(&html_path).map(|m| m.len()).unwrap_or(0);
        eprintln!(
            "KIRA_STATS output html={} bytes={}",
            html_path.display(),
            html_size
        );
    }

    if !args.no_zip {
        let t_zip = Instant::now();
        report::zip::write_zip(&args.out, &sample_name)
            .with_context(|| "failed to create zip output")?;
        stage_done(stats, "zip", t_zip);
        if stats {
            let zip_path = args.out.join(format!("{}_fastqc.zip", sample_name));
            let zip_size = fs::metadata(&zip_path).map(|m| m.len()).unwrap_or(0);
            eprintln!(
                "KIRA_STATS output zip={} bytes={}",
                zip_path.display(),
                zip_size
            );
        }
    }

    if let Some(export) = args.export_latex {
        let t_latex = Instant::now();
        let mode = match export {
            LatexExportArg::Summary => report::latex::LatexMode::Summary,
            LatexExportArg::Supplement => report::latex::LatexMode::Supplement,
        };
        report::latex::write(&out_dir, &output, mode)
            .with_context(|| "failed to write LaTeX export")?;
        stage_done(stats, "latex", t_latex);
    }

    if stats {
        eprintln!("KIRA_STATS output_dir={}", out_dir.display());
        eprintln!("KIRA_STATS total={}", fmt_dur(t0.elapsed()));
    }

    Ok(())
}

fn stats_enabled() -> bool {
    matches!(env::var("KIRA_STATS").as_deref(), Ok("1"))
}

fn stage<F>(stats: bool, name: &str, f: F) -> Result<()>
where
    F: FnOnce() -> Result<()>,
{
    let t = Instant::now();
    let res = f();
    if stats {
        eprintln!("KIRA_STATS stage={} time={}", name, fmt_dur(t.elapsed()));
    }
    res
}

fn stage_done(stats: bool, name: &str, t: Instant) {
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
