use crate::cli::args::{Cli, FormatArg, ModeArg, PhredOffsetArg, RunArgs};
use crate::core::engine::{self, Input, PhredOffsetConfig, RunConfig};
use crate::core::metrics::{AdapterSet, ContaminantList};
use crate::core::model::{Encoding, Mode};
use crate::core::stats::StatsLog;
use crate::report;
use anyhow::{Context, Result, bail};
use clap::Parser;
use std::collections::HashSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

pub fn entry() -> Result<()> {
    let cli = Cli::parse();
    run(cli.into_run_args())
}

fn run(args: RunArgs) -> Result<()> {
    let destinations = preflight(&args)?;

    let adapter_mode = match args.mode {
        ModeArg::Short => Mode::Short,
        ModeArg::Long => Mode::Long,
    };
    let adapters = match &args.adapters {
        Some(path) => AdapterSet::from_file(path, adapter_mode == Mode::Long)?,
        None => AdapterSet::default_for(adapter_mode),
    };
    let contaminants = match &args.contaminants {
        Some(path) => ContaminantList::from_file(path)?,
        None => ContaminantList::builtin(),
    };
    let ignored_modules = match &args.limits {
        Some(path) => load_limits(path, args.quiet)?,
        None => Vec::new(),
    };

    for (reads, destination) in args.reads.iter().zip(destinations) {
        run_one(
            reads,
            &args,
            destination,
            adapter_mode,
            adapters.clone(),
            contaminants.clone(),
            &ignored_modules,
        )?;
    }
    Ok(())
}

fn preflight(args: &RunArgs) -> Result<Vec<OutputDestination>> {
    if args.reads.is_empty() {
        bail!("no input files given; run `kira-qc --help` for usage");
    }
    if args.threads == 0 {
        bail!("--threads must be >= 1");
    }
    if args.sample_name.is_some() && args.reads.len() > 1 {
        bail!("--sample-name cannot be used with more than one input file");
    }
    if let Some(name) = &args.sample_name
        && (name.trim().is_empty()
            || matches!(name.as_str(), "." | "..")
            || name
                .chars()
                .any(|c| c.is_control() || "<>:\"/\\|?*".contains(c)))
    {
        bail!("--sample-name must be a single valid filename component");
    }
    if args.no_zip && args.noextract {
        bail!("--no-zip and --noextract together would leave no output");
    }
    if args.format != FormatArg::Fastq {
        bail!(
            "input format {:?} is not supported; convert it first, for example \
             `samtools fastq in.bam | kira-qc - -o qc/`",
            args.format
        );
    }

    let mut seen_stdin = false;
    let mut seen_outputs = HashSet::new();
    let mut destinations = Vec::with_capacity(args.reads.len());
    for reads in &args.reads {
        let is_stdin = reads.as_os_str() == "-";
        if is_stdin {
            if seen_stdin {
                bail!("stdin can only be read once");
            }
            seen_stdin = true;
        } else if !reads.is_file() {
            bail!("input file not found: {}", reads.display());
        }
        let sample = match &args.sample_name {
            Some(name) => name.clone(),
            None if is_stdin => "stdin".to_string(),
            None => sample_name_from_path(reads)
                .context("failed to determine sample name from input file")?,
        };
        let root = match &args.out {
            Some(dir) => dir.as_path(),
            None if is_stdin => Path::new("."),
            None => reads
                .parent()
                .filter(|parent| !parent.as_os_str().is_empty())
                .unwrap_or(Path::new(".")),
        };
        let root = resolve_output_root(root)?;
        let path = root.join(format!("{sample}_fastqc"));
        #[cfg(windows)]
        let key = path.to_string_lossy().to_lowercase();
        #[cfg(not(windows))]
        let key = path.clone();
        if !seen_outputs.insert(key) {
            bail!(
                "multiple input files target the same report: {}",
                path.display()
            );
        }
        destinations.push(OutputDestination { root, sample });
    }
    Ok(destinations)
}

struct OutputDestination {
    root: PathBuf,
    sample: String,
}

/// Resolve existing directories (including aliases) without creating output.
fn resolve_output_root(path: &Path) -> Result<PathBuf> {
    match fs::canonicalize(path) {
        Ok(path) => Ok(path),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
            let absolute = std::path::absolute(path)?;
            let (Some(parent), Some(name)) = (absolute.parent(), absolute.file_name()) else {
                return Err(error).with_context(|| format!("cannot resolve {}", path.display()));
            };
            Ok(resolve_output_root(parent)?.join(name))
        }
        Err(error) => Err(error).with_context(|| format!("cannot resolve {}", path.display())),
    }
}

fn run_one(
    reads: &Path,
    args: &RunArgs,
    destination: OutputDestination,
    mode: Mode,
    adapters: Arc<AdapterSet>,
    contaminants: Arc<ContaminantList>,
    ignored_modules: &[String],
) -> Result<()> {
    let log = StatsLog::from_env();
    let t0 = Instant::now();
    let is_stdin = reads.as_os_str() == "-";

    let input = if is_stdin {
        Input::Stdin
    } else {
        Input::Path(reads.to_path_buf())
    };

    let OutputDestination {
        root: out_root,
        sample: sample_name,
    } = destination;
    let phred_offset = match args.phred_offset {
        PhredOffsetArg::Auto => PhredOffsetConfig::Auto,
        PhredOffsetArg::P33 => PhredOffsetConfig::Fixed(33),
        PhredOffsetArg::P64 => PhredOffsetConfig::Fixed(64),
    };

    let config = RunConfig {
        input,
        sample_name: sample_name.clone(),
        threads: args.threads,
        phred_offset,
        mode,
        kmer_enabled: args.kmers,
        dup_window: args.dup_window,
        nogroup: args.nogroup,
        adapters,
        contaminants,
        ignored_modules: ignored_modules.to_vec(),
    };

    let t_engine = Instant::now();
    let output =
        engine::run(config).with_context(|| format!("failed to analyse {}", reads.display()))?;
    log.stage("engine", t_engine);

    if !args.quiet {
        warn_about_encoding(&output.ctx.encoding, args);
        warn_about_mode(&output.agg, mode);
    }

    let t_finalize = Instant::now();
    let metrics = output.agg.finalize(&output.ctx);
    log.stage("finalize", t_finalize);

    let temp_dir = if args.noextract {
        fs::create_dir_all(&out_root)?;
        Some(
            tempfile::Builder::new()
                .prefix(".kira-qc-")
                .tempdir_in(&out_root)?,
        )
    } else {
        None
    };
    let out_dir = temp_dir.as_ref().map_or_else(
        || out_root.join(format!("{sample_name}_fastqc")),
        |dir| dir.path().to_owned(),
    );
    fs::create_dir_all(&out_dir)
        .with_context(|| format!("failed to create output dir {}", out_dir.display()))?;

    report::fastqc_txt::write(&out_dir.join("fastqc_data.txt"), &metrics, &output.ctx)?;
    report::summary_txt::write(&out_dir.join("summary.txt"), &metrics, &output.ctx)?;
    report::html::write(&out_dir.join("fastqc_report.html"), &metrics, &output.ctx)?;

    let mut extra_files: Vec<String> = Vec::new();
    if args.json {
        report::json::write(&out_dir.join("kira_qc.json"), &metrics, &output.ctx)?;
        extra_files.push("kira_qc.json".to_string());
    }

    #[cfg(feature = "latex")]
    if let Some(export) = args.export_latex {
        use crate::cli::args::LatexExportArg;
        let latex_mode = match export {
            LatexExportArg::Summary => report::latex::LatexMode::Summary,
            LatexExportArg::Supplement => report::latex::LatexMode::Supplement,
        };
        let t_latex = Instant::now();
        report::latex::write(&out_dir, &metrics, &output.ctx, latex_mode)
            .with_context(|| "failed to write LaTeX export")?;
        log.stage("latex", t_latex);
        for name in report::zip::list_files_recursive(&out_dir.join("latex"), "latex")? {
            // The Overleaf archive is a packaging of the same directory; adding
            // it would duplicate every figure inside the outer bundle.
            if name.ends_with("kira_qc_latex.zip") {
                continue;
            }
            extra_files.push(name);
        }
    }

    if !args.no_zip {
        let t_zip = Instant::now();
        report::zip::write_zip_from(&out_root, &out_dir, &sample_name, &extra_files)
            .with_context(|| "failed to create zip output")?;
        log.stage("zip", t_zip);
    }

    if let Some(temp_dir) = temp_dir {
        temp_dir
            .close()
            .context("failed to remove temporary report directory")?;
    } else if !args.quiet {
        eprintln!("kira-qc: wrote {}", out_dir.display());
    }

    log.line(format_args!(
        "input={} reads={} bases={} total={}",
        reads.display(),
        output.agg.total_reads,
        output.agg.total_bases,
        crate::core::stats::fmt_dur(t0.elapsed())
    ));
    Ok(())
}

/// FastQC's output naming: strip the compression suffix, then the FASTQ
/// suffix. `sample.fastq.gz` becomes `sample_fastqc`, not `sample.fastq_fastqc`.
pub fn sample_name_from_path(path: &Path) -> Option<String> {
    let mut name = path.file_name()?.to_str()?.to_string();
    for suffix in [".gz", ".bz2", ".zst", ".xz", ".bgz", ".bgzf"] {
        if name.to_ascii_lowercase().ends_with(suffix) {
            name.truncate(name.len() - suffix.len());
            break;
        }
    }
    for suffix in [
        ".fastq", ".fq", ".fasta", ".fa", ".fna", ".txt", ".sam", ".bam", ".cram",
    ] {
        if name.to_ascii_lowercase().ends_with(suffix) {
            name.truncate(name.len() - suffix.len());
            break;
        }
    }
    if name.is_empty() { None } else { Some(name) }
}

/// FastQC limits file. Only `ignore` is actionable here; threshold overrides
/// are reported as unsupported instead of being silently dropped.
fn load_limits(path: &Path, quiet: bool) -> Result<Vec<String>> {
    let text = fs::read_to_string(path)
        .with_context(|| format!("failed to read limits file {}", path.display()))?;
    let mut ignored = Vec::new();
    let mut unsupported = Vec::new();

    for (lineno, line) in text.lines().enumerate() {
        let line = line.split('#').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 3 {
            bail!(
                "{}:{}: expected '<module> <warn|error|ignore> <value>'",
                path.display(),
                lineno + 1
            );
        }
        let value: f64 = fields[fields.len() - 1]
            .parse()
            .with_context(|| format!("{}:{}: value is not a number", path.display(), lineno + 1))?;
        let directive = fields[fields.len() - 2].to_ascii_lowercase();
        let module = fields[..fields.len() - 2].join(" ");

        match directive.as_str() {
            "ignore" => {
                if value != 0.0 {
                    let anchor = report::resolve_limits_key(&module).ok_or_else(|| {
                        anyhow::anyhow!(
                            "{}:{}: unknown module '{}'",
                            path.display(),
                            lineno + 1,
                            module
                        )
                    })?;
                    ignored.push(anchor.to_string());
                }
            }
            "warn" | "error" => unsupported.push(module),
            other => bail!(
                "{}:{}: unknown directive '{}'",
                path.display(),
                lineno + 1,
                other
            ),
        }
    }

    if !quiet && !unsupported.is_empty() {
        unsupported.sort();
        unsupported.dedup();
        eprintln!(
            "kira-qc: warning: threshold overrides in {} are not applied ({}); \
             only `ignore` directives are honoured",
            path.display(),
            unsupported.join(", ")
        );
    }
    Ok(ignored)
}

fn warn_about_encoding(encoding: &Encoding, args: &RunArgs) {
    if *encoding == Encoding::Illumina15 && matches!(args.phred_offset, PhredOffsetArg::Auto) {
        eprintln!(
            "kira-qc: warning: quality scores look like Phred+64 (Illumina 1.5). \
             If this is modern data, re-run with --phred-offset 33."
        );
    }
}

fn warn_about_mode(agg: &crate::core::metrics::Agg, mode: Mode) {
    let max_len = agg.length_hist.max_len();
    match mode {
        Mode::Short if max_len > 1_000 => eprintln!(
            "kira-qc: warning: reads up to {} bp in short-read mode. \
             Per-base modules assume uniform short reads; consider --mode long.",
            max_len
        ),
        Mode::Long if max_len > 0 && max_len < 500 => eprintln!(
            "kira-qc: warning: longest read is {} bp in long-read mode. \
             Consider --mode short for Illumina data.",
            max_len
        ),
        _ => {}
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sample_names_drop_compression_and_fastq_suffixes() {
        let cases = [
            ("sample.fastq.gz", "sample"),
            ("sample.fq.gz", "sample"),
            ("sample.fastq", "sample"),
            ("sample.fq", "sample"),
            ("sample.txt", "sample"),
            ("reads_R1_001.fastq.gz", "reads_R1_001"),
            ("sample", "sample"),
            ("sample.other", "sample.other"),
        ];
        for (input, expect) in cases {
            assert_eq!(
                sample_name_from_path(Path::new(input)).as_deref(),
                Some(expect),
                "input {input}"
            );
        }
    }

    #[test]
    fn sample_names_are_case_insensitive() {
        assert_eq!(
            sample_name_from_path(Path::new("Sample.FASTQ.GZ")).as_deref(),
            Some("Sample")
        );
    }
}
