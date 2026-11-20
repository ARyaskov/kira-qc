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

    fs::create_dir_all(&out_root)?;
    let staging = tempfile::Builder::new()
        .prefix(".kira-qc-")
        .tempdir_in(&out_root)?;
    let report_name = format!("{sample_name}_fastqc");
    let out_dir = staging.path().join(&report_name);
    fs::create_dir(&out_dir)?;

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
        report::zip::write_zip_from(staging.path(), &out_dir, &sample_name, &extra_files)
            .with_context(|| "failed to create zip output")?;
        log.stage("zip", t_zip);
    }

    let html_name = format!("{report_name}.html");
    fs::copy(
        out_dir.join("fastqc_report.html"),
        staging.path().join(&html_name),
    )?;
    publish_report(staging, &out_root, &report_name, &html_name, args.noextract)?;
    if !args.noextract && !args.quiet {
        eprintln!("kira-qc: wrote {}", out_root.join(&report_name).display());
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

fn publish_report(
    staging: tempfile::TempDir,
    root: &Path,
    report_name: &str,
    html_name: &str,
    noextract: bool,
) -> Result<()> {
    use std::io::{BufRead, BufReader, ErrorKind};

    // Keep the lock inode stable so simultaneous publishers always lock the same file.
    let lock = fs::OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(false)
        .open(root.join(".kira-qc.lock"))?;
    lock.lock()?;

    let existing = root.join(report_name);
    for directory in [
        &existing,
        &existing.join("latex"),
        &existing.join("latex/tables"),
        &existing.join("latex/figures"),
    ] {
        match fs::symlink_metadata(directory) {
            Ok(meta) if meta.is_dir() => {}
            Ok(_) => bail!("expected a report directory: {}", directory.display()),
            Err(error) if error.kind() == ErrorKind::NotFound => {}
            Err(error) => return Err(error.into()),
        }
    }
    let staged = staging.path().join(report_name);
    let mut paths = vec![
        PathBuf::from(html_name),
        PathBuf::from(format!("{report_name}.zip")),
    ];
    let files = [
        "fastqc_data.txt",
        "summary.txt",
        "fastqc_report.html",
        "kira_qc.json",
    ];
    let is_kira_report = match fs::File::open(existing.join(files[0])) {
        Ok(file) => {
            let mut lines = BufReader::new(file).lines();
            lines
                .nth(1)
                .transpose()?
                .is_some_and(|line| line.starts_with("##kira-qc\t"))
        }
        Err(error) if error.kind() == ErrorKind::NotFound => false,
        Err(error) => return Err(error.into()),
    };
    if noextract && !is_kira_report && existing.join(files[0]).exists() {
        bail!(
            "cannot publish ZIP beside a foreign report: {}",
            existing.display()
        );
    }

    if !noextract || is_kira_report {
        paths.extend(files.map(|name| Path::new(report_name).join(name)));
        paths.extend(
            managed_latex_files()
                .into_iter()
                .map(|name| Path::new(report_name).join(name)),
        );
    }
    let backup = staging.path().join("backup");
    fs::create_dir(&backup)?;
    let mut retired = Vec::new();
    let mut installed = Vec::new();
    let result = (|| -> Result<()> {
        for path in &paths {
            let dest = root.join(path);
            match fs::symlink_metadata(&dest) {
                Ok(meta) if meta.file_type().is_file() => {
                    let saved = backup.join(path);
                    fs::create_dir_all(saved.parent().unwrap())?;
                    fs::rename(&dest, &saved)?;
                    retired.push((saved, dest));
                }
                Ok(_) => bail!("expected a regular report file: {}", dest.display()),
                Err(error) if error.kind() == ErrorKind::NotFound => {}
                Err(error) => return Err(error.into()),
            }
        }
        for path in &paths {
            let source = staging.path().join(path);
            if (noextract && source.starts_with(&staged)) || !source.is_file() {
                continue;
            }
            let dest = root.join(path);
            fs::create_dir_all(dest.parent().unwrap())?;
            fs::rename(&source, &dest)?;
            installed.push((dest, source));
        }
        Ok(())
    })();
    if let Err(error) = result {
        let mut rollback_errors = Vec::new();
        for (from, to) in installed.iter().rev().chain(retired.iter().rev()) {
            if let Err(error) = fs::rename(from, to) {
                rollback_errors.push(format!("{}: {error}", to.display()));
            }
        }
        if !rollback_errors.is_empty() {
            let retained = staging.keep();
            bail!(
                "failed to publish reports: {error:#}; recovery files in {}: {}",
                retained.display(),
                rollback_errors.join(", ")
            );
        }
        return Err(error).context("failed to publish reports");
    }
    staging
        .close()
        .context("failed to remove report staging directory")
}

fn managed_latex_files() -> Vec<PathBuf> {
    let mut paths = [
        "latex/kira_qc.tex",
        "latex/README.md",
        "latex/kira_qc_latex.zip",
        "latex/tables/basic_statistics.tex",
    ]
    .map(PathBuf::from)
    .to_vec();
    for name in [
        "per_base_quality",
        "per_sequence_quality",
        "duplication_levels",
        "adapter_content",
        "sequence_length_distribution",
        "per_tile_quality",
        "per_base_content",
        "per_sequence_gc",
        "per_base_n",
        "overrepresented_sequences",
        "kmer_content",
        "per_sequence_n",
    ] {
        paths.push(PathBuf::from(format!("latex/figures/{name}.svg")));
        paths.push(PathBuf::from(format!("latex/figures/{name}.pdf")));
    }
    paths
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
