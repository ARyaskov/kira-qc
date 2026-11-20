use clap::{Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

fn default_threads() -> usize {
    std::thread::available_parallelism()
        .map(usize::from)
        .unwrap_or(1)
}

#[derive(Parser)]
#[command(
    name = "kira-qc",
    version,
    about = "FastQC-compatible quality control for FASTQ files",
    long_about = "FastQC-compatible quality control for FASTQ, FASTQ.GZ and BGZF input.\n\
                  Reads can also be piped in by passing '-' as the file name.",
    args_conflicts_with_subcommands = true
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Option<Commands>,

    /// Accepts `kira-qc reads.fastq -o qc/` without the `run` subcommand, the
    /// way `fastqc` itself is invoked.
    #[command(flatten)]
    pub run: RunArgs,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Run QC on one or more FASTQ files.
    Run(RunArgs),
}

#[derive(Parser, Clone)]
pub struct RunArgs {
    /// Input FASTQ / FASTQ.GZ / BGZF files, or `-` to read from stdin.
    pub reads: Vec<PathBuf>,

    /// Output directory. Defaults to the directory holding the input file.
    #[arg(short = 'o', long = "outdir", alias = "out", value_name = "DIR")]
    pub out: Option<PathBuf>,

    /// Number of worker threads.
    #[arg(short = 't', long, default_value_t = default_threads(), value_name = "N")]
    pub threads: usize,

    /// Sample name; only valid with a single input file.
    #[arg(long, value_name = "NAME")]
    pub sample_name: Option<String>,

    /// Quality encoding: auto-detect, or force an offset.
    #[arg(long, value_enum, default_value_t = PhredOffsetArg::Auto)]
    pub phred_offset: PhredOffsetArg,

    /// Sequencing technology profile.
    #[arg(long, value_enum, default_value_t = ModeArg::Short)]
    pub mode: ModeArg,

    /// Report every base position separately instead of FastQC base groups.
    #[arg(long, default_value_t = false)]
    pub nogroup: bool,

    /// Enable the K-mer content module (a heuristic, off by default).
    #[arg(short = 'k', long = "kmers", default_value_t = false)]
    pub kmers: bool,

    /// Reads used to estimate duplication and overrepresentation; 0 means all.
    #[arg(long, default_value_t = crate::core::metrics::DEFAULT_DUP_WINDOW, value_name = "N")]
    pub dup_window: u64,

    /// Do not create the ZIP bundle.
    #[arg(long, default_value_t = false)]
    pub no_zip: bool,

    /// Keep ZIP and standalone HTML without extracting the report directory.
    #[arg(long, default_value_t = false)]
    pub noextract: bool,

    /// Keep the unpacked report directory (the default).
    #[arg(long, default_value_t = false)]
    pub extract: bool,

    /// Custom adapter list (`name<TAB>sequence`, `#` comments).
    #[arg(short = 'a', long, value_name = "FILE")]
    pub adapters: Option<PathBuf>,

    /// Custom contaminant list used to label overrepresented sequences.
    #[arg(short = 'c', long, value_name = "FILE")]
    pub contaminants: Option<PathBuf>,

    /// FastQC limits file. Only `ignore` directives are honoured.
    #[arg(short = 'l', long, value_name = "FILE")]
    pub limits: Option<PathBuf>,

    /// Input format. Only FASTQ is read natively.
    #[arg(short = 'f', long, value_enum, default_value_t = FormatArg::Fastq)]
    pub format: FormatArg,

    /// Also write a machine-readable `kira_qc.json`.
    #[arg(long, default_value_t = false)]
    pub json: bool,

    /// Generate a LaTeX export alongside the HTML report.
    #[cfg(feature = "latex")]
    #[arg(long, value_enum)]
    pub export_latex: Option<LatexExportArg>,

    /// Suppress progress and warning messages.
    #[arg(short = 'q', long, default_value_t = false)]
    pub quiet: bool,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum PhredOffsetArg {
    #[value(name = "auto")]
    Auto,
    #[value(name = "33")]
    P33,
    #[value(name = "64")]
    P64,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum ModeArg {
    #[value(name = "short")]
    Short,
    #[value(name = "long")]
    Long,
}

#[derive(Clone, Copy, Debug, ValueEnum, PartialEq, Eq)]
pub enum FormatArg {
    #[value(name = "fastq")]
    Fastq,
    #[value(name = "bam")]
    Bam,
    #[value(name = "sam")]
    Sam,
}

#[cfg(feature = "latex")]
#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum LatexExportArg {
    #[value(name = "summary")]
    Summary,
    #[value(name = "supplement")]
    Supplement,
}

impl Cli {
    /// Arguments from either invocation style.
    pub fn into_run_args(self) -> RunArgs {
        match self.command {
            Some(Commands::Run(args)) => args,
            None => self.run,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_definition_is_valid() {
        Cli::command().debug_assert();
    }

    #[test]
    fn bare_invocation_matches_fastqc_style() {
        let cli = Cli::try_parse_from(["kira-qc", "reads.fastq", "-o", "qc"]).unwrap();
        let args = cli.into_run_args();
        assert_eq!(args.reads, vec![PathBuf::from("reads.fastq")]);
        assert_eq!(args.out, Some(PathBuf::from("qc")));
    }

    #[test]
    fn run_subcommand_still_works() {
        let cli = Cli::try_parse_from(["kira-qc", "run", "reads.fastq", "--out", "qc"]).unwrap();
        let args = cli.into_run_args();
        assert_eq!(args.reads, vec![PathBuf::from("reads.fastq")]);
        assert_eq!(args.out, Some(PathBuf::from("qc")));
    }

    #[test]
    fn multiple_inputs_are_accepted() {
        let cli = Cli::try_parse_from(["kira-qc", "a.fastq", "b.fastq.gz", "-t", "4"]).unwrap();
        let args = cli.into_run_args();
        assert_eq!(args.reads.len(), 2);
        assert_eq!(args.threads, 4);
    }

    #[test]
    fn kmers_are_off_by_default() {
        let cli = Cli::try_parse_from(["kira-qc", "a.fastq"]).unwrap();
        assert!(!cli.into_run_args().kmers);
    }
}
