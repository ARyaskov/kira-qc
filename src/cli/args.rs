use clap::{Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "kira-qc", version, about = "FastQC-style QC for plain FASTQ")]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    Run(RunArgs),
}

#[derive(Parser)]
pub struct RunArgs {
    pub reads1: PathBuf,

    #[arg(long)]
    pub out: PathBuf,

    #[arg(long, default_value_t = num_cpus::get())]
    pub threads: usize,

    #[arg(long)]
    pub sample_name: Option<String>,

    #[arg(long, value_enum, default_value_t = PhredOffsetArg::Auto)]
    pub phred_offset: PhredOffsetArg,

    #[arg(long, default_value_t = false)]
    pub no_zip: bool,

    #[arg(long, value_enum, default_value_t = ModeArg::Short)]
    pub mode: ModeArg,

    #[arg(long, value_enum)]
    pub export_latex: Option<LatexExportArg>,
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

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum LatexExportArg {
    #[value(name = "summary")]
    Summary,
    #[value(name = "supplement")]
    Supplement,
}
